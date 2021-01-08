#######################################################################################
# Functions used for scRNAseq analysis of Levesque’s data
# Autors: Leonor Schubert, Jonathan Haab, Flavio Rump
# Autumn 2020
# Course STA426, UZH
#######################################################################################


#read 10x data that has empty drops filtered but nothing else(?)
read_10x <- function() {
  sample_names <- list("patient1_HS", "patient1_SCC", "patient2_HS", "patient2_AK")
  patient1_HS.path <- file.path("data", "patient1_HS")
  patient1_SCC.path <- file.path("data", "patient1_SCC")
  patient2_HS.path <- file.path("data", "patient2_HS")
  patient2_AK.path <- file.path("data", "patient2_AK")
  
  
  # read in filtered data and create list of SingleCellExperiment objects
  paths <- list(patient1_HS.path, patient1_SCC.path, patient2_HS.path, patient2_AK.path)
  sces <- lapply(paths, function(i) read10xCounts(file.path(i, "outs/filtered_feature_bc_matrix"), col.names = TRUE))
  return(sces)
}
#read previusly stored data
read_previous_data <- function() {
  sces <- list()
  sces[[1]] <- readRDS('data/patient1_HS/clean_data/patient1_HS_cleanData_sce.RDS')
  sces[[2]] <- readRDS('data/patient1_SCC/clean_data/patient1_SCC_cleanData_sce.RDS')
  sces[[3]] <- readRDS('data/patient2_HS/clean_data/patient2_HS_cleanData_sce.RDS')
  sces[[4]] <- readRDS('data/patient2_AK/clean_data/patient2_AK_cleanData_sce.RDS')
  return(sces)
}

#read raw data and conduct QC
read_and_qc <- function(cores=4)
{
  sces <- list()
  sces[[1]] <- processCellRangerOutput("patient1_HS", 100000, TRUE, TRUE, TRUE, cores)
  sces[[2]] <- processCellRangerOutput("patient1_SCC", 100000, TRUE, TRUE, TRUE, cores)
  sces[[3]] <- processCellRangerOutput("patient2_HS", 100000, TRUE, TRUE, TRUE, cores)
  sces[[4]] <- processCellRangerOutput("patient2_AK", 100000, TRUE, TRUE, TRUE, cores)
  return(sces)
}


split_data <- function(sceo) {
  # Split the data, store ADT in alternative experiment
  sceo <- splitAltExps(sceo, rowData(sceo)$Type)
  
  # Coerce sparse matrix for ADT into a dense matrix
  counts(altExp(sceo)) <- as.matrix(counts(altExp(sceo)))
  sceo
}

# get rid of seldom detected genes
remove_rare_genes <- function(sceo, num_genes) {
  sceo <- sceo[(rowSums(counts(sceo) > 0) > num_genes),]
}

variance_stabilization <- function(sceo, num_genes) {
  
  sceo <- remove_rare_genes(sceo, 4)
  vsto <- suppressWarnings(sctransform::vst(counts(sceo)))
  logcounts(sceo, withDimnames=FALSE) <- vsto$y
  
  # Check that new assay was added to sce
  #assays(sce.patient1_HS)
  
  # get highly-variable genes
  hvgo <- row.names(sceo)[order(vsto$gene_attr$residual_variance, decreasing=TRUE)[1:num_genes]]
  results <- list("sceo" <- sceo, "hvgo" <- hvgo)
  return(results)
}

sc_PCA <- function(sceo, hvgo, known_markers) {
  # only gives index of the first encountered
  known_genes <- rownames(sceo)[which(rowData(sceo)$Symbol %in% known_markers)]
  
  # Check that the genes we know are also part of the highly variable genes
  idx_notfound <- which(!(known_genes %in% hvgo))
  
  # print the markers that were not found in the data
  #known_markers[idx_notfound]
  
  # if some of the known markers were not selected, add them
  if (length(idx_notfound) > 0) {
    cat("Adding", known_markers[idx_notfound], "to the gene set used for PCA in", attr(sceo, "name"), "\n")
    hvgo <- c(hvgo, known_genes[idx_notfound])
  }
  
  # using highly variable genes
  sceo <- runPCA(sceo, subset_row=hvgo)
  
  # check the variance explained by the PCs:
  pc.var <- attr(reducedDim(sceo),"percentVar")
  plot(pc.var, xlab="PCs", ylab="% variance explained", main=paste("Variance explained across PCs for", attr(sceo, "name")))
  
  # restrict to the first 20 components:
  reducedDim(sceo) <- reducedDim(sceo)[,1:20]
  
  # run TSNE and UMAP based on the PCA:
  sceo <- runTSNE(sceo, dimred="PCA")
  sceo <- runUMAP(sceo, dimred="PCA")
}

sc_cluster <- function(sceo, k=30) {
  go <- buildKNNGraph(sceo, BNPARAM=AnnoyParam(), use.dimred="PCA", k=k)
  
  sceo$cluster <- as.factor(cluster_louvain(go)$membership)
  
  #plots <- plot_grid( plotTSNE(sceo, colour_by="cluster", text_by="cluster") + ggtitle(attr(sceo, "name")), plotUMAP(sceo, colour_by="cluster", text_by="cluster") + ggtitle(attr(sceo, "name")) ) 
  
  plots <- plot_grid( 
    plotTSNE(sceo, colour_by="cluster", text_by="cluster"), 
    plotUMAP(sceo, colour_by="cluster", text_by="cluster"))
  title <- ggdraw() + draw_label(attr(sceo, "name"), fontface='bold')
  plots <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1))
  print(plots)
  
  return(list(sceo, go))
}

convert_rownames <- function(sceo) {
  return(paste(rownames(sceo), rowData(sceo)$Symbol, sep = "."))
}

paste_rownames <- function(sceo) {
  return(paste(rowData(sceo)$ID, rowData(sceo)$Symbol, sep = "."))
}

annotate_cells <- function(sceo, go, genes) {
  kmo <- lapply(genes, FUN=function(go) grep(paste0(go, "$", collapse="|"), rownames(sceo), value=TRUE))
  
  # TODO: try to modify this so that we don't have to convert the rownames anymore and use
  # the marker names save in rowData(sceo)$Symbol directly
  #kmo <- lapply(genes, FUN=function(go) grep(paste0(go, "$", collapse="|"), rowData(sceo)$Symbol, value=TRUE))
  #print(kmo)
  
  return(list(sceo, kmo))
}

# mean logcounts by cluster
pseudobulk <- function(sceo, kmo) {
  pbo <- aggregateData(sceo, "logcounts", by=c("cluster"), fun="mean")
  
  # TODO : find a way to delete the left annotation which is unreadable
  
  # build a heatmap of the mean logcounts of the known markers:
  h <- pheatmap(assay(pbo)[unlist(kmo),], annotation_row=data.frame(row.names=unlist(kmo), type=rep(names(kmo), lengths(kmo))), split=rep(names(kmo), lengths(kmo)), annotation_names_row=F, cluster_rows=FALSE, scale="row", main=paste(attr(sceo, "name"),"before markers aggregation"), fontsize_row=6, fontsize_col=10, angle_col = "45")
  print(h)
  
  #--- aggregation markers
  # we will assign clusters to the cell type whose markers are the most expressed
  
  # we extract the pseudo-bulk counts of the markers for each cluster
  mato <- assay(pbo)[unlist(kmo),]
  
  # we aggregate across markers of the same type
  mato <- aggregate(t(scale(t(mato))), by=list(type=rep(names(kmo), lengths(kmo))), FUN=sum)
  
  # for each column (cluster), we select the row (cell type) which has the maximum aggregated value
  cl2o<- mato[,1][apply(mato[,-1], 2, FUN=which.max)]
  # we convert the cells' cluster labels to cell type labels:
  sceo$cluster2 <- cl2o[sceo$cluster]
  
  # we aggregate again to pseudo-bulk using the new clusters
  pbo <- aggregateData(sceo, "logcounts", by=c("cluster2"), fun="mean")
  
  # we plot again the expression of the markers as a sanity check
  
  # If we want to hide the gene markers show_rownames = FALSE
  h1 <- pheatmap(assay(pbo)[unlist(kmo),], annotation_row=data.frame(row.names=unlist(kmo), type=rep(names(kmo), lengths(kmo))), split=rep(names(kmo), lengths(kmo)), annotation_names_row=F, cluster_rows=FALSE, scale="row", main=paste(attr(sceo, "name"),"after markers aggregation"), fontsize_row=6, fontsize_col=10, angle_col = "45")
  print(h1)
  
  # UMAP plot
  p <- plotUMAP(sceo, colour_by="cluster2", text_by="cluster2", point_size=1) + ggtitle(attr(sceo, "name")) + theme(legend.title=element_blank())
  print(p)
  
  return(list(sceo, pbo, mato, cl2o))
}

cluster_subtype <- function(sceo, cell_type, known_markers, subtype_markers) {
  # select the cell labeled as 'cell type', i.e. keratinocytes in the previous step
  sceo.sub <- sceo[,sceo$cluster2==cell_type]
  
  sceo.sub <- remove_rare_genes(sceo.sub, 4)
  results <- variance_stabilization(sceo.sub, 2000)
  
  sceo.sub <- results[[1]]
  hvgo.sub <- results[[2]]
  
  #---- run the pipeline again on that subdataset
  sceo.sub <- sc_PCA(sceo.sub, hvgo.sub, known_markers)
  
  #--- clustering
  
  results <- sc_cluster(sceo.sub)
  sceo.sub <- results[[1]]
  go.sub <- results[[2]]
  
  results <- annotate_cells(sceo.sub, go.sub, subtype_markers)
  
  sceo.sub <- results[[1]]
  kmo.sub <- results[[2]]
  
  #---
  
  results <- pseudobulk(sceo.sub, kmo.sub)
  
  sceo.sub <- results[[1]]
  pbo.sub <- results[[2]]
  mato.sub <- results[[3]]
  cl2o.sub <- results[[4]]
  
  return(sceo.sub)
  
}

dynamic_plot <- function(sceo, go, genes) {
  kmo <- lapply(genes, FUN=function(go) grep(paste0(go, "$", collapse="|"), rownames(sceo), value=TRUE))
  
  # mean logcounts by cluster:
  pbo <- aggregateData(sceo, "logcounts", by=c("cluster"), fun="mean")
  lengths(kmo)
  h <- pheatmap(assay(pbo)[unlist(kmo),], annotation_row=data.frame(row.names=unlist(kmo), type=rep(names(kmo), lengths(kmo))), split=rep(names(kmo), lengths(kmo)), cluster_rows=FALSE, scale="row", main="Before markers aggregation", fontsize_row=6, fontsize_col=10)
  h
  
  # we will assign clusters to the cell type whose markers are the most expressed
  # we extract the pseudo-bulk counts of the markers for each cluster
  mato <- assay(pbo)[unlist(kmo),]
  
  # we aggregate across markers of the same type
  mato <- aggregate(t(scale(t(mato))), by=list(type=rep(names(kmo), lengths(kmo))), FUN=sum)
  
  # for each column (cluster), we select the row (cell type) which has the maximum aggregated value
  cl3o <- mato[,1][apply(mato[,-1], 2, FUN=which.max)]
  
  # we convert the cells' cluster labels to cell type labels:
  sceo$cluster3 <- cl3o[sceo$cluster]
  
  print(sceo)
  
  # we aggregate again to pseudo-bulk using the new clusters
  pbo <- aggregateData(sceo, "logcounts", by=c("cluster3"), fun="mean")
  # we plot again the expression of the markers as a sanity check
  h1 <- pheatmap(assay(pbo)[unlist(kmo),], annotation_row=data.frame(row.names=unlist(kmo), type=rep(names(kmo), lengths(kmo))), split=rep(names(kmo), lengths(kmo)), cluster_rows=FALSE, scale="row", main="After markers aggregation", fontsize_row=6, fontsize_col=10) #, scale="row", main="After markers aggregation", fontsize_row=6, fontsize_col=10)
  h1
  plotUMAP(sceo, colour_by="cluster3", text_by="cluster3", point_size=1)
  
  
  #---- Histograms
  total <- ncol(sceo)
  basal <- ncol(sceo[,sceo$cluster3=="keratinocyte_basal"])
  cycling <- ncol(sceo[,sceo$cluster3=="keratinocyte_cycling"])
  diff <- ncol(sceo[,sceo$cluster3=="keratinocyte_differentiating"])
  counts <- c(basal, cycling, diff)
  percentage <- counts/total
  kera.dyn <- data.frame(names=c("Basal", "Cycling", "Differentiating"), percentage=percentage)
  kera.dyn
  
  #p <- ggplot(data=kera.dyn, aes(x=names, y=percentage, fill=names)) + geom_bar(stat="identity") + labs(title=paste("Proportion of Keratinocyte states in", attr(sceo, "name")),x="Subtypes", y="Percentage") + ylim(0,1) + theme(legend.position='none')
  #print(p)
  
  return(kera.dyn)
  
}

# create data frame to make barplot of celltypes
df_barplot_celltypes <- function(sceo){
  
  total <- ncol(sceo)
  
  keratinocyte <- ncol(sceo[, sceo$cluster2=="keratinocyte"])
  fibroblast <- ncol(sceo[, sceo$cluster2=="fibroblast"])
  endothelial <- ncol(sceo[, sceo$cluster2=="endothelial"])
  myeloid <- ncol(sceo[, sceo$cluster2=="myeloid"])
  lymphocyte <- ncol(sceo[, sceo$cluster2=="lymphocyte"])
  melanocyte <- ncol(sceo[, sceo$cluster2=="melanocyte"])
  schwann <- ncol(sceo[, sceo$cluster2=="schwann"])
  mitotic <- ncol(sceo[, sceo$cluster2=="mitotic"])
  
  
  counts <- c(keratinocyte, fibroblast, endothelial, myeloid, lymphocyte, melanocyte, schwann, mitotic)
  percentage <- counts/total
  
  df.celltypes <- data.frame(names=c("keratinocyte", "fibroblast", "endothelial", "myeloid", "lymphocyte", "melanocyte", "schwann", "mitotic"), percentage=percentage)
  
  return(df.celltypes)
  
}

# create comparable barplot of cell composition in healthy vs. disease
dynamic_barplot <- function(df, name, labels, title){
  
  p <- ggplot2.barplot(data=df, xName="names", yName="percentage",
                       groupName="Type", 
                       position=position_dodge(),
                       #background and line colors
                       backgroundColor="white", color="black", 
                       xtitle="Subtypes", ytitle="Proportion", 
                       mainTitle=paste(title, name),
                       removePanelGrid=TRUE, removePanelBorder=TRUE,
                       axisLine=c(0.5, "solid", "black"),
                       groupColors=c('#999999','#E69F00')
  ) 
  if (isTRUE(labels)){
    p <- p + theme(text = element_text(size=7),
                   axis.text.x = element_text(angle=90, hjust=1)) 
  }
  print(p)
}


