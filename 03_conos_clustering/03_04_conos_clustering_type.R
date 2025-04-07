
# author: Amparo Roig Adam
# date: 2021/06/24
# description: dimensionality reduction (UMAP), clustering of neuronal types

# inputs: seurat object containing clustering of cell classes
# outputs: conos clustering object, UMAP plots, clustering stability plots, seurat objects containing computed UMAP, clusters and manual annotation of cell classes


renv::activate()

library(Seurat)
library(dplyr)
library(ggplot2)
library(conos)
library(sccore)

#Load precalculated conos cell class clustering & data prep-----------------
load("results/class_stxbp1_conos/conos_clust2.RData")

stxbp1.gluta <- stxbp1.seurat[,stxbp1.seurat$clusters%in%c(1,5,6)]
stxbp1.gaba <- stxbp1.seurat[,stxbp1.seurat$clusters%in%c(9)]

#  change to run clustering on either gaba or glutamatergic neurons
stxbp1.seurat0<-stxbp1.seurat
stxbp1.seurat <- stxbp1.gluta

##Separate into seurat objects by sample and keep in a list
panel.preprocessed <- sapply(as.character(unique(stxbp1.seurat$SampleID)), function(id) {subset(stxbp1.seurat,subset= SampleID ==id)@assays[["RNA"]]@counts}, USE.NAMES = T)

###check data is correct - should be a list with one sample per element
str(panel.preprocessed, 1)

### check cell names are unique
any(duplicated(unlist(lapply(panel.preprocessed,colnames))))

#Basic Seurat preprocessing - include generation of UMAP
panel.preprocessed.seurat <- lapply(panel.preprocessed, function(x){basicSeuratProc(x, verbose = F, n.pcs = 50, tsne = F, umap = T)})#very low n cells in gaba generates error in npc 50, 40 works

###Add metadata
for(id in as.character(unique(stxbp1.seurat$SampleID))){
  print(id)
  panel.preprocessed.seurat[[id]]@meta.data <- cbind(panel.preprocessed.seurat[[id]]@meta.data,subset(stxbp1.seurat, subset= SampleID ==id)@meta.data[,4:17])
}

#Intergate datasets with conos----------
con <- Conos$new(panel.preprocessed.seurat, embedding = 'UMAP',n.cores=1)

##Plot sample-specifc clusters
con$plotPanel( use.local.clusters=TRUE, embedding = 'umap', title.size=3)
ggsave('results/type_stxbp1_conos/sample_cluster_panel.png')

#build graph
gc(gc)
con$pairs$PCA <- NULL
con$pairs$genes <- NULL
con$buildGraph(k=30, k.self=10, space='PCA', n.odgenes= 1000,  alignment.strength = 0.3, verbose=TRUE)

##graph embedding
con$embedGraph(method="UMAP", min.dist=0.5, spread=5 , min.prob.lower=1e-7)

con$plotGraph(embedding = 'UMAP',color.by = 'sample', mark.groups = F, show.legend = T,alpha=0.5)
ggsave(filename = paste('results/type_stxbp1_conos/umap_ sampleID', '.png', sep =''), width = 9, height = 8, device = 'png')

##Get global clusters
con$findCommunities(method=leiden.community, resolution=1, test.stability = T)
con$plotGraph(embedding = 'UMAP',color.by = 'cluster', mark.groups = F, show.legend = T,alpha=0.5)

ggsave(filename = 'results/type_stxbp1_conos/leiden_subcluster.png', width = 8, height = 8, device = 'png')

#Plot cluster stability 
con$plotClusterStability(clustering = 'leiden', what = 'all')
ggsave(paste('results/type_stxbp1_conos/leiden1stability', '.png', sep =''), height = 6, width = 9)

#Plot cluster composition
con.plt <- con$clone(T)
con.plt$clusters <- con$clusters['leiden']
plotClusterBarplots(con)
ggsave('results/type_stxbp1_conos/leiden1barplot.png', height = 11, width = 9)


#Add embedding to Seurat object---------
#set conos computed umap in seurat obj
stxbp1.seurat[["umap"]]<-  CreateDimReducObject(embeddings = con[["embeddings"]][["UMAP"]], key = "UMAP_", assay = DefaultAssay(stxbp1.seurat))
#set cluster names
stxbp1.seurat$`leiden1` <- con[["clusters"]][["leiden"]][['groups']]

DimPlot(stxbp1.seurat, reduction = "umap",group.by = 'SampleID', shuffle = T,label = F, pt.size = 0.1)
ggsave(filename = 'results/type_stxbp1_conos/umap_seurat.png', width = 9, height = 8, device = 'png')

DimPlot(stxbp1.seurat, reduction = "umap", group.by = 'leiden1', shuffle = T, pt.size = 0.1 )
ggsave(filename = 'results/type_stxbp1_conos/leiden1_sub_seurat.png', width = 9, height = 8, device = 'png')

#Plot all meta data
ftrs <- colnames(stxbp1.seurat@meta.data)[colnames(stxbp1.seurat@meta.data)%in%c('orig.ident', 'doublet_thershold', 'obs_names', 'RNA_snn_res.0.8', 'predicted_doublet', 'doublet_threshold', 'nCount_RNA', 'nFeature_RNA') == F ]

ftrs_dimpl <- c()
ftrs_ftrpl <- c()

for(ftr in ftrs){
  if (is.numeric(stxbp1.seurat@meta.data[[ftr]]) == T){
    ftrs_ftrpl <-c( ftrs_ftrpl, ftr)
  }
  if (is.numeric(stxbp1.seurat@meta.data[[ftr]]) == F){
    ftrs_dimpl <- c(ftrs_dimpl, ftr)
  }
}

ggsave(FeaturePlot(stxbp1.seurat, reduction = "umap", features = ftrs_ftrpl, pt.size = 0.1, ncol = 3), 
       filename = 'results/type_stxbp1_conos/ftrs_num.png', width = 15)
ggsave(DimPlot(stxbp1.seurat, reduction = "umap", group.by = ftrs_dimpl, shuffle = T, pt.size = 0.1), 
       filename = 'results/type_stxbp1_conos/ftrs_cat.png', width = 17)

#save
save(con, stxbp1.seurat, file = 'results/type_stxbp1_conos/conos_gluta_clust2.RData')
