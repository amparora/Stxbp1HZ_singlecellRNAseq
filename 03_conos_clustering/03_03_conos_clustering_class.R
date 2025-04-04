
# author: Amparo Roig Adam
# date: 2021/06/24
# description: dimensionality reduction (UMAP), clustering and manual annotation of cell classes

# inputs: seurat object containing full dataset
# outputs: conos clustering object, UMAP plots, clustering stability plots, seurat objects containing computed UMAP, clusters and manual annotation of cell classes


renv::activate()

library(Seurat)
library(dplyr)
library(ggplot2)
library(conos)
library(sccore)

#Input & data prep-----------------
load("results/20210621_conos_inputs.RData")

#stxbp1.gluta <- stxbp1.seurat[,stxbp1.seurat$clusters%in%c(1,5,6)]
#stxbp1.gaba <- stxbp1.seurat[,stxbp1.seurat$clusters%in%c(9)]

#use gluta gells
stxbp1.seurat0 <- stxbp1.seurat
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
ggsave('results/class_stxbp1_conos/cell_classes/sample_cluster_panel.png')

#build graph
gc(gc)
con$pairs$PCA <- NULL
con$pairs$genes <- NULL
con$buildGraph(k=30, k.self=10, space='genes', n.odgenes= 2000,  alignment.strength = 0.3, verbose=TRUE)

##graph embedding
con$embedGraph(method="UMAP", min.dist=0.5, spread=5 , min.prob.lower=1e-7)

con$plotGraph(embedding = 'UMAP',color.by = 'sample', mark.groups = F, show.legend = T,alpha=0.5)
ggsave(filename = paste('results/class_stxbp1_conos/umap_ sampleID', '.png', sep =''), width = 9, height = 8, device = 'png')

##Get global clusters
con$findCommunities(method=leiden.community, resolution=1, test.stability = T)
con$plotGraph(embedding = 'UMAP',color.by = 'cluster', mark.groups = F, show.legend = T,alpha=0.5)

#Renfine specific clusters
con[["clusters"]][["leiden1_org"]] <- con[["clusters"]][["leiden"]] #save original clusters
  #subcluster 2
con[["clusters"]][["leiden"]][["groups"]] <- as.factor(findSubcommunities(con, target.clusters = '2'))
con[["clusters"]][["leiden"]][["groups"]] <- as.factor(findSubcommunities(con, target.clusters = '2_1'))

  #join too many subdiviaions
levels(con[["clusters"]][["leiden"]][["groups"]])[levels(con[["clusters"]][["leiden1"]][["groups"]])%in%c('2_2','2_0')] <- '2_0'
levels(con[["clusters"]][["leiden"]][["groups"]])[levels(con[["clusters"]][["leiden1"]][["groups"]])%in%c('2_1_0','2_1_2', '2_1_3', '2_1_4')] <- '2_1_0'

levels(con[["clusters"]][["leiden"]][["groups"]])[levels(con[["clusters"]][["leiden1"]][["groups"]]) == '10'] <- '1'

#Plot clusters and stats
con$plotGraph(embedding = 'UMAP',clustering = 'leiden', mark.groups = T, show.legend = F,alpha=1)
ggsave(filename = 'results/class_stxbp1_conos/cell_classes/leiden_subcluster.png', width = 8, height = 8, device = 'png')

#Plot cluster stability 
con$plotClusterStability(clustering = 'leiden', what = 'all')
ggsave(paste('results/class_stxbp1_conos/cell_classes/leiden1stability', '.png', sep =''), height = 6, width = 9)

#Plot cluster composition
con.plt <- con$clone(T)
con.plt$clusters <- con$clusters['leiden']
plotClusterBarplots(con)
ggsave('results/class_stxbp1_conos/cell_classes/leiden1barplot.png', height = 11, width = 9)


#Add embedding to Seurat global object---------
#set conos computed umap in seurat obj
stxbp1.seurat[["umap"]]<-  CreateDimReducObject(embeddings = con[["embeddings"]][["UMAP"]], key = "UMAP_", assay = DefaultAssay(stxbp1.seurat))
#set cluster names
stxbp1.seurat$`leiden1` <- con[["clusters"]][["leiden"]][['groups']]

DimPlot(stxbp1.seurat, reduction = "umap",group.by = 'SampleID', shuffle = T,label = F, pt.size = 0.1)
ggsave(filename = 'results/class_stxbp1_conos/cell_classes/umap_seurat.png', width = 9, height = 8, device = 'png')

DimPlot(stxbp1.seurat, reduction = "umap", group.by = 'leiden1', shuffle = T, pt.size = 0.1 )
ggsave(filename = 'results/class_stxbp1_conos/cell_classes/leiden1_sub_seurat.png', width = 9, height = 8, device = 'png')

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
       filename = 'results/class_stxbp1_conos/cell_classes/ftrs_num.png', width = 15)
ggsave(DimPlot(stxbp1.seurat, reduction = "umap", group.by = ftrs_dimpl, shuffle = T, pt.size = 0.1), 
       filename = 'results/class_stxbp1_conos/cell_classes/ftrs_cat.png', width = 17)
#marker genes
FeaturePlot(stxbp1.seurat, reduction = "umap", features = c('Slc17a7', 'Gad1', 'Aqp4', 'Pdgfra', 'Mog', 'Itgam', 'Flt1', 'Pdgfrb', 'Spp1'), pt.size = 0.1)
ggsave(filename = 'results/class_stxbp1_conos/cell_classes/markers.png', device = 'png')

#Manual annotation
stxbp1.seurat$class <- stxbp1.seurat$leiden1
  #manual annotate cluster names intro annotation name according to markers above - the order of these may change in different runs
levels(stxbp1.seurat$class) <- c('Glutamatergic neurons', 'Endothelial cells', 'Pericytes', 'VLMCs', 'Astrocytes', 'Glutamatergic neurons', 'Microglia', 'OPCs', 'Oligodendrocytes', 'Glutamatergic neurons', 'GABAergic neurons')
DimPlot(stxbp1.seurat, reduction = "umap", group.by = 'leiden1', shuffle = T, label = T,repel = T, pt.size = 0.1 )
ggsave(filename = 'results/class_stxbp1_conos/cell_classes/leiden1_seuratumap.png', width = 9, height = 8, device = 'png')

#save
save(con, stxbp1.seurat, file = 'results/class_stxbp1_conos/cell_classes/conos_clust2.RData')
