
# author: Amparo Roig Adam
# date: 2021/06/24
# description: dimensionality reduction (UMAP) and clustering of neuronal cell types

# inputs: seurat object containing clustering of cell classes
# outputs: conos clustering object, UMAP plots, clustering stability plots, seurat objects containing calculated UMAP and clusters 

renv::activate()

library(Seurat)
library(dplyr)
library(ggplot2)
library(conos)
library(sccore)

#Input & data prep-----------------
load("results/20210621_conos_inputs.RData")

#from cell classes filter GABA or gluta neurons according to marker expression in clusters
#stxbp1.gluta <- stxbp1.seurat[,stxbp1.seurat$clusters%in%c(1,5,6)]
#stxbp1.gluta <- stxbp1.seurat[,stxbp1.seurat$clusters%in%c(9)]

#use gluta gells
stxbp1.seurat0 <- stxbp1.seurat
stxbp1.seurat <- stxbp1.gluta #change this to use gluta or gaba cells only

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
ggsave('results/20220301_gluta_subclass_stxbp1_conos/cell_classes/sample_cluster_panel.png')

con.test <- con

#build graph
gc(gc)
con.test$pairs$PCA <- NULL
con.test$pairs$genes <- NULL
con.test$buildGraph(k=30, k.self=10, space='genes', n.odgenes= 2000,  alignment.strength = 0.3, verbose=TRUE)

##graph embedding
con.test$embedGraph(method="UMAP", min.dist=0.5, spread=5 , min.prob.lower=1e-7)

con.test$plotGraph(embedding = 'UMAP',color.by = 'sample', mark.groups = F, show.legend = T,alpha=0.5)
ggsave(filename = paste('results/20220301_gluta_subclass_stxbp1_conos/umap_ sampleID', '.png', sep =''), width = 9, height = 8, device = 'png')

##Get global clusters
#con.test$testStability(method=leiden.community, leidres=2, name = 'leiden2',test.stability = T, stability.subsamples = 100)
con.test$findCommunities(method=leiden.community, resolution=1, test.stability = T)

#Renie specific clusters
con.test[["clusters"]][["leiden1_org"]] <- con.test[["clusters"]][["leiden1"]]
con.test[["clusters"]][["leiden1"]][["groups"]] <- as.factor(findSubcommunities(con, target.clusters = '2'))
con.test[["clusters"]][["leiden1"]][["groups"]] <- as.factor(findSubcommunities(con, target.clusters = '2_1'))
con.test[["clusters"]][["leiden1_sub"]] <- con.test[["clusters"]][["leiden1"]]

levels(con.test[["clusters"]][["leiden1"]][["groups"]])[levels(con.test[["clusters"]][["leiden1"]][["groups"]])%in%c('2_2','2_0')] <- '2_0'
levels(con.test[["clusters"]][["leiden1"]][["groups"]])[levels(con.test[["clusters"]][["leiden1"]][["groups"]])%in%c('2_1_0','2_1_2', '2_1_3', '2_1_4')] <- '2_1_0'

levels(con.test[["clusters"]][["leiden1"]][["groups"]])[levels(con.test[["clusters"]][["leiden1"]][["groups"]]) == '10'] <- '1'

#Plot clusters and stats
con.test$plotGraph(embedding = 'UMAP',clustering = 'leiden1', mark.groups = T, show.legend = F,alpha=1)
ggsave(filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/leiden2.png', width = 8, height = 8, device = 'png')

#Plot cluster stability 
con.test$plotClusterStability(clustering = 'leiden1', what = 'all')
ggsave(paste('results/20220301_gluta_subclass_stxbp1_conos/cell_classes/leiden1stability', '.png', sep =''), height = 6, width = 9)

#Plot cluster composition
con.plt <- con.test$clone(T)
con.plt$clusters <- con.test$clusters[5]
plotClusterBarplots(con.test)
ggsave('results/20220301_gluta_subclass_stxbp1_conos/cell_classes/leiden1barplot.png', height = 11, width = 9)


#Add embedding to Seurat global object---------
stxbp1.seurat <- NormalizeData(object = stxbp1.seurat)
#stxbp1.seurat <- FindVariableFeatures(object = stxbp1.seurat)
stxbp1.seurat <- ScaleData(object = stxbp1.seurat)
# stxbp1.seurat <- RunPCA(object = stxbp1.seurat)
# stxbp1.seurat <- FindNeighbors(object = stxbp1.seurat)
# stxbp1.seurat <- FindClusters(object = stxbp1.seurat)
# stxbp1.seurat <- RunUMAP(object = stxbp1.seurat, dims = 1:30)
# DimPlot(object = stxbp1.seurat, reduction = "umap")
#change class of metadata Live Cells percent
stxbp1.seurat$LiveCellsPercent <- as.numeric(as.character(stxbp1.seurat$LiveCellsPercent))

stxbp1.seurat[["umap"]]<-  CreateDimReducObject(embeddings = con.test[["embeddings"]][["UMAP"]], key = "UMAP_", assay = DefaultAssay(stxbp1.seurat))
#set cluster names
stxbp1.seurat$`leiden1` <- con.test[["clusters"]][["leiden1"]][['groups']]

DimPlot(stxbp1.seurat, reduction = "umap",group.by = 'SampleID', shuffle = T,label = F, pt.size = 0.1)
ggsave(filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/umap_seurat.png', width = 9, height = 8, device = 'png')


DimPlot(stxbp1.seurat, reduction = "umap", group.by = 'leiden1', shuffle = T, pt.size = 0.1 )
ggsave(filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/leiden1_sub_seurat.png', width = 9, height = 8, device = 'png')

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
       filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/ftrs_num.png', width = 15)
ggsave(DimPlot(stxbp1.seurat, reduction = "umap", group.by = ftrs_dimpl, shuffle = T, pt.size = 0.1), 
       filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/ftrs_cat.png', width = 17)
#marker genes
FeaturePlot(stxbp1.seurat, reduction = "umap", features = c('Slc17a7', 'Gad1', 'Aqp4', 'Pdgfra', 'Mog', 'Itgam', 'Flt1', 'Pdgfrb', 'Spp1'), pt.size = 0.1)
ggsave(filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/markers.png', device = 'png')

#Manual annotation
stxbp1.seurat$class <- stxbp1.seurat$leiden1
levels(stxbp1.seurat$class) <- c('Glutamatergic neurons', 'Endothelial cells', 'Pericytes', 'VLMCs', 'Astrocytes', 'Glutamatergic neurons', 'Microglia', 'OPCs', 'Oligodendrocytes', 'Glutamatergic neurons', 'GABAergic neurons')
DimPlot(stxbp1.seurat, reduction = "umap", group.by = 'leiden1', shuffle = T, label = T,repel = T, pt.size = 0.1 )
ggsave(filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/leiden1_seuratumap.png', width = 9, height = 8, device = 'png')

#markers = c('Snap25', 'Gabrb2','Gad1', 'Vip', 'Sst', 'Pvalb','Slc17a7','Aqp4','Gja1','Pdgfra', 'Grm5', 'Ccnb1','Mog', 'Plp1', 'Cnp','Itgam', 'Cd74', 'Aif1','Flt1', 'Cldn5', 'Icam2','Osr1', 'Spp1','Pdgfrb', 'Notch3')
#DoHeatmap(stxbp1.seurat, features = markers, size = 3, group.bar = T,group.by = 'class', slot = 'data')
#ggsave(filename = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/annotated_heatmap_seurat.png', width = 9, height = 8, device = 'png')

##Number of cells/class wt vs hz
ggplot(data = stxbp1.seurat@meta.data, aes(x = class, fill = Genotype))+
  geom_bar()+
  scale_y_continuous(n.breaks =15)+
  ylab('n cells') + theme(axis.text.x =element_text(angle = 45, hjust = 1, size = 11))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pie(table(stxbp1.seurat$class), labels = paste(names(table(stxbp1.seurat$class)), table(stxbp1.seurat$class), sep =' - '), border="white", col = gg_color_hue(9))



#save for reload later
save(con.test, con, stxbp1.seurat, file = 'results/20220301_gluta_subclass_stxbp1_conos/cell_classes/conos_clust2.RData')
