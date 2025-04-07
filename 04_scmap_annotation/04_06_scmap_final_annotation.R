# author: Amparo Roig Adam with support from Lisa Bast
# date:2022/03/04
# description: Run scmap reference cluster annotation of GABA and Glutamatergic neurons in Stxbp1+/- dataset (inc wt samples)

# inputs: filtered Allen Brain dataset dataset in SCE file format; HVGs files precalculated 04_04_scmap_HCGs.py for either gaba or glutamatergic neurons
# outputs: seurat object with annotated neuronal types in metadata slot 'scmap_subclass'
############
library(scmap)
library(SummarizedExperiment)
library(Seurat)
library(ggplot2)
#############
#load either GABA or Glutamatergic subset for annotation
load('results/20220301_gluta_subclass_stxbp1_conos/cell_classes/conos_clust_1000pca.RData')
load('results/20220304_gaba_subclass_stxbp1_conos/cell_classes/conos_clust_1000pca.RData')

DimPlot(stxbp1.seurat, group.by = 'clusters')

test_stxbp1 <- as.SingleCellExperiment(stxbp1.seurat)
rowData(test_stxbp1)$feature_symbol <- rownames(test_stxbp1)

#reference
ref_ABM_S1_10x <- readRDS("data/ALLEN_RefData/20220217_lognorm_ref_ABM_S1_10x.rds")
rowData(ref_ABM_S1_10x)$feature_symbol <- rownames(ref_ABM_S1_10x)
train_ref_ABM_S1_10x <- ref_ABM_S1_10x[,colnames(ref_ABM_S1_10x)[ref_ABM_S1_10x$class_label =='Glutamatergic']]

#all classes together
nftrs = 2000
selmethod = 'seurat_v3_HVGs'
mapping = 'scmapCell2Cluster'
      
print(paste(nftrs, selmethod, mapping))

annotation.level <- 'subclass'
      
#read genes from reference
sce <- train_ref_ABM_S1_10x[!duplicated(rownames(train_ref_ABM_S1_10x)), ]
      
ftrs.file <- paste('results/20220307_gaba_subclass_stxbp1_scmap/', nftrs, '_', selmethod,'_ref_ABM_S1_10x.txt', sep = '')
      
select.ftrs <- read.table(ftrs.file)$V1

    #genes to add       
genes_priority = c('Slc17a7', 'Rasgrf2', 'Cux2', 'Rorb', 'Rspo1', 'PLcxd2', 'Thsd7a', 'Kcnk2', 'Sulf2', 'Foxp2', 'Syt6', 'Rprm', 'Cplx3', 'Lgr5')
            
genes_priority <- unique(c(genes_priority,suptyp.genes)[c(genes_priority,suptyp.genes)%in%rownames(stxbp1.seurat)])
      
sce<-setFeatures(sce, features = c(select.ftrs, genes_priority))
      
#propagate selected ftrs to test set
test_stxbp1 <- setFeatures(test_stxbp1, features = c(select.ftrs, genes_priority))
      
### ANNOTATION
mapping = 'scmapCell2Cluster'
#index
sce <- indexCell(sce)
#projection
scmap_results <- scmapCell(projection = test_stxbp1, list(train = metadata(sce)$scmap_cell_index))

#annotation
scmap_results_cluster <- scmapCell2Cluster(
    scmapCell_results = scmap_results, threshold = 0.2,
    cluster_list = list(
    as.character(colData(sce)$subclass_label))
    )

#saveRDS(scmap_results_cluster, 'results/20220307_gaba_subclass_stxbp1_scmap/subclass_annotation.rds')

plot(
  getSankey(
     
    scmap_results_cluster$scmap_cluster_labs[,'train'],
    colData(test_stxbp1)$leiden1,
    plot_height = 500, 
    #colors = ct_cols$ct_cols[ct_cols$cts%in%unique(scmap_results_cluster$scmap_cluster_labs[,'train'])]
  ))

stxbp1.seurat$scmap_subclass <- NA
stxbp1.seurat$scmap_subclass[colnames(test_stxbp1)] <- scmap_results_cluster$scmap_cluster_labs[,'train']

stxbp1.seurat$scmap_supertype <- NA
stxbp1.seurat$scmap_supertype[colnames(test_stxbp1)] <- scmap_results_cluster$scmap_cluster_labs[,'train']

stxbp1.seurat$scmap_un_subclass <- 'annotated'
stxbp1.seurat$scmap_un_subclass[stxbp1.seurat$scmap_subclass == 'unassigned'] <- 'unassigned'

DimPlot(object = stxbp1.seurat, reduction = 'umap', group.by = 'scmap_subclass', label = F,shuffle = T)
DimPlot(object = stxbp1.seurat, reduction = 'umap', group.by = 'scmap_subclass', label = F,shuffle = T, cells = colnames(stxbp1.seurat[,stxbp1.seurat$scmap_subclass == 'unassigned']))
DimPlot(object = stxbp1.seurat, reduction = 'umap', group.by = 'scmap_un_subclass', label = F,shuffle = T, cols = c('lightgrey', 'red'))

ggsave(filename = 'results/20220307_gaba_subclass_stxbp1_scmap/umap_scmap_subclass.png', width = 8, height = 6)

saveRDS(stxbp1.seurat, 'results/20220307_gaba_subclass_stxbp1_scmap/stxbp1_fullannotation.rds')

png('results/20220307_gaba_subclass_stxbp1_scmap/scmap_subclass_similarity_hist.png')
hist(scmap_results_cluster[["scmap_cluster_siml"]], main = 'Similarity to assigned annotation - subclass')
dev.off()

