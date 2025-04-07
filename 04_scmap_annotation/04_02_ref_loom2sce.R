# author: Amparo Roig Adam
# date: 2021/02/17
# description: transform loom file to SCE object

# inputs: filtered Allen Brain dataset dataset in loom file format
# outputs: RData with SingleCellExperiment object to use in scmap. Plots tSNE to double check that the data matrix and metadata are properly matched

##Packages
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(loomR)

#Input data prep-----------------
##Connect with loom file and transform into Seurat object
AMB_S1_10x<- Connect(filename ='STXBP1/data/ALLEN_RefData/20220217_ref_ABM_S1_10x.loom', mode = "r")
AMB_S1_10x

seurat.ref_AMB_S1_10x <- as.Seurat(AMB_S1_10x)

  #log normalisation
seurat.ref_AMB_S1_10x <- NormalizeData(seurat.ref_AMB_S1_10x)
sce.ref_AMB_S1_10x <- as.SingleCellExperiment(seurat.ref_AMB_S1_10x)
AMB_S1_10x$close()

saveRDS(sce.ref_AMB_S1_10x, 'STXBP1/results/20220217_ref_norm_ABM_S1_10x.rds')

##CHECK TSNE AND CELL ANNOTATION LABELS
#read data

tsne.loadings <- read.table("STXBP1/data/ALLEN_RefData/tsne.csv", sep = ',', header = T, row.names = 1)
colnames(tsne.loadings) <- paste0("TSNE_", 1:2)

#add tsne loadings to seurat object
seurat.ref_AMB_S1_10x@reductions[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(tsne.loadings[colnames(seurat.ref_AMB_S1_10x),]), key = "TSNE_", assay = DefaultAssay(seurat.ref_AMB_S1_10x))

DimPlot(object = seurat.ref_AMB_S1_10x, reduction = 'tsne', group.by = 'subclass_label')+NoLegend()

##check metadata is correct --- cells selecteed seem to be wrong - at least the metadata and its match with the cell names
table(colnames(seurat.ref_AMB_S1_10x)==seurat.ref_AMB_S1_10x@meta.data$sample_name)
