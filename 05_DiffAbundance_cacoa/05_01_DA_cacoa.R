# author: Amparo Roig Adam

# description: Run cacoa to measure and do statistics in cell type compositional changes between genotypes

# inputs: cell class and neuron type annotated seurat object 
# outputs: cell loading FC and p value of each cell type comparing HZ to WT samples

library(cacoa)
library(Seurat)
library(coda.base)
library(ggplot2)

#load previosly annotated data (each dataset has to be run throu the code below)
stxbp1.seurat <-  readRDS('results/annotationsmerged_seurat.rds')
stxbp1.gluta.seurat <- readRDS('results/20220309_gluta_subclass_stxbp1_scmap/stxbp1_fullannotation.rds')
stxbp1.gaba.seurat <- readRDS('results/20220307_gaba_subclass_stxbp1_scmap/stxbp1_fullannotation.rds')

#metadata
cao.dataset<- stxbp1.gluta.seurat  #do this for full dataset and class annotation, gaba and glutamateric only datasets and scmap_subclass annotationa
sample.groups=as.character(cao.dataset@meta.data[['Genotype']])
cell.groups=as.character(cao.dataset@meta.data[['scmap_subclass']])
sample.per.cell=as.character(cao.dataset@meta.data[['SampleID']])

names(sample.groups) = names(cell.groups) = names(sample.per.cell) = cao.dataset@meta.data[['SampleID']]

#composition analysis
cao.stxbp1 <- Cacoa$new(NULL, sample.groups, cell.groups, 
                        sample.per.cell, ref.level='WT', target.level='HZ' 
                        , embedding =as.data.frame(cao.dataset@reductions$umap@cell.embeddings))

cao.stxbp1$estimateCellLoadings(n.boot = 1000)
cao.stxbp1$plotCellLoadings(show.pvals=T, alpha = 0.2)
cao.stxbp1$plotCellGroupSizes(proportions = T, show.significance = T)
