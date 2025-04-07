# author: Amparo Roig Adam with support from Lisa Bast
# date:2022/02/17
# description: Calculate HDGs for scmap automated cell type annotation prameter screen

# inputs: filtered Allen Brain dataset dataset in SCE file format 
# outputs: text files with cells used as 'training set' and HDO genes calculated for parameter screen

library(SingleCellExperiment)
library(scmap)
set.seed(2022)

######include options to used reference for testing the methods or big classes to test reference in or data
#input

ref_ABM_S1_10x <- readRDS("data/ALLEN_RefData/20220217_lognorm_ref_ABM_S1_10x.rds")
rowData(ref_ABM_S1_10x)$feature_symbol <- rownames(ref_ABM_S1_10x)

# NOTE - Use only GABA or Glutamatergic neurons when runnning the final analysis of neuron type annotation - keep commented for parameter screen
#GABA_ref_ABM_S1_10x <- ref_ABM_S1_10x[,colnames(ref_ABM_S1_10x)[ref_ABM_S1_10x$class_label =='GABAergic']]
#table(GABA_ref_ABM_S1_10x$subclass_label)
#ref_ABM_S1_10x <- GABA_ref_ABM_S1_10x

#trainset/testset
test_ref_ABM_S1_10x <- ref_ABM_S1_10x[,sample(colnames(ref_ABM_S1_10x), size = ncol(ref_ABM_S1_10x)/5)]
test_ref_ABM_S1_10x 

train_ref_ABM_S1_10x <- ref_ABM_S1_10x[,colnames(ref_ABM_S1_10x)%in%colnames(test_ref_ABM_S1_10x) == F]
train_ref_ABM_S1_10x 

  #write file with train set cell IDs for scanpy HVGs calculation in the exact same train set
write.table(colnames(train_ref_ABM_S1_10x), file = "results/20220217_ALL_ref_scmap_lognorm/cellID_train_ref_ABM_S1_10x.txt", sep = '\n', quote = F, row.names = F, col.names = F)

#set dataset to calculate genes
sce <- train_ref_ABM_S1_10x[!duplicated(rownames(train_ref_ABM_S1_10x)), ]

outfile <- 'HDOgenes_ref_ABM_S1_10x'

for(nftrs in c(500, 1000, 2000, 3000, 5000)){
  print(nftrs)
  sce <- selectFeatures(sce, n_features = nftrs, suppress_plot = F)
  
  #main marker gene forced to be included
  genes_priority = c('Snap25','Gabrb2', 'Gad1', 'Vip', 'Sst', 'Pvalb', 'Slc17a7', 'Aqp4','Gja1', 'Siglech','Serpinb1a', 'Pdgfra', 'Grm5', 'Ccnb1', 
                     'Cspg4', 'Enpp6', 'Mog', 'Plp1', 'Cnp', 'Itgam', 'Cd74', 'Aif1', 'Ccl4','Cldn5','Flt1', 'Cldn5', 'Icam2', 
                     'Acta2', 'Bgn', 'Mrc1', 'Kcnj8', 'Ctla2a', 'Mbp', 'Ogn','Osr1', 'Spp1','Pdgfrb', 'Notch3')
  
  sce@rowRanges@elementMetadata[sce@rowRanges@elementMetadata$feature_symbol%in%genes_priority,]$scmap_features <- TRUE
  
  ## save list of selected genes
  nftrs.outfile = paste('results/20220217_ALL_ref_scmap_lognorm/', nftrs, '_',outfile,'.txt', sep = '')
  write.table(as.character(sce@rowRanges@elementMetadata[sce@rowRanges@elementMetadata$scmap_features==T,]$feature_symbol), file = nftrs.outfile, sep = '\n', quote = F, row.names = F, col.names = F)
  
  gc(gc())
}
