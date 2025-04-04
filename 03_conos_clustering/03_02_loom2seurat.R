
# author: Amparo Roig Adam
# date: 14/2/2021
# description: transform loom file to seurat object

# inputs: filtered dataset in loom file format
# outputs: RData with seurat object ready for conos use


##Packages
library(Seurat)
library(SeuratDisk)

#Input data prep-----------------
##Connect with loom file and transform into Seurat object
stxbp1.loom<- Connect(filename ='results/20210604_conos_filt_norm.loom', mode = "r")
stxbp1.loom

stxbp1.seurat <- as.Seurat(stxbp1.loom)

### close loom files when done
stxbp1.loom$close_all()
rm(stxbp1.loom)

#save
save(stxbp1.seurat, file = 'results/20210621_conos_inputs.RData')
