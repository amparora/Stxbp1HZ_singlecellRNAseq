# author: Amparo Roig Adam with support from Lisa Bast
# date:2022/02/23
# description: Run scmap reference cluster annotation prameter screen

# inputs: filtered Allen Brain dataset dataset in SCE file format; HDO and HVGs files precalculated in 04_03_scmap_HDOgenes.R and 04_04_scmap_HCGs.py 
# outputs: file and barplots with performance metrics of all tested parameter combination
############
library(SingleCellExperiment)
library(scmap)
library(Seurat)
set.seed(2022)
############
##opts - make F/T to test a subset of the reference set or check overlap of mannual annotation of Stxbp1 +/- with scmap annotation 
ref_selftest = FALSE #test classes and clusters in reference set 1/5 test and 5/4 train set
ref_stxbp1test = TRUE #test class annotation in Stxbp1 +/- dataset with full reference test

if(ref_selftest){
    annotation_levels=c('class', 'subclass', 'cluster')
    #directory to read and save
    outdir <- 'results/ALL_ref_scmap_lognorm/'
}
if(ref_stxbp1test){
    annotation_levels=c('class')
    #directory to read and save
    outdir <- 'results/class_stxbp1_scmap/'
}

#input

ref_ABM_S1_10x <- readRDS("data/ALLEN_RefData/20220217_lognorm_ref_ABM_S1_10x.rds")
rowData(ref_ABM_S1_10x)$feature_symbol <- rownames(ref_ABM_S1_10x)

#trainset/testset
## test annotation in reference dataset
if (ref_selftest){
    test_ref_ABM_S1_10x <- ref_ABM_S1_10x[,sample(colnames(ref_ABM_S1_10x), size = ncol(ref_ABM_S1_10x)/5)]
    test_annotation_set <- test_ref_ABM_S1_10x 

    train_ref_ABM_S1_10x <- ref_ABM_S1_10x[,colnames(ref_ABM_S1_10x)%in%colnames(test_ref_ABM_S1_10x) == F]
    train_ref_ABM_S1_10x 
}

##test class annotation in our data
if (ref_stxbp1test){
  #trainset is full reference dataset
  train_ref_ABM_S1_10x <- ref_ABM_S1_10x
  #test set is full our dataset
  load('results/20210624_conos_final/cell_classes/conos_clust2.RData')
  #Manual annotation
  stxbp1.seurat$class_label <- stxbp1.seurat$class
  #fix names to coincide with classes in reference (Glutamatergic, GABAergic, Non-Neuronal)
  levels(stxbp1.seurat$class_label)[levels(stxbp1.seurat$class)=='Glutamatergic neurons']<- 'Glutamatergic'
  levels(stxbp1.seurat$class_label)[levels(stxbp1.seurat$class)=='GABAergic neurons']<- 'GABAergic'
  levels(stxbp1.seurat$class_label)[levels(stxbp1.seurat$class)%in%c('GABAergic neurons', 'Glutamatergic neurons')==F]<- 'Non-Neuronal'
  test_annotation_set <- as.SingleCellExperiment(stxbp1.seurat)
  rowData(test_annotation_set)$feature_symbol <- rownames(test_annotation_set)
  rm(con.test, con)
  }

##store calculations with all methods
params.performance <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(params.performance) <- c("mapping", "ftr.selection", "nGenes", "annotation.level","output","percent.cells")                               

#initiate file to append results
write.table(data.frame("mapping", "ftr.selection", "nGenes", "annotation.level","output","percent.cells"), file = paste(outdir,"scmap_performance.csv", sep=''), sep = ',', col.names = F,row.names = F)

for(nftrs in c(500,1000,2000,3000,5000)){
  for(selmethod in c('HDOgenes','seurat_v3_HVGs', 'cell_ranger_HVGs')){
    for(mapping in c('scmapCell2Cluster', 'scmapCluster')){
      for(annotation.level in annotation_levels){
        
        print(paste(nftrs, selmethod, mapping, annotation.level))
        
        #read genes from reference
        sce <- train_ref_ABM_S1_10x[!duplicated(rownames(train_ref_ABM_S1_10x)), ]
        
        ftrs.file <- paste('results/20220217_ALL_ref_scmap_lognorm/', nftrs, '_', selmethod,'_ref_ABM_S1_10x.txt', sep = '')
        select.ftrs <- read.table(ftrs.file)$V1
        
        genes_priority = c('Snap25','Gabrb2', 'Gad1', 'Vip', 'Sst', 'Pvalb', 'Slc17a7', 'Aqp4','Gja1', 'Siglech','Serpinb1a', 'Pdgfra', 'Grm5', 'Ccnb1', 
                           'Cspg4', 'Enpp6', 'Mog', 'Plp1', 'Cnp', 'Itgam', 'Cd74', 'Aif1', 'Ccl4','Cldn5','Flt1', 'Cldn5', 'Icam2', 
                           'Acta2', 'Bgn', 'Mrc1', 'Kcnj8', 'Ctla2a', 'Mbp', 'Ogn','Osr1', 'Spp1','Pdgfrb', 'Notch3')
        
        sce<-setFeatures(sce, features = c(select.ftrs, genes_priority))
        
        #propagate selected ftrs to test set
        test_annotation_set <- setFeatures(test_annotation_set, features = c(select.ftrs, genes_priority))
        
        ### ANNOTATION
        if(mapping == 'scmapCluster'){
          #index
          sce <- indexCluster(sce, cluster_col = paste(annotation.level,'_label', sep = ''))
          
          test_annotation_set <- indexCluster(test_annotation_set, cluster_col =  paste(annotation.level,'_label', sep = ''))
          
          #projection
          scmap_results <- scmapCluster(
            projection = test_annotation_set, 
            index_list = list(
              train = metadata(sce)$scmap_cluster_index
            ), threshold = 0.5
          )
          #head(scmapCluster_results$scmap_cluster_labs)
          
        }
        
        if(mapping == 'scmapCell2Cluster'){
          #index
          sce <- indexCell(sce)
          #projection
          scmap_results <- scmapCell(projection = test_annotation_set, list(train = metadata(sce)$scmap_cell_index))
          #annotation
          scmap_results <- scmapCell2Cluster(
            scmapCell_results = scmap_results, 
            cluster_list = list(
              as.character(colData(sce)[, paste(annotation.level,'_label', sep = '')])
            )
          )
        }
        #save performance for posterior comparison
        nolab <- 100*sum(scmap_results$scmap_cluster_labs == 'unassigned')/length(scmap_results$scmap_cluster_labs)
        oklab<- 100*sum(scmap_results$scmap_cluster_labs ==test_annotation_set@colData@listData[[paste(annotation.level,'_label', sep = '')]])/length(scmap_results$scmap_cluster_labs)
        wronglab <- 100*sum(scmap_results$scmap_cluster_labs != test_annotation_set@colData@listData[[paste(annotation.level,'_label', sep = '')]] & scmap_results$scmap_cluster_labs != 'unassigned')/length(scmap_results$scmap_cluster_labs)
        
        
        #store performance in data frame
        params.performance <- rbind(params.performance, data.frame(mapping = mapping, ftr.selection = selmethod, nGenes = nftrs, annotation.level = annotation.level,
                                                                   output = 'unassigned',percent.cells = nolab))
        params.performance <- rbind(params.performance, data.frame(mapping = mapping, ftr.selection = selmethod, nGenes = nftrs, annotation.level = annotation.level,
                                                                   output = 'correct',percent.cells = oklab))
        params.performance <- rbind(params.performance, data.frame(mapping = mapping, ftr.selection = selmethod, nGenes = nftrs, annotation.level = annotation.level,
                                                                   output = 'wrong',percent.cells = wronglab))
        
        write.table(data.frame(mapping = mapping, ftr.selection = selmethod, nGenes = nftrs, annotation.level = annotation.level,
                                output = 'unassigned',percent.cells = nolab), file = paste(outdir,"scmap_performance.csv", sep =''), sep = ',', col.names = F,row.names = F, append = T)
        write.table(data.frame(mapping = mapping, ftr.selection = selmethod, nGenes = nftrs, annotation.level = annotation.level,
                                output = 'correct',percent.cells = oklab), file = paste(outdir,"scmap_performance.csv", sep =''), sep = ',', col.names = F,row.names = F, append = T)
        write.table(data.frame(mapping = mapping, ftr.selection = selmethod, nGenes = nftrs, annotation.level = annotation.level,
                                output = 'wrong',percent.cells = wronglab), file = paste(outdir,"scmap_performance.csv", sep =''), sep = ',', col.names = F,row.names = F, append = T)
        # 
        rm(scmap_results)
        
 }}}
  gc(gc(verbose = F))
  }

#plot performance
library(ggplot2)
library(dplyr)
library(ggthemes)
theme_set(theme_minimal())
theme_update(plot.background = element_rect(fill = "transparent", color = NA))

scmap.test.performance <- read.table(paste(outdir,"scmap_performance.csv", sep =''), header = T, sep = ',')
scmap.test.performance$mod <- paste(scmap.test.performance$mapping, scmap.test.performance$ftr.selection)


plot<- scmap.test.performance[,]%>% arrange(nGenes) %>% mutate(nGenes = factor(as.character(nGenes), levels=c('500', '1000', '2000', '3000', '5000'))) %>%
ggplot( aes(x = nGenes, y = percent.cells, fill = mod))+
  geom_bar(stat = 'identity', position = 'dodge') +ylim(c(0,100))+facet_wrap(annotation.level~output)

plot

ggsave(paste(outdir,"plot_results.svg", sep =''), plot = plot, bg='transparent', width = 18, height = 15)

plot(
  getSankey(
    colData(test_annotation_set)$class_label, 
    scmap_results$scmap_cluster_labs[,'train'],
    plot_height = 500
  )
)


