#####
#authtor: Amparo Roig Adam
#created: 2022/05/10
#last update: 2024/04/29
#description: Downsample celltypes for equal power in DEA
#####

library(dplyr)
library(Seurat)
library(ggplot2)
library(forcats)
library(Libra)
library(DESeq2)
library(ggpubr)
set.seed(123)

#plotting settings
theme_set(theme_minimal(base_size = 12))
theme_update(axis.line=element_line(linetype = 1, size = 0.25), axis)

#load
stxbp1.seurat <- readRDS('results/annotationsmerged_seurat.rds')

## paths and opts
PATH_OUT = 'results/20220616_deseq2_subclass_downsampling/'

#annotation level to perform downsampling
level = 'class03' #options: cell classes: 'class03'; neuronal types: 'scmap_subclass'

#thresholds
pthr = 0.05
fcthr = 0.05

#number of cells to downsample to - 250 for classes and 40 for neuronal types
table(as.factor(stxbp1.seurat@meta.data[level][,1]), stxbp1.seurat$Genotype)   # print total cell numbers 

if (level == 'class03'){
  n = 250
}
if (level == 'scmap_subclass'){
  n = 40

        # filter dataset to only neurons if level is subclass
  stxbp1.seurat <- stxbp1.seurat[,stxbp1.seurat$class03 %in% c('Glutamatergic neurons', 'GABAergic neurons')]
}
#####

# Keep only cts that do not reach minimum n in both genotypes
celltypes = table(as.character(stxbp1.seurat@meta.data[level][,1])) > n*2
celltypes = names(celltypes)[celltypes==T]

# set up results table
nDEGs <- data.frame(celltype = NA, nCells = NA, nDEGs = NA, nUpreg = NA, nDownreg = NA, nGtested = NA, round = NA )

## Run DESeq2 in all cell types
for(ct in celltypes){
  print(ct)
  
  #get all HZ and WT cells
  nHZ <- ncol(stxbp1.seurat[,stxbp1.seurat@meta.data[level] == ct & stxbp1.seurat$Genotype == 'HZ'])
  nWT <- ncol(stxbp1.seurat[,stxbp1.seurat@meta.data[level] == ct & stxbp1.seurat$Genotype == 'WT'])
  
  #Generate 100 random subsamples of n cells
  sampleHZ <- replicate(100,sample(1:nHZ, size = n))
  sampleWT <- replicate(100,sample(1:nWT, size = n))
  
    #iterate throught the samples and save DEGs
  for(r in 1:100){
    print(r)
    print(ct)
    
    #select cells
    cells <- colnames(stxbp1.seurat[,stxbp1.seurat@meta.data[level] == ct & stxbp1.seurat$Genotype == 'HZ'])[sampleHZ[,r]]
    cells <- c(cells, colnames(stxbp1.seurat[,stxbp1.seurat@meta.data[level] == ct & stxbp1.seurat$Genotype == 'WT'])[sampleWT[,r]])
    
    ## DESeq2 
    DDS_full <-to_pseudobulk(stxbp1.seurat@assays$RNA@counts[,cells], meta = stxbp1.seurat@meta.data[cells,],replicate_col = "SampleID", cell_type_col = level,label_col = "Genotype")[[1]]
    
    #filter genes with low expression
    p_samples_larger_10_counts<-70
    keep <- rowSums(DDS_full) >= 10
    DDS_full <- DDS_full[keep,]
    n_samples_min <- dim(DDS_full)[2]*(p_samples_larger_10_counts/100)
    DDS_full <- DDS_full[rowSums(DDS_full>=10)>n_samples_min,]

    metadata_full = stxbp1.seurat@meta.data[cells,]%>% group_by(SampleID) %>% 
      transmute(mean_nCount = mean(nCount_RNA), mean_nFeature = mean(nFeature_RNA), orig.ident = orig.ident, 
                mean_pct_mito =  mean(percent_mito), mean_total_mito =  mean(total_counts_mt), 
                AgeWeek, Genotype, LiveCellsPercent, Sex, TissueExtraction
      ) %>% 
      distinct()
    metadata_full<-as.data.frame(metadata_full)
    row.names(metadata_full) <- paste0(metadata_full$SampleID, ':',metadata_full$Genotype)
    metadata_full <- metadata_full[colnames(DDS_full),] ##metadata needs to be ordered as samples are in count matrix
    
    
    #Use RUV calculated with full cell type set previously
    RUVs <- read.csv(paste0('results/20220511_pseudobulkDESeqRUV_subclass/2RUVfactor',gsub(' ','',gsub('/', '', ct)),'.csv'), row.names = 1)
    metadata_full$W1 <- RUVs[row.names(metadata_full),'W_1']
    metadata_full$W2 <- RUVs[row.names(metadata_full),'W_2']
    
    #change full desing and run DESeq2
    design.formula = ~ W1+W2+Genotype
    reduced = ~ W1 + W2
    
    dds_full <- DESeqDataSetFromMatrix(DDS_full, 
                                       colData = metadata_full,
                                       design = design.formula)
    
    ### run DESEq2 
    dds_full <- DESeq(object=dds_full,parallel=F, test="LRT", reduced = reduced, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
    res_full <- results(dds_full,alpha=pthr, contrast=c("Genotype", "HZ", "WT"), pAdjustMethod = 'fdr', independentFiltering = F)
    
    ## Add results to table  
    
    topgenes <- res_full[res_full$padj<pthr & abs(res_full$log2FoldChange)>fcthr,]
    
    exp = c(celltype = ct, nCells = length(cells), nDEGs = nrow(topgenes), nUpreg =  nrow(topgenes[topgenes$log2FoldChange>0,]),
            nDownreg = nrow(topgenes[topgenes$log2FoldChange<0,]), nGtested =  nrow(res_full), round = r) 
    
    nDEGs <- rbind(nDEGs, exp)

    #save results table and DEGs
    write.csv2(nDEGs, file = paste0(PATH_OUT, level,"_results_n",n,"_iterative_DESeq2.csv"), quote = F, row.names = F)
    if(length(rownames(topgenes))>0){
      write(rownames(topgenes), paste0(PATH_OUT , level,"_results_n",n,"_DEGs_",gsub(' ','',gsub('/', '', ct)),'_', r, '.txt'), sep = '\n')
    }
  }

  write.csv2(sampleHZ,paste0(PATH_OUT,'GABA_iteration_results_n35/',gsub(' ','',gsub('/', '', ct)), "randomHZcells.txt"), quote = F, row.names = F)
  write.csv2(sampleWT,paste0(PATH_OUT,'GABA_iteration_results_n35/',gsub(' ','',gsub('/', '', ct)), "randomWTcells.txt"), quote = F, row.names = F)
    
}

nDEGs<- nDEGs[2:nrow(nDEGs),]
nDEGs$percent.degs <- 100*as.numeric(nDEGs$nDEGs)/as.numeric(nDEGs$nGtested)

## Plot results
down_degs <- nDEGs

    #get cell type colors and add to results data.frame
ct_cols <- read.csv2('data/figure_colors/class_colors.csv')[,2:3]
rownames(ct_cols) <- ct_cols$cts
unique_celltypes <- unique(down_degs$celltype)
color_df <- data.frame(celltype = unique_celltypes, color = ct_cols[unique_celltypes, 'ct_cols'])
color_df$color[color_df$celltype=='unassigned'] <- '#D2D2D2'

down_degs <- merge(down_degs, color_df, by = "celltype", all.x = TRUE)

#mean percentage of DEGs per cell type
mean_percent <- down_degs %>%
  filter(nCells == n) %>%
  group_by(celltype, color) %>%
  summarize(mean_percent = mean(percent.degs)) %>%
  arrange(mean_percent)
down_degs$celltype <- factor(down_degs$celltype, levels = mean_percent$celltype)

#plot
ggplot(down_degs, aes(y = log10(1+percent.degs), x =celltype))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(width = 0.15, alpha = 0.3, aes(color = color), size = 0.5, shape = 16)+
  geom_point(data = mean_percent, aes(y = log10(1+mean_percent), x =celltype),
             shape = 4, size=1, color="black")+
  coord_flip(ylim=c(NA,0.25))+scale_color_identity()+
  theme( axis.text = element_text(size = 5), axis.title = element_text(size = 5), 
         axis.title.y = element_blank(), panel.grid.minor = element_line(linewidth = 0.25), panel.grid.major = element_line(linewidth = 0.25),
         panel.grid.major.y = element_blank(), axis.ticks = element_line(linewidth = 0.25),
         panel.border = element_rect(linewidth = 0.25, colour = 'black', fill = NA),
         axis.line =element_blank())

ggsave(filename =  paste0(PATH_OUT,'downsamplig_',level, '.svg'), width = 1.75, height = unique(down_degs$celltype)*0.125)
