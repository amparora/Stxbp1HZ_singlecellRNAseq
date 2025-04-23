#####
#authtor: Amparo Roig Adam
#created: 2022/05/03
#last update: 2024/04/18
#description: DESeq2 analysis including RUV factor and normalisation - specific for cell classes and celltypes
####################

library("Libra")
library("Seurat")
library("ggplot2")
library('ggrepel')
library("dplyr")
library("DESeq2")
library("RUVSeq")

#plotting settings
theme_set(theme_minimal(base_size = 12))
theme_update(axis.line=element_line(linetype = 1, size = 0.25))

## paths and opts
write_folder = 'results/20240503_pseudobulkDESeq2_class03/'

#annotation level to calculate DEGs
level = 'class03' #'scmap_subclass' #CHANGE ACCORDINGLY TO CALCULATE DEGs class03 == cell classes; scmap_subclass == all neuronal types

#thresholds
pthr = 0.05
fcthr = 0.05

####################
#read data
stxbp1.seurat <- readRDS('results/annotationsmerged_seurat.rds')
stxbp1.seurat$scmap_subclass[is.na(stxbp1.seurat$scmap_subclass)]<- 'NA'

## Make Pseudobulk matrix from single cell counts for each cell class/neuronal type
ctcol = level
DDS <-to_pseudobulk(stxbp1.seurat@assays$RNA@counts, meta = stxbp1.seurat@meta.data, min_features = -Inf ,replicate_col = "SampleID", cell_type_col = ctcol, label_col = "Genotype")

##Run rest of the analysis for all cell classes/neuronal types
ct.list = names(DDS)[names(DDS)!='NA']
for(ct in ct.list){
    print(ct)
    ctDDS <- DDS[[ct]]

    # remove low count genes - genes with less than 10 total counts
    keep <- rowSums(ctDDS) >= 10
    ctDDS <- ctDDS[keep,]

        #at least p% of samples need to have more than 10 counts for a gene to be kept
    p_samples_larger_10_counts<-70
    n_samples_min <- dim(ctDDS)[2]*(p_samples_larger_10_counts/100)
    ctDDS <- ctDDS[rowSums(ctDDS>=10)>n_samples_min,]

    # save pseudobulk matrix
    write.csv(ctDDS, paste0(write_folder,'pseudobulkmatrix_fulldata.csv'))

    ## Format metadata for DESeq2
    metadata = stxbp1.seurat@meta.data[stxbp1.seurat@meta.data[ctcol] == ct,] %>% group_by(SampleID) %>% 
    transmute(mean_nCount = mean(nCount_RNA), mean_nFeature = mean(nFeature_RNA), orig.ident = orig.ident,
                mean_pct_mito =  mean(percent_mito), mean_total_mito =  mean(total_counts_mt), 
                AgeWeek, Genotype, LiveCellsPercent, Sex, TissueExtraction
    ) %>% 
    distinct()
    metadata<-as.data.frame(metadata)
    row.names(metadata) <- paste0(metadata$SampleID, ':',metadata$Genotype)
    metadata <- metadata[colnames(ctDDS),] ##metadata needs to be ordered as samples are in count matrix

    ## DESeq2 - run with only the factor of interest in design (Genotype)
        # set up DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(countData = ctDDS, 
                                    colData = metadata,
                                    design = ~ Genotype)
        #Run basic design DESeq2
    dds <- DESeq(object=dds ,parallel=F, test="LRT", reduced = formula(~1), useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)

        #plot and save results
    res <- results(dds,alpha=pthr, contrast=c("Genotype", "HZ", "WT"), pAdjustMethod = 'fdr', independentFiltering = F)
    topgenes <- res[res$padj<pthr & is.na(res$padj) == F,]
    topgenes<- topgenes[order(topgenes$padj),]

    print(ggplot(data.frame(res), aes(x=log2FoldChange, y=-log10(padj), color = abs(log2FoldChange)>fcthr&padj<pthr  )) +
            geom_point(show.legend = T, alpha = 0.7) + NoLegend() + scale_color_manual(values = c('darkgrey', 'red'))+
            ggtitle(label = paste0(ct, ' - DESeq2'), subtitle = paste0(nrow(res[res$padj<pthr & is.na(res$padj) == F,]), ' DEGs - ', nrow(res), 'genes tested'))+
            geom_text_repel(data= data.frame(topgenes[topgenes$log2FoldChange<0,]),show.legend = FALSE,aes(label=rownames(topgenes[topgenes$log2FoldChange<0,]))))+
            geom_text_repel(data= data.frame(topgenes[topgenes$log2FoldChange>0,]), show.legend = FALSE,aes(label=rownames(topgenes[topgenes$log2FoldChange>0,])))
    
    ggsave(paste0(write_folder,"basicdesign_volcano_",gsub(' ','',gsub('/', '', ct)),".png"),  width =6, height = 5)
    
    write.csv2(res[order(res$padj),], paste0(write_folder,"basicdesign_DEGs_",gsub(' ','',gsub('/', '', ct)),".csv"))

    ## Calculate RUV factors in cell type
    ct.set <- newSeqExpressionSet(counts(dds, normalized = F), phenoData = metadata)
    idx  <- rowSums(counts(ct.set) > 5) >= 2
    ct.set  <- ct.set[idx, ]
    par(mar = c(7, 3, 1, 1), mfrow = c(1, 3))
    plotRLE(ct.set, col = as.numeric(ct.set$Genotype)+1, outline = F, las = 2, cex.axis = 0.75)
    
    #Perform normalisation
    ct.set <- betweenLaneNormalization(ct.set, which = 'upper')
    plotRLE(ct.set, col = as.numeric(ct.set$Genotype)+1, outline = F, las = 2, cex.axis = 0.75)
    plotPCA(ct.set,col = as.numeric(ct.set$Genotype))

    #Select empirical control genes (genes with high p-val in basic design DESeq2)
    not.sig <- rownames(res)[which(res$pvalue > 0.7)]
    empirical <- rownames(ct.set)[rownames(ct.set) %in% not.sig ]
    
    #run RUVg
    ct.set <- RUVg(ct.set, empirical, k = 2)
    pData(ct.set)
    
    ## Rerun DESeq2 including RUVs
    dds_RUV <- dds
    dds_RUV$W1 <- ct.set@phenoData@data[row.names(dds_RUV@colData),'W_1']
    dds_RUV$W2 <- ct.set@phenoData@data[row.names(dds_RUV@colData), 'W_2']

    #run DESeq2 - 1 RUV
    design.formula = ~ W1+Genotype
    reduced = ~ W1
    
    dds_1RUV <- DESeqDataSetFromMatrix(ctDDS, 
                                        colData = colData(dds_RUV)[,colnames(colData(dds_RUV))!='sizeFactor'], 
                                        design = design.formula)
    dds_1RUV <- DESeq(object=dds_1RUV,parallel=F, test="LRT", full = design.formula, reduced = reduced, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
    
    res_1RUV <- results(dds_1RUV,alpha=pthr, contrast=c("Genotype", "HZ", "WT"), pAdjustMethod = 'fdr', independentFiltering = F)
    topgenes <- res_1RUV[abs(res_1RUV$log2FoldChange)>fcthr&res_1RUV$padj<pthr & is.na(res_1RUV$padj) == F,]
    
    print(
        ggplot(data.frame(res_1RUV), aes(x=log2FoldChange, y=-log10(padj), color = abs(log2FoldChange)>fcthr&padj<pthr &is.na(padj)==F)) + 
        geom_point(show.legend = T, alpha = 0.7)+NoLegend() + scale_color_manual(values = c('darkgrey', 'red'))+
        ggtitle(label = paste0(ct, ' - ctRUV_DESeq2', '\n ~ ',design(dds_1RUV)[2]), subtitle = paste0(nrow(res_1RUV[res_1RUV$padj<pthr & is.na(res_1RUV$padj) == F,]), ' DEGs - ', nrow(res_1RUV), 'genes tested'))+
        geom_text_repel(data= data.frame(topgenes[topgenes$log2FoldChange<0,]),show.legend = FALSE,aes(label=rownames(topgenes[topgenes$log2FoldChange<0,]))))+
        geom_text_repel(data= data.frame(topgenes[topgenes$log2FoldChange>0,]), show.legend = FALSE,aes(label=rownames(topgenes[topgenes$log2FoldChange>0,])))
    
    ggsave(paste0(write_folder,"ct1RUV_volcano_",gsub(' ','',gsub('/', '', ct)),".png"),  width =6, height = 5)
    
    write.csv2(res_1RUV[order(res_1RUV$padj),], paste0(write_folder,"ct1RUV_DEGs_",gsub(' ','',gsub('/', '', ct)),".csv"))
    
    #run DESeq2 - 2 RUV
    design.formula = ~ W1+W2+Genotype
    reduced = ~ W1+W2
    
    dds_2RUV <- DESeqDataSetFromMatrix(ctDDS, 
                                        colData = colData(dds_RUV)[,colnames(colData(dds_RUV))!='sizeFactor'], 
                                        design = design.formula)
    dds_2RUV <- DESeq(object=dds_2RUV,parallel=F, test="LRT", full = design.formula, reduced = reduced, useT=TRUE, minmu = 1e-6, minReplicatesForReplace = Inf)
    
    res_2RUV <- results(dds_2RUV,alpha=pthr, contrast=c("Genotype", "HZ", "WT"), pAdjustMethod = 'fdr', independentFiltering = F)
    topgenes <- res_2RUV[abs(res_2RUV$log2FoldChange)>fcthr&res_2RUV$padj<pthr & is.na(res_2RUV$padj) == F,]
    
    print(
        ggplot(data.frame(res_2RUV), aes(x=log2FoldChange, y=-log10(padj), color = abs(log2FoldChange)>fcthr&padj<pthr &is.na(padj)==F)) + 
        geom_point(show.legend = T, alpha = 0.7)+NoLegend() + scale_color_manual(values = c('darkgrey', 'red'))+
        ggtitle(label = paste0(ct, ' - ctRUV_DESeq2', '\n ~ ',design(dds_2RUV)[2]), subtitle = paste0(nrow(res_2RUV[res_2RUV$padj<pthr & is.na(res_2RUV$padj) == F,]), ' DEGs - ', nrow(res_2RUV), 'genes tested'))+
        geom_text_repel(data= data.frame(topgenes[topgenes$log2FoldChange<0,]),show.legend = FALSE,aes(label=rownames(topgenes[topgenes$log2FoldChange<0,]))))+
        geom_text_repel(data= data.frame(topgenes[topgenes$log2FoldChange>0,]), show.legend = FALSE,aes(label=rownames(topgenes[topgenes$log2FoldChange>0,])))
    
    ggsave(paste0(write_folder,"ct2RUV_volcano_",gsub(' ','',gsub('/', '', ct)),".png"),  width =6, height = 5)
    
    write.csv2(res_2RUV[order(res_2RUV$padj),], paste0(write_folder,"ct2RUV_DEGs_",gsub(' ','',gsub('/', '', ct)),".csv"))
    
    #check empirical 'control' genes used for RUV do not overlap with final result DEGs
    print(paste0("Overlap ct empirical genes with DEGs in ct"))
    print(table(empirical%in%row.names(topgenes)))

    ## Plot Relative Log Expression (RLE)
    png(paste0(write_folder,"RLE_ctRUV_",gsub(' ','',gsub('/', '', ct)),".png"), width = 800, height = 400)
    par(mfrow = c(1, 2), mar = c(7, 3, 1, 1))
    plotRLE(counts(dds_2RUV, normalized = F), outline=F,  
            col=as.numeric(dds_2RUV$Genotype)+1, 
            main = 'Raw Counts', las = 2)
    plotRLE(counts(dds_2RUV, normalized = T), outline=F, 
            col = as.numeric(dds_2RUV$Genotype)+1, 
            main = 'Normalized Counts', las = 2)
    dev.off()  
    
    ## Plot RUV PCA
    png(paste0(write_folder,"RUVk_",gsub(' ','',gsub('/', '', ct)),"_data.png"), width = 800, height = 400)
    par(mfrow = c(1, 2))
    for(k in 1:2) {
        set0 <- newSeqExpressionSet(counts(dds, normalized = F), phenoData = metadata)
        set0 <- betweenLaneNormalization(set0, which = 'upper')
        set_g <- RUVg(x = set0, cIdx = empirical, k = k)
        plotPCA(set_g, col=as.numeric(ct.set$Genotype), cex = 0.9, adj = 0.5, 
                main = paste0('with RUVg, k = ',k), 
                ylim = c(-1, 1), xlim = c(-1, 1), )
    }
    dev.off()
    
    ##Save RUV factor and results 
    write.csv(ct.set@phenoData@data[,c('SampleID', 'W_1', 'W_2')], paste0(write_folder,"2RUVfactor",gsub(' ','',gsub('/', '', ct)),".csv" ))
    
    save(dds, dds_1RUV, dds_2RUV, res, res_1RUV, res_2RUV, empirical, file = paste0(write_folder,gsub(' ','',gsub('/', '', ct)),'DESeq2.RData'))
}