#####
#authtor: Amparo Roig Adam - adapted from hdWGCNA vignettes
#created: 2024/10/01
#last modified: 2025/02/16
#description: use hdWGCNA to find gene modules and correlation analysis with recent seizure (twitches and jumps) history

#input: seurat object with cell annotations; table with manual annotation of seizere events in Hz mice
#output: RData file with results from WGCNA, module eigengenes, module-trait correlation and GO enrichment in modules 
#####

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(enrichR)
theme_set(theme_cowplot())
set.seed(12345)

## Load Seurat object
seurat_obj <- readRDS('results/annotationsmerged_seurat.rds')

## Read manual annotation of videos 
annot_events <-read.table('results/20230314_videoannot/annots_table.csv', header = T, sep = ',')

# make summary table with total counts per event type 
sum_annot_events <- annot_events %>%
  group_by(SampleID, Genotype)%>%
  summarise(twitches = sum(event%in%c('twitch','twitch ')), jumps = sum(event%in%c('jump', 'jump ')),
            sleep = sum(event%in%c('sleep', 'sleep ')))
rm(annot_events)

#add to seurat metadata
seurat_obj <- SetIdent(seurat_obj, value = seurat_obj$orig.ident)
seurat_obj$twitches = NA
seurat_obj$jumps = NA

for(smp in unique(sum_annot_events$SampleID)){
  seurat_obj$twitches[seurat_obj$SampleID==smp] = sum_annot_events$twitches[sum_annot_events$SampleID==smp]
  seurat_obj$jumps[seurat_obj$SampleID==smp] = sum_annot_events$jumps[sum_annot_events$SampleID==smp]
}
seurat_obj$totevent <- seurat_obj$twitches+seurat_obj$jumps

## Prep data
    #get only cell types of interest and run recommended seurat preprocessing
seurat_obj<- seurat_obj[,seurat_obj$class03%in%c('Astrocytes', 'GABAergic neurons', 'Glutamatergic neurons')&
                          seurat_obj$Genotype=='HZ']

DimPlot(seurat_obj, group.by='class03', label=TRUE) +
  umap_theme()+ NoLegend()

    #setup for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.1, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Stxbp1-WGCNA" # the name of the hdWGCNA experiment
)

    #standard seurat preprocess
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(object = seurat_obj)
seurat_obj <- ScaleData(object = seurat_obj)
seurat_obj <- RunPCA(object = seurat_obj)
seurat_obj <- RunUMAP(object = seurat_obj, dims = 1:30, reduction.name = 'seurat_umap')
DimPlot(seurat_obj, group.by='class03',reduction = 'seurat_umap', label=TRUE, shuffle = T) +
  umap_theme()+ NoLegend()

## Construct metacells  in each group

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("class03", "SampleID"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', 
  k = 15, # number of cells to aggregate
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'class03', # set the Idents of the metacell seurat object,
  min_cells = 20
)

table (seurat_obj@misc$`Stxbp1-WGCNA`$wgcna_metacell_obj$class03) #check n metacells per cell class

# normalize metacell expression matrix: 
seurat_obj <- NormalizeMetacells(seurat_obj)

## Run WGCNA per cell class
for(ct in unique(seurat_obj$class03)){
      #set up exp matrix
  wcgna_obj <- SetDatExpr(
    seurat_obj,
    group_name = ct, # the name of the group of interest in the group.by column
    group.by='class03', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'data')

    # Test different soft powers:
  wcgna_obj <- TestSoftPowers(
  wcgna_obj,
  networkType = 'signed' )

    # plot the results:
  plot_list <- PlotSoftPowers(wcgna_obj)
  print(wrap_plots(plot_list, ncol=2))

    # construct co-expression network:
  wcgna_obj <- ConstructNetwork(
  wcgna_obj,
  tom_name = paste0('stxbp1-',ct,'ModsHZ'),# name of the topological overlap matrix written to disk
  tom_outdir = 'results/20250216_hcwgcna/',
  overwrite_tom = T )

  print(PlotDendrogram(wcgna_obj, main=paste0('Stxbp1 hdWGCNA ',ct,' - HZ')))

    # compute Module Eigengenes in all cells
  wcgna_obj <- ModuleEigengenes(
    wcgna_obj,
    group.by.vars="SampleID",
    exclude_grey = TRUE
    )

    # harmonized module eigengenes:
  hMEs <- GetMEs(wcgna_obj)
    # module eigengenes:
  MEs <- GetMEs(wcgna_obj, harmonized=FALSE)

    # compute eigengene-based connectivity (kME): #########. NO
  wcgna_obj <- ModuleConnectivity(
      wcgna_obj,
      group.by = 'class03', 
      group_name = ct, # the name of the group of interest in the group.by column
    )
    # rename the modules
  wcgna_obj <- ResetModuleNames(
      wcgna_obj,
      new_name ="M"
    )
    # plot genes ranked by kME for each module
  p <- PlotKMEs(wcgna_obj, ncol=3, text_size = 2.5)
  p
  ggsave(paste0('results/20250216_hcwgcna/hdwgcna_',gsub(' ','', gsub('/','',ct)),'_HZonly_kMEs.png'), p)

    # get the module assignment table:
  modules <- GetModules(wcgna_obj) %>% subset(module != 'grey')

    ## TRAIT CORRELATION ANALYSIS
        # list of traits to correlate
  wcgna_obj$class03 <- as.character(wcgna_obj$class03)
  wcgna_obj$scmap_subclass <- as.character(wcgna_obj$scmap_subclass)

  cur_traits <- c('twitches', 'jumps', 'totevent')

        # correlation of traits to MEs
  wcgna_obj <- ModuleTraitCorrelation(
  wcgna_obj,features = 'hMEs',
  traits = cur_traits,
  group.by='class03')
    # get the module-correlation results
  mt_cor <- GetModuleTraitCorrelation(wcgna_obj)

  if(ct%in%c('Glutamatergic neurons', 'GABAergic neurons')){
    # Correlation to MEs in neuronal types
    cur_traits <- c('twitches', 'jumps', 'totevent')

    wcgna_obj$scmap_subclass[is.na(wcgna_obj$scmap_subclass)] <- 'unassigned'
    wcgna_obj <- ModuleTraitCorrelation(
    wcgna_obj, features = 'hMEs',
    traits = cur_traits,
    group.by='scmap_subclass'
    )
    # get the mt-correlation results
    mt_cor_cts <- GetModuleTraitCorrelation(wcgna_obj)
  }
  ## Module GO term enrichment
    # define the enrichr databases to test
  dbs <- c('GO_Biological_Process_2023')
  # perform enrichment tests
  wcgna_obj <- RunEnrichr(
    wcgna_obj,
    dbs=dbs,
    max_genes = Inf )
  # retrieve the output table
  enrich_df <- GetEnrichrTable(wcgna_obj)

  #save all data
  save(wcgna_obj, MEs, hMEs, mt_cor,modules, enrich_df,mt_cor_cts,
     file = paste0('results/20250216_hcwgcna/hdwgcna_',ct,'_HZonly.RData'))
}
