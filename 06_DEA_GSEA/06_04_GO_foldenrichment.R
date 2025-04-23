#####
#authtor: Amparo Roig Adam
#created: 2024/10/15
#description: calculate fold enrichment of go terms as in Mi, Huaiyu, et al. Nature protocols (2019) and plot grouped/ungrouped GO term enrichment results

#input: aggregated ClueGO output (fuunann and GOgroups objects generated in 06_03_GO_format_plotshared.R); genes in GO terms used by ClueGO (can be retrieved from ClueGOConfiguration folder)
#output: heatmaps and barplots showing ClueGO GSEA results
#####

library('biomaRt')
library('ggplot2')
library('dplyr')
library('reshape2')
library('Seurat')

## opt
lev = 'class' #'ctypes' #change this option to run for cell classes or neuronal celltypes
groupedGO = FALSE # change T/F for aggregating GO group data or keeping one line per GO term

## Read aggregated ClueGO output (generated in 06_03_GO_format_plotshared.R)
load(paste0(lev,'_isgrouped',groupedGO, '_aggGOtables.RData'))

## Read ClueGO data
# get all ensembl genes from a go term from same source as cluego
cluegosource<-'/Users/amparoroigadam/ClueGOConfiguration/v2.5.9/ClueGOSourceFiles/Organism_Mus Musculus/Mus Musculus_GO_BiologicalProcess-EBI-UniProt-GOA-ACAP-ARAP_10.03.2023_00h00.txt.gz'
cluegosource <- read.delim(cluegosource)

#filter to GOs in celltypes
cluegosource <- cluegosource[cluegosource$Ontology_GO_BiologicalProcess.EBI.UniProt.GOA.ACAP.ARAP_10.03.2023_00h00%in%funann$ID,]
cluegosource$AllTheGenesInTheNodeENS <- NA

ensembl <- useMart("ensembl", dataset = 'mmusculus_gene_ensembl')

for(go in cluegosource$Ontology_GO_BiologicalProcess.EBI.UniProt.GOA.ACAP.ARAP_10.03.2023_00h00){
  all.associated.genes = cluegosource$AllTheGenesInTheNode[cluegosource$Ontology_GO_BiologicalProcess.EBI.UniProt.GOA.ACAP.ARAP_10.03.2023_00h00==go]
  all.associated.genes = gsub("[A-Z|]", "", all.associated.genes)
  all.associated.genes <- unique(unlist(strsplit(all.associated.genes, ":")))
  
  
  #transform entrez to entrez ensembl
  result <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id","external_gene_name"),
                  filters = "entrezgene_id",
                  values = all.associated.genes,
                  mart = ensembl)
  
  cluegosource$AllTheGenesInTheNodeENS[cluegosource$Ontology_GO_BiologicalProcess.EBI.UniProt.GOA.ACAP.ARAP_10.03.2023_00h00==go] <- paste(result$ensembl_gene_id, collapse=',')
}

## Fold enrichment function
GOFoldEnrichment = function(background = c(), degs = c(), GOgenes = c(), get.ngenes = F){
  
  #filter go genes to only those in background
  GOgenes = GOgenes[GOgenes%in%background]
  
  observed = sum(degs%in%GOgenes)
  expected = length(degs)*length(GOgenes)/length(background)
  
  FE = observed/expected
  
  if(get.ngenes==T){
    return(c('FoldEnrichment'=FE, 'expected' = expected, 'observed'= observed ))
  }else{return(FE)}
}

## Calculate Fold Enrichment of each GO in all celltypes
funann.heatFE <- matrix(nrow = length(unique(funann$ID)), ncol = length(unique(funann$celltype)))
funann.heatpv <- matrix(nrow = length(unique(funann$ID)), ncol = length(unique(funann$celltype)))
rownames(funann.heatpv)<-unique(funann$ID)
rownames(funann.heatFE)<-unique(funann$ID)
colnames(funann.heatpv)<- unique(funann$celltype)
colnames(funann.heatFE)<- unique(funann$celltype)

for(ct in colnames(funann.heatFE)){
  print(ct)
  if(ct%in%c('Sncg', 'Vip')){class = 'gaba';ct.deg=ct}else{class = 'gluta'}
  if(class=='gluta'){ct.deg<- paste0(ct,'CTX')}else{ct.deg=ct}
  if(ct == 'astro'){ctdeg='Astrocytes'}
  if(ct == 'gaba'){ctdeg='GABAergicneurons'}
  if(ct == 'gluta'){ctdeg='Glutamatergicneurons'}
  
  if(ctdeg%in%c('Astrocytes', 'GABAergicneurons', 'Glutamatergicneurons')){
    dea <- read.table(paste0('results/20240503_pseudobulkDESeq2_class03/ct2RUV_DEGs_',ctdeg,'.csv'),sep = ';', header = T )
  }
  else{
    dea <- read.table(paste0('results/20220511_pseudobulkDESeq1RUV_',class,'_subclass/ct2RUB_',class,'_subclass/ct2RUV_DEGs_',ct.deg,'.csv'), sep =';', header = T)
  }
  dea$padj <- as.numeric(gsub(",", ".", dea$padj))
  
  #get ensembl id of genes
  bkgnd <- stxbp1.seurat@assays$RNA@meta.features[dea$X,'Accession']
  degs <- stxbp1.seurat@assays$RNA@meta.features[dea$X[as.numeric(dea$padj)<0.1],'Accession']
  
  for(go in rownames(funann.heatFE)){
    go.all.genes <- unlist(strsplit(cluegosource$AllTheGenesInTheNodeENS[cluegosource$`Ontology_GO_BiologicalProcess-EBI-UniProt-GOA-ACAP-ARAP_10.03.2023_00h00`==go], split = ','))
    
    fe = GOFoldEnrichment(background = bkgnd, degs = degs, GOgenes = go.all.genes)    
    
    #write fold enrichment into table
    funann.heatFE[go,ct] <- fe
    
    #get pval table if go in funann
    if(go %in% funann$ID[funann$celltype==ct]){
      funann.heatpv[go,ct] <- funann$Term.PValue.Corrected.with.Benjamini.Hochberg[funann$ID==go & funann$celltype==ct]
    }
  }}

## Plot Fold enrichment heatmap
  #make pval of non enriched terms 1
funann.heatpv[is.na(funann.heatpv)] <- 1

    #prep for plotting
heatmap.fe <- melt(funann.heatFE, varnames = c('GOid', 'celltype'), value.name = 'fe')
heatmap.pv <- melt(funann.heatpv, varnames = c('GOid', 'celltype'), value.name = 'pval')

all.heat <- left_join(heatmap.fe, heatmap.pv)
all.heat$pvalNA <- all.heat$pval
all.heat$pvalNA[all.heat$pval==1]<-NA

    #plot
ggplot( all.heat ,aes(x = celltype, y = GOid))+
  geom_tile(aes(fill = fe, alpha = pval<=0.05 ), color = 'grey')+ 
  scale_alpha_manual(values = c(0,1))+
  scale_fill_distiller(palette = 'YlOrRd', direction = 1)+
  theme_void()+theme(axis.text.x = element_text())
ggsave(paste0('results/20241018_GOfoldenrichment/',lev,'_heatmap_fe.svg'))

ggplot( all.heat ,aes(x = celltype, y = GOid))+
  geom_tile(aes(fill = pvalNA ), color = 'grey')+
  scale_alpha_manual(values = c(0,1))+
  scale_fill_distiller(palette = 'PuRd', direction = 1)+
  theme_void()+theme(axis.text.x = element_text())
ggsave(paste0('results/20241018_GOfoldenrichment/',lev,'_heatmap_pv.svg'))

## Calculate average Fold Enrichment of each GO group in all celltypes and plot - this section also plots non grouped GOs (see opts above)
GOgroups$avgFE <- NA

if(groupedGO){
  for(group in paste0(GOgroups$GOgroup,'-',GOgroups$celltyp )){   
    grp = strsplit(group, '-')[[1]][1]
    ct = strsplit(group, '-')[[1]][2]
    
    #get GO group genes - need to have all GO IDs to get all GOs in the group
    groupGOs <- funann.long$ID[funann.long$manual.groups==grp&funann.long$celltype==ct]
    
    #average of all FE from gos in the group - DOING THIS
    GOgroups$avgFE[GOgroups$celltyp==ct &GOgroups$GOgroup==grp] <- mean(all.heat$fe[all.heat$GOid%in%groupGOs&all.heat$celltype==ct])
    
    rm(grp, ct, groupGOs)
  }
  # Prepare data for plotting
  plot_data <- melt(GOgroups, id.vars = c('celltyp', 'GOgroup', 'group.name','ndegs', 'avgFE'),  
                    measure.vars = c('nupdegs', 'ndowndegs'))
  plot_data$groupz <- paste0(plot_data$celltyp,plot_data$group.name)
  }
if(groupedGO==F){
  ## if individual go terms instead of go groups in G0groups, get FE from there and add to column avgFE for plotting supplemental
  for(ct in GOgroups$celltyp){
    for(term in GOgroups$GOid[GOgroups$celltyp==ct]){
      GOgroups$avgFE[GOgroups$GOid==term&GOgroups$celltyp==ct] <- all.heat$fe[all.heat$GOid==term&all.heat$celltype==ct]
    }}
  # Prepare data for plotting
  plot_data <- melt(GOgroups, id.vars = c('celltyp', 'GOgroup', 'group.name','ndegs', 'avgFE', 'Term'),  
                    measure.vars = c('nupdegs', 'ndowndegs'))
  plot_data$groupz <- paste0(plot_data$celltyp,plot_data$Term) 
  }

plot_data <- plot_data %>%
  group_by(celltyp) %>%
  mutate(groupz = factor(groupz, levels = unique(groupz[order(avgFE)])))

plot_data0 <- plot_data

for(ct in unique(plot_data0$celltyp)){
  print(ct)
  plot_data <-plot_data0[plot_data0$celltyp==ct,]
  p1<-ggplot(plot_data, aes(x = groupz, y =  as.numeric(value), fill = variable))+
    geom_bar(stat = 'identity') +
    scale_fill_manual(values = c("nupdegs" = "red", "ndowndegs" = "blue")) +
    labs( y = 'Number of genes', x='') +
    {if(groupedGO)ggforce::facet_col(celltyp~., scales = 'free', space = 'free')} +
    {if(groupedGO==F)ggforce::facet_col(celltyp~group.name, scales = 'free', space = 'free')} +
    coord_flip()+
    theme(legend.position = "right",  panel.grid.major.y = element_line())
 
  data <- plot_data[!duplicated(plot_data[,c('groupz', 'group.name')]),c('celltyp','groupz', 'avgFE', 'group.name')] 
  
  p2 <- plot_data[!duplicated(plot_data[,c('groupz', 'group.name')]),c('celltyp','groupz', 'avgFE', 'group.name')] %>% 
    ggplot( aes(x = groupz, y =  avgFE))+
    geom_bar(stat = 'identity') +
    labs(y = 'Fold Enrichment', x='') +
    {if(groupedGO)ggforce::facet_col(celltyp~., scales = 'free', space = 'free')} +
    {if(groupedGO==F)ggforce::facet_col(celltyp~group.name, scales = 'free', space = 'free')} +
    coord_flip()+
    theme(legend.position = "none", panel.grid.major.y = element_line())
 
  pdf(paste0('results/20241018_GOfoldenrichment/',lev,'_isgrouped',groupedGO,'suppl_all_GOs2.pdf'), width = 20, height = nrow(plot_data)*0.15+length(unique(plot_data$group.name))*0.15)
  print(p1+p2)
  dev.off()
}