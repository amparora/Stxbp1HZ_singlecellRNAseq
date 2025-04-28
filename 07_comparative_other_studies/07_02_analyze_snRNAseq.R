#####
#authtor: PF Sullivan; adapted by Annemiek A van Berkel and Amparo Roig Adam
#created: 2022/07
#last modified: 2022/11/23
#description: analyze brain scRNAseq data from  experiment contrasts mice +/- vs +/+ for Stxbp1. Needs functions in 07_03_analysis.code.R
#####

#=== initialize
library(tidyverse)
library(data.table)
library(skimr)
library(rstatix)
library(ComplexUpset)
library(readxl)
library(rio)
setwd("results/otherstudies/Overlap Genes")

gene <- fread("Stxbp1-gene-metadata.tsv")
dseq <- fread("Stxbp1-results-long.tsv")
    ##remove OPC, Oligodendr and micrglia 
gene <-gene[,!c('OPCs', 'Oligodendrocytes', 'Microglia')]
dseq <- dseq[dseq$celltype%in%c('OPCs', 'Oligodendrocytes', 'Microglia')==F,]

#first complete genesets by merging from the different groups
genesets <- list.files(path = "genesets", pattern="*.tsv", full.names = TRUE) %>%
  map_df(function(x) fread(x, header = T))

  #Manually reduce granularity of celltypes from Pfisterer for more fair comparison - merge degs from subtypes starting with same name
for (cty in unique(genesets$geneset[genesets$group=='Pfisterer'])){
   genesets$geneset[genesets$geneset==cty] <- str_split(cty, '_')[[1]][1]
 }
 genesets$geneset[genesets$geneset=='L5']<- 'L5_6'
  
  #fix group/subgroups names
genesets$subgroup[genesets$geneset%in%c("Cid_deep_sig.txt","Cid_sup_sig.txt")] <- 'Kainate'
genesets$subgroup[genesets$geneset%in%c("Dugger_L24_sig.txt","Dugger_L56_sig.txt","Dugger_lge_sig.txt","Dugger_radial_glia_sig.txt" ,"Dugger_spn_sig.txt",        
                                       "Dugger_sst_sig.txt","Dugger_vip_sig.txt")] <- 'HNRNPU'
genesets$subgroup[genesets$geneset%in%c("Hawkins_sig.txt")] <- 'SCN1A'
genesets$subgroup[genesets$geneset%in%c("Kim_sig.txt")] <- 'Slc6a20a'
genesets$subgroup[genesets$geneset%in%c("Lee_sig.txt")] <- 'Shank2'
genesets$subgroup[genesets$geneset%in%c("Renthal_sig.txt")] <- 'Mecp2'

#=== save processed data
fwrite(genesets, "genesets.tsv", sep = "\t")
#genesets = fread('genesets_mouseStudies.tsv')

#=== analysis :: main cell types :: differential expression

a <- gene %>% 
  filter(gene_type == "protein_coding") %>% 
  inner_join(dseq)

#=== now we have to convert to human
g <- a %>% 
  filter(MouseHuman1to1 == TRUE, 
         ensgid.human != "") %>% 
  mutate(down = log2FoldChange < 0, 
         sig10 = padj < 0.10) %>% 
  select(ensgid=ensgid.human, mus_name, celltype, padj, down, sig10)
# g %>% group_by(celltype, ensgid) %>% filter(n()>1)  # no dups

freq_table(g, celltype, sig10) %>% filter(sig10 == TRUE) %>%  arrange(desc(n)) %>% select(-prop)
freq_table(g, sig10, down, na.rm = F)

### GSA (all have background_gene == T)
alloverlap<- c()
for(ct in unique(g$celltype)){
  h <- g %>% filter(celltype == ct)
  i <- GSA(h, sig10)
  GSA2excel(i, paste0("modif-gsa-",ct,"-sig10-up-down.xlsx"))
  fwrite(i, paste0("modif-gsa-", ct, "-sig10-merge.txt"), sep = "\t")
  i$celltype = ct
    #calc fold enrich
  i$fold.enrichment <- i$overlap.TestVar.geneset/((i$genes.TestVar.true*i$genes.in.geneset)/i$genes.in.backround)
  
  alloverlap<- rbind(alloverlap, i)
} 

h <- g %>% filter(celltype == "Astrocytes")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-astrocytes-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-astrocytes-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "Glutamatergicneurons")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-glutamatergicneurons-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-glutamatergicneurons-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "GABAergicneurons")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-GABAergicneurons-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-GABAergicneurons-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "Microglia")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-Microglia-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-Microglia-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "Vip")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-Vip-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-Vip-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "Sncg")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-Sncg-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-Sncg-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "L23ITCTX")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-L23ITCTX-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-L23ITCTX-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "L45ITCTX")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-L45ITCTX-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-L45ITCTX-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "L6ITCTX")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-L6ITCTX-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-L6ITCTX-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "L56NPCTX")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-L56NPCTX-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-L56NPCTX-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "L6CTCTX")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-L6CTCTX-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-L6CTCTX-sig10-merge.txt", sep = "\t")

h <- g %>% filter(celltype == "L6bCTX")
i <- GSA(h, sig10)
GSA2excel(i, "modif-gsa-L6bCTX-sig10-up-down.xlsx")
fwrite(i, "modif-gsa-L6bCTX-sig10-merge.txt", sep = "\t")

#now merge and pivot
j <- list.files(pattern=glob2rx("*modif*merge*"), full.names = TRUE) %>%####modified to get only the modif files
  map_df(function(x) fread(x, header = T) %>% mutate(fname=basename(x)))
fwrite(f, "modif-gsa-overlapgenes-allgroup-bonfgroup.txt", sep = "\t")


ja <- j %>% 
  filter(group == "mouseStudies")
jb <- j %>% 
  filter(group == "Proteomics")
jc <- j %>% 
  filter(group == "Pfisterer")
jd <- j %>% 
  filter(group == "postmortem")

#make pivot and save
ka <- ja %>% 
  mutate(cell_type = str_replace_all(fname, "-sig10-merge.txt", ""),
         tmp = TRUE) %>% 
  select(geneset, cell_type, P.bonf.group) %>% 
  pivot_wider(id_cols = geneset, names_from = cell_type, values_from = P.bonf.group, values_fill = FALSE)
fwrite(ka, "modif-gsa-overlapgenes-mouseStudies-bonfgroup.txt", sep = "\t")

kb <- jb %>% 
  mutate(cell_type = str_replace_all(fname, "-sig10-merge.txt", ""),
         tmp = TRUE) %>% 
  select(geneset, cell_type, P.bonf.group) %>% 
  pivot_wider(id_cols = geneset, names_from = cell_type, values_from = P.bonf.group, values_fill = FALSE)
fwrite(kb, "modif-gsa-overlapgenes-proteomics-bonfgroup.txt", sep = "\t")

kc <- jc %>% 
  mutate(cell_type = str_replace_all(fname, "-sig10-merge.txt", ""),
         tmp = TRUE) %>% 
  select(geneset, cell_type, P.bonf.group) %>% 
  pivot_wider(id_cols = geneset, names_from = cell_type, values_from = P.bonf.group, values_fill = FALSE)
fwrite(kc, "modif-gsa-overlapgenes-pfisterer-bonfgroup.txt", sep = "\t")

kd <- jd %>% 
  mutate(cell_type = str_replace_all(fname, "-sig10-merge.txt", ""),
         tmp = TRUE) %>% 
  select(geneset, cell_type, P.bonf.group) %>% 
  pivot_wider(id_cols = geneset, names_from = cell_type, values_from = P.bonf.group, values_fill = FALSE)
fwrite(kd, "modif-gsa-overlapgenes-postmortem-bonfgroup.txt", sep = "\t")

#all
k <- j %>% 
  mutate(cell_type = str_replace_all(str_replace_all(fname, "-sig10-merge.txt", ""),"modif-gsa-",""),
         tmp = TRUE) %>% 
  select(group, geneset, cell_type, overlap.TestVar.geneset, P.bonf05.group) %>% 
  pivot_wider(id_cols = c(group, geneset, P.bonf05.group), names_from = cell_type, values_from = c(overlap.TestVar.geneset), values_fill = FALSE)
fwrite(k, "modif-gsa-overlapgenes-allgroup-bonfgroup.txt", sep = "\t")

# Plot final figure
alloverlap$class = 'gluta cts'
alloverlap$class[alloverlap$celltype%in%c('Vip', 'Sncg')]<- 'gaba cts'
alloverlap$class[alloverlap$celltype%in%c('Astrocytes', 'GABAergicneurons', 'Glutamatergicneurons', 'Microglia', 'OPCs', 'Oligodendrocytes')]<- 'class'

alloverlap%>% 
  mutate(fold.enrichment = if_else(P.bonf.group>0.05, NA_real_, fold.enrichment)) %>%
  group_by(geneset) %>%  # Group by geneset
  filter(!all(is.na(fold.enrichment))) %>%  
  ungroup() %>%
ggplot(aes(x = geneset, y = celltype, fill = fold.enrichment))+#, alpha = is.na(value)
  geom_tile(color = 'white')+facet_grid(class~group+subgroup, scales = 'free', space = 'free')+
  theme(axis.text.x = element_text(angle = 90, h = 1))+
  scale_fill_viridis_b(na.value = 'lightgrey')+theme(panel.background = element_blank())

filename = 'genesets_foldenrich_heatmap'
ggsave(paste0('results/',filename,'.png'))
ggsave(paste0('results/',filename,'.svg'))


rm(list=ls(pattern="^[a-z]"))



