#####
#authtor: PF Sullivan; adapted by Annemiek A van Berkel and Amparo Roig Adam
#created: 2022/07
#last modified: 2024/12/18
#description: process brain scRNAseq data from  experiment contrasts mice +/- vs +/+ for Stxbp1
#####

#=== initialize
library(tidyverse)
library(data.table)
library(skimr)
library(rstatix)
library(ComplexUpset)
setwd("results/otherstudies/Overlap Genes")

stxbp1.seurat <- readRDS('annotationsmerged_seurat.rds') 

#=== read full DSeq2 DEG results for all cell types 

  #for each cell type/class, need a file in cell_classes/ for : ct_background.txt; ct_DESeq2_results.csv
a <- list.files(path = "cell_classes", pattern="background", full.names = TRUE) %>%
  map_df(function(x) fread(x, header = F) %>% mutate(fname=basename(x)))
b <- a %>% 
  mutate(celltype = str_replace_all(fname, "_background.txt", ""),
         tmp = TRUE) %>% 
  select(mus_name=V1, celltype, tmp) %>% 
  pivot_wider(id_cols = mus_name, names_from = celltype, values_from = tmp, values_fill = FALSE)

# upset plot, long right tail of intersections
c <- colnames(select(b, -mus_name))
upset(b, name = "crude background (min 200)",  intersect = c, min_size = 200)
ggsave("background-upset.png")

#=== get mus_gene-to-ensgid link from her data 
d <- list.files(path = "cell_classes", pattern=".csv", full.names = TRUE) %>%
  map_df(function(x) fread(x, header = T) %>% mutate(fname=basename(x)))
e <- d %>%
  mutate(celltype = str_replace_all(fname, "_DESeq2_result.csv", ""),
         mus_name = coalesce(X, V1)) %>% 
  select(ensgid.mouse=ENSID, mus_name, celltype, baseMean:padj)
  
  #add missing ENSIDs
meta_features <- stxbp1.seurat@assays$RNA@meta.features %>%
  as.data.frame() %>%  # Ensure it is a data frame
  select(Accession, var_names)
e <- e %>%
  left_join(meta_features, by = c("mus_name" = "var_names")) %>%
  mutate(
    ensgid.mouse = coalesce(ensgid.mouse, Accession)  # Fill missing ENSID
  ) %>%
  select(-Accession) 

skim(e)

#=== finalize the gene meta-data file
f <- e %>% distinct(ensgid.mouse, mus_name)
# skim(f) 
g <- inner_join(b, f)

# add info, 163 not found
h <- read.csv("musGeneMatrix.csv") %>% 
  select(ensgid.mouse=ensgid, h0:hstr, ensgid.human, MouseHuman1to1, gene_type, description)
i <- left_join(g, h) %>% 
  mutate(musNameFound = is.na(h0)==F,
         mainGeneType = gene_type %in% c("lncRNA", "protein_coding")) %>% 
  select(ensgid.mouse, mus_name, mainGeneType, musNameFound, everything())

# review these, cannot be used
ichk <- subset(i, musNameFound == FALSE)

# the genes I'd focus on, N=15800
i %>% filter(musNameFound == TRUE & mainGeneType==TRUE) %>% summarise(n())

#=== save processed data
fwrite(i, "Stxbp1-gene-metadata.tsv", sep = "\t")
fwrite(e, "Stxbp1-results-long.tsv", sep = "\t")

rm(list=ls(pattern="^[a-z]"))



