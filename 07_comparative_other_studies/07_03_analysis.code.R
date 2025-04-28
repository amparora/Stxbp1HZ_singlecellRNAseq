#=== 03.analysis.code.NEW.R by PF Sullivan, UNC/KI, 07/2022
# function for hypergeometric gene set analysis
# note that you have to do the following below ::  ### CHANGE to fit your environment
# update 7/2022 :: reorganized gene sets, added test var, added Bonf and FDR for all and by group
#
# NB - be extremely cautious if LD is involved, violates assumptions of this "counting" approach
#
# sampling without replacement from an urn with white and black balls eg:
# RNA-seq:  N=15000 genes, R=3000 diff expressed
# ChIP-seq: n=400 genes with peak
# observe   x=100 genes that are diff expressed and ChIP-seq
# goal: prob of >= 100 white balls when pick 400 balls from urn 15000 balls & 3000 white ?
#
# N = Ngene = # genes in geneset_background (ie brain expressed)
# R = Rflag = # genes implicated in the test set
# n = nVar  = genes in some geneset
# x = xVar  = overlap of R and n
# Phyper    = phyper(x-1, R, N-R, n, lower.tail=FALSE)


#=== gene set analyses
# USAGE :: mygsa <- GSA(mydata, logicalFlag)
# mydata (aa below) is a dataframe with ensgid. All rows in mydata from the gene set background
# mydata contains a logical column (bb below), TRUE=query genes for gene set analysis, FALSE=otherwise
# cc is the minimum gene overlap to report (default 0), does not affect counting

library(tidyverse)
library(data.table)
library(rstatix)
library(readxl)
library(writexl)

GSA <- function(aa, bb, cc = 0) {
  require(tidyverse)
  require(data.table)
  require(rstatix)
  message("GSA :: hypergeometric gene set analysis by pfs 07/2022")
  message("GSA :: recommendation, ensure clean input df (no missings, all ensgid are in geneMatrix)")
  
  # sanitize input (drop missing and duplicates)
  aa <- aa %>% select(ensgid, {{bb}})
  TestVar <-  colnames(aa)[2]
  zz <- aa %>% 
    filter(is.na(ensgid)==FALSE, 
           str_sub(ensgid, 1, 4) == "ENSG",
           is.na({{bb}})==FALSE
          , {{bb}} %in% c(TRUE, FALSE)) %>% 
    mutate(Test = {{bb}}) %>% 
    select(ensgid, Test)
  zz <- unique(zz, by = "ensgid")
  zzsum <- sum(zz$Test)
  message(paste("GSA ::", nrow(aa), "rows input,", nrow(zz), "after cleaning"))
  message(paste("GSA :: for", TestVar, "fraction TRUE", zzsum/nrow(zz)))

  # read gene sets (are unique)
  yy <- fread(paste0("genesets.tsv"))  ### CHANGE to fit your environment
  
  # left join to preserve all genesets
  # assume :: all rows in input are valid (no missing or blank values) :: if not, nrow(zz) will be too large
  # note :: 1-2% of geneset ensgid will not be in input (eg input is PC, genesets have lncRNA)
  # compute worst case Phyper too
  Genesets <- merge(yy, zz, by.x = "ensgid", by.y = "ensgid", all.x = TRUE, all.y = FALSE)
  
  # counting
  Test <- Genesets %>%
    mutate(Test2 = Test) %>% 
    replace_na(list(Test2=FALSE)) %>% 
    group_by(group, subgroup, geneset) %>% 
    summarise(nVar  = sum(is.na(Test)==FALSE),
              xVar  = sum(Test, na.rm = T),
              nVarX = n(),
              xVarX = sum(Test2)) %>%
    ungroup %>%
    filter(xVar >= cc) %>%   # minimum size
    mutate(TestVar = TestVar,
           Ngene = nrow(zz),   # background is intended post-QC, not nrow(aa)
           Rflag = zzsum,
           Phyper = phyper(xVar-1, Rflag, Ngene-Rflag, nVar, lower.tail = FALSE),
           phCheck = phyper(xVarX-1, Rflag, nrow(aa)-Rflag, nVarX, lower.tail = FALSE))
  message(paste("GSA ::", nrow(Test), "gene set tests, minimum overlap", cc, "genes"))
  
  # add signif corrections
  xx <- Test %>% 
    adjust_pvalue(p.col = "Phyper", output.col = "P.bonf.all", method = "bonferroni") %>% 
    adjust_pvalue(p.col = "Phyper", output.col = "P.fdr.all", method = "fdr") %>% 
    mutate(tests.all = nrow(Test))
  ww <- xx %>% 
    group_by(group) %>% 
    summarise(tests.group = n()) %>% 
    ungroup
  Results <- xx %>% 
    group_by(group) %>% 
    adjust_pvalue(p.col = "Phyper", output.col = "P.bonf.group", method = "bonferroni") %>% 
    adjust_pvalue(p.col = "Phyper", output.col = "P.fdr.group", method = "fdr") %>% 
    ungroup %>% 
    left_join(ww) %>% 
    rename(P.hyper = Phyper) %>% 
    mutate(P.bonf05.all = P.bonf.all < 0.05,
           P.bonf05.group = P.bonf.group < 0.05) %>% 
    select(TestVar, group, subgroup, geneset,
           genes.in.backround=Ngene, genes.TestVar.true=Rflag, genes.in.geneset=nVar, overlap.TestVar.geneset=xVar,
           P.hyper, 
           P.bonf05.all, P.bonf.all, P.fdr.all, tests.all, 
           P.bonf05.group, P.bonf.group, P.fdr.group, tests.group, 
           phCheck, xVarX, nVarX) %>% 
    arrange(P.hyper)
  
  return(Results)
  rm(zz,zzsum,yy,xx,ww,TestVar,Genesets,Test)
}


#=== GSA2excel
# USAGE ::  GSA2excel(mygsa, "test.xlsx")
# write GSA output to excel spreadsheet, readme plus one worksheet per group
# aa :: df containing GSA output :: eg mygsa for mygsa <- GSA(mydata, someFlag)
# bb :: excel file name in quotes
GSA2excel <- function(aa, bb) {
  require(tidyverse)
  require(data.table)
  require(writexl)
  require(readxl)
  message("GSAexcel :: results to excel by pfs 07/2022")
  readme <- read_xlsx("00metadata-df.xlsx")    ### CHANGE to fit your environment
  # maka a list of the df to write
  zz <- c(list(readme), split(aa, f = aa$group))  
  write_xlsx(zz, path = as.character(bb), col_names = TRUE, format_headers = TRUE)
  rm(readme,zz)
}

#=== Fold enrichment
FoldEnrichment = function(background = c(), degs = c(), geneset = c(), get.ngenes = F){
  
  #filter go genes to only those in background
  geneset = geneset[geneset%in%background]
  
  observed = sum(degs%in%geneset)
  expected = length(degs)*length(geneset)/length(background)
  
  FE = observed/expected
  
  if(get.ngenes==T){
    return(c('FoldEnrichment'=FE, 'expected' = expected, 'observed'= observed ))
  }else{return(FE)}
}

