#####
#authtor: Amparo Roig Adam
#created: 2024/07/10
#description: Format ClueGO output table for custom plotting. Plots genes in shared GOterms between neuronal classes.
#####
library('ggplot2')
library('dplyr')
library('reshape2')
library('Seurat')

theme_set(theme_minimal(base_size = 12))
theme_update(axis.line=element_line(linetype = 1, size = 0.25), panel.border = element_rect(colour = "black", size=1, fill = NA))

## opt
lev = 'class' #'ctypes' #change this option to run for cell classes or neuronal celltypes
groupedGO = FALSE # change T/F for aggregating GO group data or keeping one line per GO term

#read seurat dataset
stxbp1.seurat <- readRDS('results/20220510_seurat_ctannotated.rds')

## Read and aggregate ClueGO outputs
funann <- c()
funann.long <- c() #make a long version with repeated terms in shared groups

#cell types to use according to lev
if(lev == 'class'){
  cts = c('gaba', 'gluta', 'astro')
}
if(lev == 'ctypes'){
  cts = c('Sncg', 'Vip', 'L6IT', 'L6CT', 'L56NP', 'L5IT', 'L45IT', 'L23IT')
}

for(ct in cts){ 

  if(lev == 'class'){
    ct.funann <- read.table(paste0('results/20240816_ClueGO_final/',ct,'Kappa_manualgroup.csv'), sep = ';', header = T)
    ct.funann <- cbind(ct.funann, celltype = ct)
  }
  if(lev =='ctypes'){
    if(ct%in%c('Sncg', 'Vip')){class = 'gaba'}else{class = 'gluta'}
    ct.funann <- read.table(paste0('results/20240913_ClueGO_cts_final/',ct,'/',ct,' NodeAttributeTable.txt'), sep = '\t', header = T)
    ct.funann <- cbind(ct.funann, celltype = ct, class = class)
  }

  print(nrow(ct.funann))
  
  funann <- rbind(funann, ct.funann)
  
  #make multiple lines for terms in more than one group
  for(rw in c(1:nrow(ct.funann))){
    if(length(unlist(strsplit(ct.funann$GOGroups[rw], ', ')))>1){
      #print(ct.funann$Term[rw])
      for(grp in unlist(strsplit(ct.funann$GOGroups[rw], ', '))){
        grp <- grp%>%
          gsub('\\[','', .)%>%
          gsub('\\]', '', .)
        #join each new line individualy with the individual group name
        grp.line <- ct.funann[rw,]
        grp.line$GOGroups<- paste0('[', grp, ']')
        funann.long <- rbind(funann.long, grp.line)
      }}
    else{
      funann.long <- rbind(funann.long, ct.funann[rw,])}
    }}

# Format output for grouped or not grouped plotting (see opts above)
GOgroups <- c()

if(lev =='ctypes'){
  funann.long$manual.groups <- funann.long$GOGroups #only for celltypes 
}

if(groupedGO==F){
  grouping = unique(paste(sep='/', funann.long$celltype,funann.long$manual.groups, funann.long$ID, funann.long$Term))
}
if(groupedGO){
  grouping = unique(paste(sep='/', funann.long$celltype,funann.long$manual.groups))
}

for (grpct in grouping){ 

  ct = strsplit(grpct, '/')[[1]][1]
  grp = strsplit(grpct, '/')[[1]][2]
  
  grp.degs = funann.long$Associated.Genes.Found[funann.long$manual.groups==grp&funann.long$celltype==ct ] 
  
    #include term ID and name if not grouping GO terms
  if(groupedGO==F){
    GOid = strsplit(grpct, '/')[[1]][3]
    term = strsplit(grpct, '/')[[1]][4]
    grp.degs = funann.long$Associated.Genes.Found[funann.long$manual.groups==grp&funann.long$celltype==ct&funann.long$Term==term ] 
  }
  
  if(ct%in%c('Sncg', 'Vip')){class = 'gaba'}else{class = 'gluta'}
  
  grp.degs = unique(unlist(lapply(grp.degs, function(input_str) {
    # Remove the square brackets
    cleaned_str <- gsub("\\[|\\]", "", input_str)
    # Split the string by commas and trim any extra whitespace
    elements <- trimws(strsplit(cleaned_str, ",")[[1]])
    return(elements)
  })))
  
  ndegs = length(grp.degs)
  
  grp.all.associated.genes = funann.long$All.Associated.Genes[funann.long$manual.groups==grp&funann.long$celltype==ct]
  grp.all.associated.genes = unique(unlist(lapply(grp.all.associated.genes, function(input_str) {
    # Remove the square brackets
    cleaned_str <- gsub("\\[|\\]", "", input_str)
    # Split the string by commas and trim any extra whitespace
    elements <- trimws(strsplit(cleaned_str, ",")[[1]])
    return(elements)
  })))
  
  #count up and downregulated DEGs
  if(ct == 'astro'){
    ctdeg='Astrocytes'
  }
  if(ct == 'gaba'){
    ctdeg='GABAergicneurons'
  }
  if(ct == 'gluta'){
    ctdeg='Glutamatergicneurons'
  }
  
    #for classes
  if(lev=='class'){
    ct.degs <- read.table(paste0('results/20240503_pseudobulkDESeq2_class03/ct2RUV_DEGs_',ctdeg,'.csv'), sep = ';', header = T )
      #get group name
    grp.name = funann.long$Term[funann.long$manual.groups==grp&funann.long$celltype==ct&funann.long$manual.groups.name==T][1]
  }
  
    #for cell types
  if(lev =='ctypes'){
    if(class=='gluta'){ct.deg<- paste0(ct,'CTX')}else{ct.deg=ct}
    ct.degs <- read.table(paste0('results/20220511_pseudobulkDESeq1RUV_',class,'_subclass/ct2RUB_',class,'_subclass/ct2RUV_DEGs_', ct.deg,'.csv'), sep = ';', header = T, dec = ',')
      # get GO with lowest pval as name of group 
    group00_data <- funann.long[funann.long$manual.groups == grp&funann.long$celltype==ct,]
    grp.name <- group00_data$Term[order(group00_data$Term.PValue.Corrected.with.Benjamini.Hochberg, decreasing = F)][1]
    rm(group00_data)
    }

  ct.degs$padj <-  as.numeric(gsub(",", ".", ct.degs$padj))
  ct.degs$log2FoldChange <-  as.numeric(gsub(",", ".", ct.degs$log2FoldChange))
  
  updegs <- ct.degs$X[ct.degs$log2FoldChange>0&ct.degs$padj<0.1]
  downdegs <- ct.degs$X[ct.degs$log2FoldChange<0&ct.degs$padj<0.1]
  
  grp.updegs <- sum(grp.degs%in%updegs==T)
  grp.downdegs <- sum(grp.degs%in%downdegs==T)
  

  grouped = c(celltyp = ct,GOgroup = grp, group.name = grp.name, grp.degs=paste(grp.degs, collapse = ", "), 
              ndegs = ndegs, nupdegs = grp.updegs, ndowndegs = 0-grp.downdegs, grp.all.associated.genes = paste(grp.all.associated.genes, collapse = ", "))
  
  if(groupedGO==F){       #optional include; GOid = GOid ,Term = termfor ungruped GOs
    grouped = c(GOid = GOid ,Term = term,celltyp = ct,GOgroup = grp, group.name = grp.name, grp.degs=paste(grp.degs, collapse = ", "), 
                ndegs = ndegs, nupdegs = grp.updegs, ndowndegs = 0-grp.downdegs, grp.all.associated.genes = paste(grp.all.associated.genes, collapse = ", "))
    rm(GOid, term)
  }
  GOgroups <- rbind(GOgroups, grouped)
  rm(ct, grp, grp.name, grp.degs, ndegs, grp.updegs, grp.downdegs, grp.all.associated.genes)
}
GOgroups = as.data.frame(GOgroups)

save(funann, funann.long, GOgroups, lev, groupedGO, file = paste0(lev,'_isgrouped',groupedGO, '_aggGOtables.RData'))

## Plot Shared genes in neuronal classes
#get number of shared genes
sharedGOs <- c()   
for(GOname in funann.long$Term){ 

  if(sum(funann$celltype[funann$Term==GOname]%in%c('gaba', 'gluta'))==2){ #get only neuronal classes
    cts = funann$celltype[funann$Term==GOname]

    genes.gaba = funann$Associated.Genes.Found[funann$Term==GOname & funann$celltype=='gaba']
    genes.gaba = unique(unlist(lapply(genes.gaba, function(input_str) {
      # Remove the square brackets
      cleaned_str <- gsub("\\[|\\]", "", input_str)
      # Split the string by commas and trim any extra whitespace
      elements <- trimws(strsplit(cleaned_str, ",")[[1]])
      return(elements)
    })))
    
    genes.gluta = funann$Associated.Genes.Found[funann$Term==GOname & funann$celltype=='gluta']
    genes.gluta = unique(unlist(lapply(genes.gluta, function(input_str) {
      # Remove the square brackets
      cleaned_str <- gsub("\\[|\\]", "", input_str)
      # Split the string by commas and trim any extra whitespace
      elements <- trimws(strsplit(cleaned_str, ",")[[1]])
      return(elements)
    })))
    
    #count genes per cell class
    ngaba.onlyneurons = sum(genes.gaba%in%genes.gluta==F)
    ngluta.onlyneurons = sum(genes.gluta%in%genes.gaba==F)
    
    ngaba.gluta.only = length(intersect(genes.gaba, genes.gluta))
    
    shared = c(GOname=GOname, gaba = ngaba.onlyneurons, gluta = ngluta.onlyneurons, neuron.shared = ngaba.gluta.only)

    sharedGOs <- rbind(sharedGOs, shared)
  }
}

sharedGOs = unique(as.data.frame(sharedGOs))
sharedGOsall = unique(as.data.frame(sharedGOsall))

write.csv(sharedGOs,'results/20240816_ClueGO_final/shared_gos.csv')

# Prepare data for plotting
plot_data <- melt(sharedGOs, id.vars = c('GOname'), variable.name = 'celltype', value.name = 'ngenes')

# Create the barplot
ggplot(plot_data, aes(x = GOname, y = as.numeric(ngenes), fill = celltype)) +
  geom_bar(stat = 'identity') +
  labs(title = "Number of Shared Genes",
       y = "Number of Genes") +
  coord_flip()+scale_fill_manual(values = c('#E4629B', '#49BF51', 'lightblue')) 
ggsave('results/20240816_ClueGO_final/shared_GOs.svg', width = 5, height = nrow(plot_data)*0.025, limitsize = F)
