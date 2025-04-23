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

save(funann, funann.long, file = paste0(lev, '_aggGOtables.RData'))


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
