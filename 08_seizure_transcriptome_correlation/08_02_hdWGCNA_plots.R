#####
#authtor: Amparo Roig Adam
#created: 2025/02/16
#description: plot hdWGCNA trait correlation result

#input: output from 08_01_hdWGCNA.R
#output: correlation plots and GO enrichment custom plots
#####

for(ct in c('Glutamatergic neurons', 'GABAergic neurons', 'Astrocytes')){
  #load cell class specific hdWGCNA results
  load(paste0('results/20250216_hcwgcna/hdwgcna_',gsub(' ', '',ct),'_HZonly_20250226.RData'))
  
  #get module colors
  plot_cols <- unique(modules[, 'color'])
  names(plot_cols) <- str_extract(unique(modules[, 'module']) , "M\\d+")
  
  ## Plot correlation to events in ct
  mod_cor <- reshape2::melt(mt_cor$cor[ct], value.name = 'cor')
  mod_fdr <-  reshape2::melt(mt_cor$fdr[ct], value.name = 'fdr')
  
  plot <- merge(mod_cor, mod_fdr, by=c("Var1","Var2", "L1"))
  colnames(plot) <- c('event', 'module','celltype', 'cor', 'fdr')
  
  plot <- plot[plot$event%in%c('jumps', 'twitches'),]
  plot$module <-str_extract(plot$module, "M\\d+")
  plot$module <- factor(plot$module, levels = paste0("M", 1:length(unique(plot$module))))
  
  ggplot(plot[plot$celltype==ct&plot$fdr<0.05,], aes(x= cor, y = module, fill = module))+
    geom_col()+
    scale_fill_manual(values = plot_cols, guide = "none")+
    scale_y_discrete(limits = rev)+
    geom_vline(xintercept = 0, linetype = 2)+ 
    facet_wrap(.~event, ncol = 2)+theme_linedraw()+xlab('correlation')+
    theme(legend.position = 'right', text = element_text(size = 6), axis.title = element_text(size = 6),
          strip.background = element_blank(), strip.text = element_text(colour = 'black'),
          axis.title.y = element_blank(), plot.background = element_blank(), panel.background = element_blank())
  
  ggsave(paste0('results/20250216_hcwgcna/figures/',ct,'_corr.eps'), width = 80, height = 4*nrow(plot[plot$celltype==ct&plot$fdr<0.05,]),units = 'mm')
  
  if(ct %in%c('Glutamatergic neurons', 'GABAergic neurons')){
    ## Plot correlation to events in celltypes
    ct_cols <- read.csv2('data/figure_colors/class_colors.csv')[,2:3]
    plot_ct_cols <- ct_cols$ct_cols
    names(plot_ct_cols)<- ct_cols$cts
    
    cts = unique(wcgna_obj$scmap_subclass[wcgna_obj$class03==ct])
    cts = cts[cts%in%c('unassigned', 'Sst Chodl')==F]
    
    mod_cor <- reshape2::melt(mt_cor_cts$cor[cts], value.name = 'cor')
    mod_fdr <-  reshape2::melt(mt_cor_cts$fdr[cts], value.name = 'fdr')
    
    plot <- merge(mod_cor, mod_fdr, by=c("Var1","Var2", "L1"))
    colnames(plot) <- c('event', 'module','celltype', 'cor', 'fdr')
    
    plot <- plot[plot$event%in%c('jumps', 'twitches'),]
    plot$module <-str_extract(plot$module, "M\\d+")
    plot$module <- factor(plot$module, levels = paste0("M", 1:length(unique(plot$module))))
    
    ggplot(plot[plot$fdr<0.05,], aes(x= cor, y = module, fill = celltype))+
      geom_col(position = position_dodge2(preserve = "single", reverse = T))+
      scale_y_discrete(limits = rev)+scale_fill_manual(values = plot_ct_cols)+
      geom_vline(xintercept = 0, linetype = 2)+
      facet_wrap(.~event,ncol = 2)+theme_linedraw()+xlab('correlation')+
      theme(legend.position = 'none', text = element_text(size = 6), axis.title = element_text(size = 6),
            strip.background = element_blank(), strip.text = element_text(colour = 'black'),
            axis.title.y = element_blank(), plot.background = element_blank(), panel.background = element_blank())
    
    ggsave(paste0('results/20250216_hcwgcna/figures/',ct,'_corr_cts.eps'), width = 80, height = 6*length(unique(plot[plot$fdr<0.05,'module'])),units = 'mm')
  }
  ##GO enrichment custom plotting
  
  #modules to use (enriched in either one type or class)
  m = colnames(mt_cor$fdr[[ct]])[colSums(mt_cor$fdr[[ct]][c('jumps', 'twitches'),]<0.05)>0] #get M number with modules enriched in at least one event type
  
  #include any enriched in celltype
  if(ct %in%c('Glutamatergic neurons', 'GABAergic neurons')){
    cts = unique(wcgna_obj$scmap_subclass[wcgna_obj$class03==ct])
    cts = cts[cts%in%c('unassigned', 'Sst Chodl')==F]
    for(c in cts){
      ms<- colnames(mt_cor_cts$fdr[[c]])[colSums(mt_cor_cts$fdr[[c]][c('jumps', 'twitches'),]<0.05)>0]
      m <- c(m,ms)
    }}
  
  m <- str_extract(unique(m), "M\\d+")
  
  enrich_df$module <-str_extract(enrich_df$module, "M\\d+")
  enrich_df$module <- factor(enrich_df$module, levels = paste0("M", 1:length(unique(enrich_df$module))))
  
  enrich_df[enrich_df$Adjusted.P.value<0.05&enrich_df$module%in%m&enrich_df$Odds.Ratio>5,]%>%
    group_by(module) %>%
    arrange(module, Adjusted.P.value) %>%
    slice_head(n = 5) %>%  
    mutate(
      Term = str_remove(Term, "\\s*\\(GO:\\d+\\)"),  # Remove (GO:XXXXX)
      Term =  fct_reorder(Term, Adjusted.P.value, .desc = F)  # Order terms by p-value *within each module*
    ) %>%
    ggplot(aes(x=module, y = Term, group = module))+
    geom_point(aes(size = Odds.Ratio, colour = -log10(Adjusted.P.value)))+
    scale_color_stepsn(colors=rev(viridis::magma(256)))+
    scale_y_discrete(limits = rev) +
    theme_linedraw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.background = element_blank(), legend.key.size = unit(2.25, units = 'mm'),
          legend.position = 'right',
          plot.background = element_blank(), panel.background = element_blank(),
          text = element_text(size = 6), axis.title = element_blank(),
          strip.background = element_blank())
  
  plot_height<- nrow(enrich_df[enrich_df$Adjusted.P.value<0.05&enrich_df$module%in%m&enrich_df$Odds.Ratio>5,]%>%
                       group_by(module) %>%
                       arrange(module, Adjusted.P.value) %>%
                       slice_head(n = 5) )
  
  ggsave(paste0('results/20250216_hcwgcna/figures/',gsub(' ', '',gsub('/','',ct)),'_signifmodsGOs.eps'),
         , width = 60+6*length(unique(enrich_df[enrich_df$Adjusted.P.value<0.05&enrich_df$module%in%m,'module'])), 
         height = 8+3.5*plot_height,units = 'mm')
}