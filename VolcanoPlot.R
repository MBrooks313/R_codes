volPlot <- function(dat, FC.nm, FDR.nm, fc, fdr, expt.nm, sym.lab = FALSE, res.path = "../data/processed/", ret=FALSE){
  
  #############################################
  # Written by MJB on Sept 19th, 2019, last modified on 
  # dat = DE results data.frame
  # FC.nm = colname for log fold-change
  # FDR.nm = colname for FDR
  # fc = fold-change threshold
  # fdr = FDR threshold
  # expt.nm = experiment name to be used in the subtitle and export file
  # sym.lab = gene symbols labels for significant genes
  # res.path = path for writing plot
  # ret = logical to return the plot
  #############################################
  
  
  library(plyr); library(dplyr)
  library(ggplot2)
  library(ggrepel)
  
  # Prep the data table
  dat <- dat %>%
    mutate(SIG = case_when(
      get(FDR.nm) < fdr & (get(FC.nm) <= -fc | get(FC.nm) >= fc) ~ 'sig', 
      TRUE ~ 'nsig')) %>% 
    mutate(FDR = -log10(get(FDR.nm)))
  
  # Volcano plot
  p <- ggplot(dat, aes_string(FC.nm, 'FDR')) +
    geom_point(aes(fill = SIG, color = SIG), alpha = 0.75, size = 1.75, pch = 21) +
    scale_fill_manual(values = c('sig' = "red", 'nsig' = "grey70")) +
    scale_color_manual(values = c('sig' = "red", 'nsig' = "grey70")) +
    xlab("Fold Change (log2)") + ylab("FDR (-log10)") +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-fc, fc), linetype = "dashed", alpha = 0.5) +
    theme_classic() +
    guides(fill = "none") +
    labs(subtitle = paste0("Volcano Plot - ", expt.nm)) +
    annotate("text", label = paste0("FDR: ", fdr),
             y = fdr + 0.205,
             x = range(dat[, FC.nm])[1] + 0.25) +
    annotate("text", label = paste0("FC: ", round(2^fc, 2)),
             y = range(dat[, FDR.nm])[2],
             x = 0)
  
  # Add significant gene names
  if (sym.lab == TRUE){
    p +
    geom_text_repel(data = filter(dat, SIG == 'sig'),
                    aes(label = external_gene_name),
                    nudge_y = 0.25, nudge_x = 1.25, 
                    min.segment.length = 0.25, 
                    direction = "both", 
                    segment.color = "grey70")
  }
  
  # Export plot
  pdf(paste0(res.path, expt.nm, "_vol.pdf"), 
      width = 8, 
      height = 6, 
      useDingbats = F)
  print(p)
  dev.off()
  
  # Return plot if needed
  if (ret == TRUE){
    return(p)
  }
  
}


# Testing variables #################
# dat.test <- et_results_grp_test$all
# FC.nm = "logFC"
# FDR.nm = "FDR" 
# fc = 1 
# fdr = 0.05
# expt.nm = "grp_test"
# sym.nm = TRUE
####################################
