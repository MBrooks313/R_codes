pca_simple <- function(dat, gp, pc.x=1, pc.y=2, clr, expt.nm, res.path = "../data/processed/", center=TRUE, scale=TRUE, ret=T){
  
  #############################################
  # Written by MJB on Sept 20th, 2019, last modified on 
  # dat = data for performing the PCA
  # gp = factor of the grouping
  # pc.x = number, PC to display on x-axis, default = 1 
  # pc.y = number, PC to display on y-axis, default = 2 
  # col = vector of colors to be used
  # expt.nm = experiment name to be used in the subtitle and export file 
  # center = logical, center the data
  # scale = logical, scale the data
  # ret = logical to return the plot
  #############################################
  
  
  library(ggplot2)
  
  #Make PCA of original data
  pca <- prcomp(t(dat), center=center, scale=scale)
  
  #Get varaince percentages
  pc.pca <- pca$sdev^2 / sum(pca$sdev^2) * 100; 
  pc.pca <- round(pc.pca, digits=2)
  
  #Plot PCA
  p <- ggplot(data.frame(pca$x), aes_string(paste0("PC", pc.x), paste0("PC", pc.y))) +
      geom_point(aes(color = gp), size = 6) +
      labs(size = NULL, color = "Group", 
           x=paste0("PC", pc.x, " (", pc.pca[pc.x],"%)"), 
           y=paste0("PC", pc.y, " (", pc.pca[pc.y],"%)")) +
      geom_text(aes(label = gsub(".+_(\\d)$", "\\1", row.names(pca$x)))) +
      scale_color_manual(values = unique(clr)) +
      labs(subtitle = paste0("PCA - ", expt.nm))
  
  # Export plot
  pdf(paste0(res.path, expt.nm, "_pca.pdf"), 
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
# dat <- dge$lcpm
# gp <- gene_dge_e2$samples$group
# clr <- colorsSamp_e2[4:12]
# expt.nm <- "Expt_2"
# dge <- gene_dge_e2
# center=T
# scale=T
# res.path = "../data/processed/"
# ret = T
####################################
