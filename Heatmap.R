
heatmap_MB <- function(expdat, sym_col, dat_col, sort_col, filenm, ret = FALSE, base="CPM", cluster = "PCA", symb = TRUE, color=NA) {
  
  #############################################
  # Written by MJB on Oct 24th, 2018, last modified on Sept 19th, 2019
  # expdat = expression data.frame with annotation
  # sym_col = column number for the gene symbol from expdat
  # dat_col = column number for the data to be plotted
  # sort_col = column numbers to be used in sorting
  # filenm = prefix for file name generation
  # ret = logical to return the ggplot for further modification in the envirn
  # base = whether plot is in log2 CPM or Z-scores, c("CPM", "Zscore"), default = "CPM"
  # cluster = gene ordering method, c("PCA", "HC"), default = "PCA"
  # symb = logical to show gene symbols on plot
  # color = color scheme to override default for use in the heatmap
  #############################################
  

  library(ggplot2)
  library(tidyr)
  library(gtools)
  library(dplyr)
  library(scales)
  
  ##############################################################################
  # Remove any rows ith duplicated gene names
  if (sum(duplicated(expdat[,5])) > 0){
    print("Duplicated gene names...generating plot...check original df.")
    expdat <- expdat[!duplicated(expdat[,5]), ]
  }
  ##############################################################################
  
  
  ##############################################################################
  # Assign scale, legend header, and color for Z-score or CPM
  if (base == "Zscore"){
    
    lim = c(-3,3)
    nm = "Z-score"
    jet <- c("#7F0000","red","#FF7F00","yellow","white","cyan", "#007FFF", "blue","#00007F")
    mb.colors <- colorRampPalette(rev(jet))(1000)
    
  } else {
    
    lim = c(0,16)
    nm = "CPM (log2)"
    mb.colors <- colorRampPalette(c("#023858","#3690c0", "#e0f3f8","#fff7bc","#fec44f","#ec7014",
                                    "#fc4e2a","#bd0026","#800026"))(n=1000)
    
  }
  
  #Override the default color schemes with user supplied
  if (!is.na(color)){
    mb.colors <- colorRampPalette(color)(n=1000)
  }
  
  ##############################################################################
  # Order the plot
  if (cluster == "HC"){
    
    #Heirarchical clustering using Ward's method with Euclidean distanced
    distfunc <- function(x) dist(x,method="euclidean")
    hclustfunc <- function(x) hclust(x, method="ward.D2")
    
    #Get distance and heirarchical clustering
    d <- distfunc(as.matrix(expdat[,sort_col]))
    hc <- hclustfunc(d)
    
    #Tidy table
    dat.tidy <- expdat[,c(sym_col, dat_col)] %>%
      gather(Sample, CPM, -external_gene_name)
    
    #Factor sort sample names and gene name
    dat.tidy$Sample <- factor(dat.tidy$Sample, levels = unique(dat.tidy$Sample))
    dat.tidy$external_gene_name <- factor(dat.tidy$external_gene_name, 
                                          levels = dat.tidy$external_gene_name[hc$order], ordered = T)
    
  } else {
    
    # PCA method of ordering
    pca <- prcomp(expdat[,sort_col], center = T, scale. = F)
    pca.order <- order(pca$x[,1])
    
    #Tidy table
    dat.tidy <- expdat[,c(sym_col, dat_col)] %>%
      gather(Sample, CPM, -external_gene_name)
    
    #Factor sort sample names and gene name
    dat.tidy$Sample <- factor(dat.tidy$Sample, levels = unique(dat.tidy$Sample))
    
    #Check to make ensure plot is decreasing from top to bottom
    dat.tmp <- dat.tidy
    dat.tmp$external_gene_name <- factor(dat.tmp$external_gene_name, levels = dat.tmp$external_gene_name[pca.order], ordered = T)
    
    first.mean <- dat.tmp %>%
      filter(external_gene_name == levels(external_gene_name)[1]) %>% dplyr::summarize(mean(CPM))
    
    last.mean <- dat.tmp %>%
      filter(external_gene_name == levels(external_gene_name)[length(pca.order)]) %>% dplyr::summarize(mean(CPM))
    
    if (first.mean > last.mean) {
      dat.tidy$external_gene_name <- factor(dat.tidy$external_gene_name, levels = dat.tidy$external_gene_name[rev(pca.order)], ordered = T)
    } else {
      dat.tidy$external_gene_name <- factor(dat.tidy$external_gene_name, levels = dat.tidy$external_gene_name[pca.order], ordered = T)
    }
    
  }
  ##############################################################################
  
  
  ##############################################################################
  # Generate plot
  p <- ggplot(dat.tidy, aes(Sample, external_gene_name)) +
        geom_tile(aes(fill = CPM), colour = "transparent") + 
        #scale_fill_gradient(high = "white", low = "steelblue", name = "CPM (log2)") +
        scale_fill_gradientn(colours = mb.colors, name = nm, limits = lim, oob=squish) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_blank(), axis.title.y = element_blank(),
              panel.background = element_rect(fill = NA))
        
  # Remove gene names and ticks if necessary
  if (symb == FALSE){
    p <- p +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  # Export plot
  pdf(paste(filenm, "_hm.pdf", sep=""), 
      width = dim(unique(dat.tidy[2]))[1] / 2 + 2, 
      height = dim(unique(dat.tidy[1]))[1] / 6 + 1, 
      useDingbats = F)
  
  print(p)
  
  dev.off()
  
  # Return plot if needed
  if (ret == TRUE){
    return(p)
  }
  
  
}


#############################################
# Change log for 1.7
# Added HC sorting
# Added including gene symbol or not in plot
# Depricated "ord" function
# 
# ***Testing variables at bottom 
#############################################

#Testing variables #################
# load("/Volumes/PEGASUS/Projects/BioInf_Scripts/R_functions/ggplots/Test_CPM-Z_datasets.Rdata")
# #expdat = gene.mstr.de.ave.z
# expdat = gene.mstr.de.ave.z
# sym_col = 5
# dat_col = c(10:12)
# sort_col = c(10:12)
# ret = F
# base = "Zscore"
# cluster = "HC"
# #ord = T
# symb = FALSE
# color = NA
# filenm = "R_plots/TF_comb-sig_HMv2"
#expdat = gene.mstr.de.ave.z; sym_col = 5; dat_col = c(10:12); sort_col = c(10:12); base="Zscore"; cluster = "HC"
####################################