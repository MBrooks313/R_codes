###################################
# Multi-dot plot function

dot_go <- function(dat, top, grp, desc, p.val, ratio, filenm){ 
  
  ##############################################
  #Written by MJB on Sept 14th, 2018, last modified on OCT 25th, 2018
  #dat - data.frame
  #top - top number of categories to display in barplot 
  #grp - columnn name containing the x-axis variable, character
  #desc - column name containing the pathway description, character
  #p.val - column name with p values, character
  #ratio - column name with the gene-pathway ratio, character
  #filenm - name for pdf file export of plot
  ##############################################
  
  library(ggplot2)
  library(gtools)
  library(tidyr)
  library(dplyr)
  
  
  #Select on the the top specified number
  if (!is.na(top)) {
    dat <- dat %>% group_by_(grp) %>% arrange_(p.val) %>% top_n(top, dplyr::desc(!!as.symbol(p.val))) %>% as.data.frame
  }
  
  #Prep the data table
  # dat[,`grp`] <- factor(dat[,`grp`], levels = mixedsort(levels(factor(dat[,`grp`])))) #Factorize the group using a natural sort
  dat <- dat[order(dat[,`grp`], dat[,`p.val`]),] #Order by group and p value
  dat[,`desc`] <- substring(dat[,`desc`], 1, 60) #Cut description name to 60 characters
  dat[,`desc`] <- factor(dat[,`desc`], levels = rev(unique(dat[,`desc`]))) #Factorize the order of descriptions
  dat[,`p.val`] <- -log10(dat[,`p.val`])
  
  #Open the pdf device
  pdf(paste(filenm, "_dotplot.pdf", sep=""), useDingbats = F,
      height = dim(dat)[1]/10 + 2.5,
      width = 8)
  
  #Make the plot
  p <- ggplot(dat, aes_string(grp, desc, colour = p.val)) + 
      geom_point(aes_string(size = ratio)) + 
      theme_bw() +
      scale_colour_gradient(low = "#3B9AB2", high = "#F21A00", name = "P Value\n-log10", limits=c(0,10), oob = scales::squish) +
      guides(size=guide_legend(title="Ratio")) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1))
  
  print(p)
  
  dev.off()
}
###################################


#####################
#Testing
# dat = go.bp
# top = 10
# grp = "Cluster"
# desc = "Description"
# p.val = "p.adjust"
# ratio = "ratio"
#####################
