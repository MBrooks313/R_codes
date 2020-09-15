###################################
#Single BarPlot function

bar_go <- function(dat, top, desc, p.val, genenum, filenm) {
  
  ##############################################
  #Written by MJB on Sept 14th, 2018, last modified on Oct 25th, 2018
  #dat - data.frame
  #top - number of categories to display in barplot 
  #desc - column name containing the pathway description, character
  #p.val - column name with p values, character
  #genenum - column name with the gene number in category, character
  #filenm - name for pdf file export of plot
  ##############################################
  
  
  library(ggplot2)
  
  id_col = "intersection"
  
  #Prep data table
  dat <- dat[order(dat[,`p.val`]),] #Order by group and p value
  dat[,`desc`] <- substring(dat[,`desc`], 1, 60) #Cut description name to 60 characters
  dat[,`desc`] <- factor(dat[,`desc`], levels = rev(unique(dat[,`desc`]))) #Factorize the order of descriptions
  dat$p.value <- -log10(dat$p.value)
  
  #Select only the category number wanted
  if (top > dim(dat)[1]){
    dat <- dat
  } else {
    dat <- dat[1:`top`,]
  }
  
  #Get intersection length
  ids <- lapply(dat[,`id_col`], function(x) {unlist(strsplit(x, split = ","))})
  ids.nm <- unlist(lapply(ids, length))
  
  #Open the pdf device
  pdf(paste(filenm, "_barplot.pdf", sep=""), useDingbats = F,
      #height = dim(dat[1:`catnum`,])[1] / 5 + 0.75,
      height = dim(dat)[1] / 5 + 1,
      width = 7)
  
  #Make the plot
  q <- ggplot(dat) +
    geom_bar(stat = 'identity', aes_string(x=desc, y = ids.nm, fill = p.val)) +
    coord_flip() +
    scale_fill_gradient(low="#3B9AB2", high="#F21A00", name="P Value \n(-log10)",limits=c(0,10), oob = scales::squish) +
    #scale_fill_gradientn(colors = c("#3B9AB2", "#F21A00"), name="P Value (log10)", trans = 'reverse', limits=c(-50,0), oob=squish) +
    labs(y="Number of Significant Genes", x=NULL) +
    theme_bw()
  
  print(q)
  
  dev.off()
}
###################################


###############
#Testing
# p.val = "p.value"
# top = 10 
# desc = "term.name" 
# p.val = "p.value" 
# genenum = "term.size"
###############