GOANA_Plot <- function(file.in="GOANA_GO_Summary.csv", ana = "GO", rlvl.x = NA, rlvl.y=NA, output = "GOANA_plot", x=12, y=12, top=NA){
  ##This program was written to generate a plot output from GOANA_Summary.R
  ##This is the 3rd of 3 programs for GO/KEGG enrichment.
  ##GOANA_Plot.R written by Matthew J. Brooks, last modified March 23rd, 2017
  
  ####USAGE####
  # file.in = file from GOANA_Summary output, ex. "GOANA-BP_GO_Summary.csv"
  # ana = either "GO" for gene ontology or "KG" for Kegg
  # rlvl.x = numberic list contaning the rearanged row numbers for the x-axis, ex. c(3,5,6,1,2,4)
  # rlvl.y = numberic list contaning the rearanged row numbers for the y-axis, ex. c(7,6,8,1,10,5,11,9,3,12,4,2)
  # output = character string of the file name for the plot, example: "GOANA_plot"
  # x = width of output pdf in inches, example: 12
  # y = height of output pdf in inches, example: 12
  # top = the maximum number of GO categories per group to display

  goana.sum <- read.csv(file.in)
  
  # Plot function
  gk.plot<-function(x){ 
    require(ggplot2)
    p<-ggplot(x, aes(x=TP, y=factor(GO.term, levels=rev(GO.term)), color = FDR, size=Ratio), ) 
    p <- p + geom_point() + scale_colour_gradient(low = "red", high = "blue")
    style <- element_text(color = "black", size = 10)
    p<-p+theme(axis.text.x = style, axis.text.y = style, 
               axis.title.x = element_blank(), axis.title.y = element_blank())
    plot (p)
  }

  #Make data.frame for plot function
  if (ana == "KG"){
  	GO.tmp <- data.frame(goana.sum[,2:3])
  	TP.tmp <- data.frame(do.call('rbind', strsplit(as.character(goana.sum[,1]), '\\.\\d+', perl=TRUE,fixed=FALSE)))
  	goana.plot <- data.frame(GOIDs=GO.tmp$X, GO.term=GO.tmp$Pathway, FDR=goana.sum$FDR, 
                                   TP=TP.tmp[,1], Ratio=goana.sum[,5]/goana.sum[,4])
  } else {
  	GO.tmp <- data.frame(goana.sum[,2:3])
  	TP.tmp <- data.frame(do.call('rbind', strsplit(as.character(goana.sum[,1]), '\\.\\d', perl=TRUE, fixed=FALSE)))
  	goana.plot <- data.frame(GOIDs=GO.tmp$X, GO.term=GO.tmp$Term, FDR=goana.sum$FDR, 
                                   TP=TP.tmp[,1], Ratio=goana.sum[,6]/goana.sum[,5])
  }
  
  #Take the top members of each cluster, if necessary
  if (!is.na(top)[1]) {
      top.tmp <- c()
      top.rows <- as.factor(goana.plot$TP)
      for (i in levels(top.rows)){
          top.tmp <- rbind(top.tmp, goana.plot[goana.plot$TP == i,][1:top,])
          top.tmp <- top.tmp[!is.na(top.tmp$GOIDs),]
      }
      goana.plot <- top.tmp
  }
  
  #Relevel the x-axis
  if (!is.na(rlvl.x)[1]){
  	ord <- rlvl.x[rlvl.x %in% goana.plot$TP]
    ord <- match(ord, levels(goana.plot$TP))
    goana.plot$TP <- factor(goana.plot$TP, levels=levels(goana.plot$TP)[ord])
    goana.plot <- goana.plot[order(goana.plot$TP),]
  }
  
  #Relevel the y-axis
  if (!is.na(rlvl.y)[1]){
    goana.plot$GO.term <- factor(goana.plot$GO.term, levels=levels(goana.plot$GO.term)[rlvl.y])
  }
  
  #Plot graph
  pdf(file=paste(output, ".pdf", sep=""), width=x, height=y, useDingbats = F)
  gk.plot(goana.plot)
  dev.off()

}
