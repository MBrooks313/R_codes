ParPlot<-function (eset, cl, ord=NA, mfrow = c(1, 1), col.cl = "blue", xlabels, ylabel, file.nm=fname) {
    ##This function creates a pdf file of parallel plots with a mean-line for your cluster analysis.
    ##Written by Matthew Brooks, last modified March 23rd, 2017.
    
    ####USAGE####
    #eset =	matrix of expression values, e.g. rows (gene) x columns (samples)
    #cl =	vector of clusters for each row in the matrix
    #ord =	order for the cluster plots to be displayed, numerical by default
    #mfrow =	rows and columns for the final plot arrangement
    #col.cl = vector of mean-line colors
    #xlabels =	vector of x-axis tick lables for the plots, e.g. c('P2', 'P4', 'P6', 'P10')
    #ylabel = character for the y-axis title label, e.g. "Z-score", "CPM (log2)" 
    #file.nm =	name for pdf
    #############
    
    clusterindex <- cl
    
    #Specify the ordering of the plots
    if(!is.na(max(ord))){
    }else{
    	ord=c(1:max(cl))
    }    		
    
    #PDF of the parallel plots
    pdf(paste(file.nm, "pdf", sep = "."), width = 3 * mfrow[2], height = 3 * mfrow[1], useDingbats = F) 
    par(mfrow = mfrow)
    for (j in ord) {
    	tmp <- eset[clusterindex == j, ] 
    	ymin <- min(tmp)
    	ymax <- max(tmp)
    	par(las=2)
      
      #Main plot command
    	plot.default(x = NA, xlim = c(1, dim(eset)[[2]]), 
                ylim = c(ymin, ymax), xlab = "Age", ylab = "Z-score", cex.lab=1.5,
                main = paste("Cluster ", j, " (n = ",dim(tmp)[[1]],")", sep=""), cex.main=2, axes = FALSE)
    	
      #Add labels to x-axis
      if (missing(xlabels)) {
                axis(1, 1:dim(eset)[[2]], c(1:dim(eset)[[2]]))
                axis(2)
    	} else {
            	 axis(1, 1:dim(eset)[[2]], xlabels)
               axis(2)
        }
      
      #Add data lines to the plot    
    	for (k in 1:dim(tmp)[[1]]) {
    		lines(c(1:dim(tmp)[[2]]), tmp[k,], col = "gray")
    	}
      
      #Add mean-line
    	lines(apply(tmp,2,mean), col = col.cl[j], lwd=2)
    }
	dev.off()

}
