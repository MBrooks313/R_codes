parPlot_embed<-function (eset, cl, ord=NA, mfrow = c(1, 1), col.cl, time.labels=NULL,fn=fname) {
    #eset=	matrix of expression values
    #cl=	cluster list for each transcript in the matrix
    #ord=	order for the cluster plots to be displayed, numerical by default
    #mfrow=	rows and columns for the plots
    #time.labels=	x-axis lables for the plots
    #fn=	name for pdf
    
    library(pryr)
  
    clusterindex <- cl
    
    if(!is.na(max(ord))){
    }else{
    	ord=c(1:max(cl))
    }    		
    #center <- cl[[2]]
    
    #pdf(paste(fn,"pdf",sep="."),width=3*mfrow[2],height=3*mfrow[1]) 
    parplot %<a-% {
      par(mfrow = mfrow, mar=c(5,4,2,1))
      for (j in ord) {
        tmp <- eset[clusterindex == j, ] 
        ymin <- min(tmp)
        ymax <- max(tmp)
        #par(las=2, mar = c(8,5,2,3))
        plot.default(x = NA, xlim = c(1, dim(eset)[[2]]), 
                     ylim = c(ymin, ymax), xlab = "", ylab = "", cex.lab=1.5,
                     main = paste("Cluster ", j, " (n = ",dim(tmp)[[1]],")", sep=""), cex.main=1, axes = FALSE)
        title(ylab = "Z-score", line=3, cex.lab = 1.5)
        if (is.null(time.labels)) {
          axis(1, 1:dim(eset)[[2]], c(1:dim(eset)[[2]]), labels = FALSE)
          text(1:dim(eset)[[2]], par("usr")[3] - 0.5, labels = c(1:dim(eset)[[2]]), srt = 45, pos = 1, xpd = TRUE)
          axis(2)
        }
        else {
          axis(1, 1:dim(eset)[[2]], labels = FALSE)# time.labels, srt=45, las=2)
          text(1:dim(eset)[[2]], par("usr")[3] - 0.5, labels = time.labels, srt = 45, pos = 1, xpd = TRUE)
          axis(2)
        }
        for (k in 1:dim(tmp)[[1]]) {
          lines(c(1:dim(tmp)[[2]]), tmp[k,], col = "gray")
        }
        lines(apply(tmp,2,mean), col = col.cl[j], lwd=3)
      }
    }
    
	#dev.off()
    #assign("parplot", parplot, envir = .GlobalEnv)
    return(parplot)
}