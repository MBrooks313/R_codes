TopTranscript <- function(data, data.cols, list, colnm, row.ord=F) {
	##This program searches a wide format data.frame containing annotation and expression values for queried gene names and reports back the highest expressing item of any rows having duplicate gene symbols. This is useful for reporting back the highest expressing transcript as the gene representative for further anaysis.
  ##This was written by Matthew Brooks, last modified March 23rd, 2017.
  
	####USAGE####
	# data = an annotated data.frame of annotation and expression data
	# data.cols = columns of expression data to evaluate the highest transcript, ex. c(11:25)
	# list = character vector of gene symbols to be searched for in the data, alternatively this could be a vector of row.names
	# colnm = column number to search for duplicates, i.e. column containing gene symbols
	# row.ord = boolean, if the final output should be ordered alphabetically by the gene symbol
	#########################
	
	goit.list <- data.frame()
	
	for (i in 1:length(list)) {
		goi <- list[i]
		
		# Pulls out all matching rows
  	goit <- data[grepl(paste("^", goi, "$", sep=""), data[,colnm], ignore.case=TRUE),]
  		
  	
  	# Keeps the highest value transcript
  	if (dim(goit)[1] != 0){
  		x <- which(goit[data.cols] == max(goit[data.cols]), arr.ind=TRUE)
      
  		if (dim(x)[1] > 1){
  			goit <- goit[1,]
  		} else {
    		goit <- subset(goit, row.names(goit) == row.names(which(goit[data.cols] == max(goit[data.cols]), arr.ind=TRUE)))
  		}
    	goit.list <- rbind(goit.list, goit)
    		
    	# Re-orders the rows alphabetically
    	if (row.ord == T){
    		goit.list <- goit.list[order(goit.list[,colnm]),]
    	}
  	}
	}
	return(goit.list)
}
