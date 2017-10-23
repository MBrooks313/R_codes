ExpFilt <- function(data, group, exp=1){
  ##This program expression filtering on a wide format data.frame. It requires all members of a group to be greater than or equal to the expression value given.
  ##Written by Matthew Brooks, last modified March 23rd, 2017
	
  ####USAGE###
  # data: data.frame of expression values
  # group: list of groupings for samples in data
  # exp: expression level to require all samples in a group to be equal to or greater
  #
  # Return an index of rows passing the filtering criteria
  ############
	
  #Define variables
  mtx <- data
  fct <- as.factor(group)
  lvl <- levels(fct)
  	
  #Loop for each level of group
  idx <- list()
  for (i in 1:length(lvl)) {
  		
        #If group has only one sample
   	if (length(grep(lvl[i], group)) == 1) {
     		idx.tmp <- row.names(mtx[mtx[,grep(lvl[i], group)] >= exp,])
   	}
    	
   	#If group has more than one sample
   	else {
     		idx.tmp <- row.names(mtx[rowSums(mtx[,grep(lvl[i], group)] >= exp) == length(grep(lvl[i], group)),])
   	} 
   	idx <- c(idx,idx.tmp)
    }
  	
   #Get unique list of indexes
   idx <- unique(unlist(idx))
  	
   #Return varaible containing the list of indexes passing the filtering
   assign("idx.filt", idx ,envir=.GlobalEnv)
  }
