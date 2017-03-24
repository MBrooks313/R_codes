DFave <- function(dat, grp, ord=NULL) {
  ##Created by Matthew Brooks, last modified March 23rd, 2017
  ##This function uses the stats::aggregate function to average a wide format data.frame based on given grouping. 
  ##In addition, it can re-order the returned, averaged data.frame.  
  
  ####USAGE####
  #dat = data.frame of numeric data
  #grp = factor of grouping to average
  #ord = vector of order from levels(group) for final ordering
  #############
  
  #Aggregate function
  agg <- aggregate(t(dat), by=list(grp), FUN="mean")
  
  #Re-order averaged data.frame
  if (is.null(ord)) {
    agg <- t(agg[,-1])
    colnames(agg)<- levels(grp)
  } else {
    agg <- t(agg[ord,-1])
    colnames(agg)<- levels(grp)[ord]
  }
  
  agg <- data.frame(agg)
}
