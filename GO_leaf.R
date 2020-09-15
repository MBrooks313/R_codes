###########################################
# Reduce GO redundancy by taking leaf term

GO_leaf <- function(godf, query = NA, id.col, ont = "BP") {
  
  #############################################
  # Written by MJB on Feb 26th, 2019, last modified on Feb 26th, 2019
  # This function reduces Gene Ontology redundancy by taking a result data.frame containing an 
  # output of significant GO terms and finds/returns a subsetted data.frame of only the "leaf terms".
  # 
  # godf = data.frame of GO terms for redundancy reduction
  # query = column number with query or cluster identifier. NA if singular enrichment was performed.
  # id.col = column number with GO IDs, IDs will have the format: GO:1234567
  # ont = which ontology the reduction is to be performed, one of c("BP", "MF", "CC")
  #############################################
  
  
  library(GO.db)
  
  ############ GO Leaf Function Definition ############
  leaf <- function(tmpdf = tmpdf){
    
    #Define GOCHILDREN
    if (ont == "BP") GOCHILDREN <- GO.db::GOBPCHILDREN
    if (ont == "CC") GOCHILDREN <- GO.db::GOCCCHILDREN
    if (ont == "MF") GOCHILDREN <- GO.db::GOMFCHILDREN
    
    tmp.GO <- as.character(tmpdf[,id.col]) #List of GO terms
    tmp.cluster <- c()
    
    #Loop for each cluster with a significant GO term 
    if (length(tmp.GO) > 0){
      
      #Set empty end-or-branch list
      tmp.eob <- c()
      
      #If only a single sig GO term
      if (length(tmp.GO) == 1){
        
        tmp.cluster <- godf 
        
        #Parse for more than one sig GO term
      } else{
        
        #Loop for each GO term
        for (k in 1:length(tmp.GO)){
          
          #Get all children of a GO term
          if (!is.null(GO.db::GOTERM[[tmp.GO[k]]])) { 
            
            tmp.offspring <- unique(unlist(AnnotationDbi::mget(tmp.GO[k], GOCHILDREN)))
            
            #Keep only End-of-Branch GO Terms
            if (sum(tmp.offspring %in% tmp.GO) == 0) {tmp.eob <- c(tmp.eob, tmp.GO[k])}
            
          }
          
        }
        
        tmp.cluster <- tmpdf[tmpdf[,id.col] %in% tmp.eob,] #Get results back
        
      }
      
    }
    
    return(tmp.cluster)
  }
  ###################################################
  
  
  #Prep output table
  output <- data.frame()
  
  #If query is NA 
  if (is.na(query)){
    
    output <- leaf(godf)
  
    #If query is for multiple clusters, etc. 
  } else {
    
    for (i in unique(godf[,query])) {
      
      tmpdf <- godf[godf[,query] == i,]
      indi.df <- leaf(tmpdf)
      output <- rbind(output, indi.df)
      
    }
    
  }
  
  return(output)
  
}



#Testing variables #################
# godf <- read.csv("~/Dropbox/John_RPEscaffold/R_output/gProfileR_PC2.tsv", sep = "\t", row.names = 1) #Single query
# godf <- read.csv("/Volumes/PEGASUS/Projects/Manju_MPSI/Compounds_Jan2019/R_output/gProfileR_Cluster.tsv", sep = "\t", row.names = 1) #Multiple query
# godf <- godf[godf$domain == "BP",]
# query = 1
# id.col = 9
# ont = "BP"
####################################

