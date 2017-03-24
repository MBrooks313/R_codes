GOANA_Summary <- function(output = "Output", name.pre=enrich.names, name.suf="_GOANA.csv", ontology = "BP", FDR=0.01) {
  ##This program was written to reduce the inherent redundancy associated with Gene Ontology and Kegg enrichment analysis.
  ##It examines the parent-child relationship of significantly enriched terms and reports back the leaf (end-of-branch) term
  ##as the representative enriched ontology for each branch enriched.
  ##This is the 2nd of 3 programs for GO/KEGG enrichment.
  ##Written by Matthew Brooks, last modified March 23rd, 2017
  
  ####USAGE####
  # output = prefix of the output file name
  # name.pre = names from the list (or list-of-lists) from gene.list of the GOANA.R function, prefix for import file
  # name.suf = suffix for import file
  # ontology = 'BP', 'CC', 'MF', or 'KG' for ontology to summarize 
  # FDR = FDR cutoff to subset prior to examining the parent-child relationship
  #
  # Examples:
  # Output from GOANA.R to be used as input: list1_GOANA.csv
  #
  # GOANA_Summary(output="GOANA-BP", name.pre = names(cl.list), name.suf = "_GOANA.csv", ontology = "BP", FDR=0.01)
  # GOANA_Summary(output="GOANA-CC", name.pre = names(cl.list), name.suf = "_GOANA.csv", ontology = "CC", FDR=0.01)
  # GOANA_Summary(output="GOANA-MF", name.pre = names(cl.list), name.suf = "_GOANA.csv", ontology = "MF", FDR=0.01)
  # GOANA_Summary(output="KG", name.pre = names(cl.list), name.suf = "_KGANA.csv", ontology = "KG", FDR=0.01)
  #############
  
  require(GO.db)
  
  #Load ontology
  if (ontology == "BP") GOCHILDREN <- GOBPCHILDREN
  if (ontology == "CC") GOCHILDREN <- GOCCCHILDREN
  if (ontology == "MF") GOCHILDREN <- GOMFCHILDREN
  
  
  global.list <- lapply(1:length(name.pre), function(x){
  
  	# Get input, filter for ontology, and calculate FDR
    input <- read.csv(paste(name.pre[x], name.suf, sep="")) #Import of the output
    
    if (ontology == "KG") {
    	input <- input[input$DE > 0,]
    } else {
    	input <- input[input$Ont == ontology & input$DE > 0,]
    }
    input <- input[input$FDR <= FDR,]
	
	  kg <- data.frame()
	  if (ontology == "KG") {
		
		  kg <- rbind(kg, input) 
		  return(kg)
		
	  } else {
	
	    ## Get Top End-of-Branch GO Term in Each Cluster
	    
	    # Loop for each cluster
	    top <- data.frame()
	      
	    tmp.GO <- as.character(input[,1]) #List of GO terms
	      
	    #Loop for each cluster with a significant GO term 
		  if (length(tmp.GO) > 0){
		    tmp.eob <- c()
      
		    #Parse for single GO term
		    if (length(tmp.GO) == 1){
		      tmp.cluster <- input[input$X == tmp.GO,] #If only a single sig GO term
		    
		      #Parse for more than one sig GO term
		    } else{
		      for (k in 1:length(tmp.GO)){
        
		        if (!is.null(GOTERM[[tmp.GO[k]]])) { tmp.offspring <- unique(unlist(mget(tmp.GO[k], GOCHILDREN)))
		      
		        #Keep only End-of-Branch GO Terms
		        if (sum(tmp.offspring %in% tmp.GO) == 0) {tmp.eob <- c(tmp.eob, tmp.GO[k])}
		        }
		    }
		    
		    tmp.cluster <- input[input$X %in% tmp.eob,] #Get results back
		    
		  }
		  
		  top <- rbind(top, tmp.cluster) 
		  
		}
	    
	    return(top)
	}
  })  
  
  names(global.list) <- name.pre
  
  #Write output file
  if (ontology == "KG"){
  	write.csv(do.call("rbind", global.list), file=paste(output ,"_KG_Summary.csv", sep=""))
  
  } else {
  	
  	write.csv(do.call("rbind", global.list), file=paste(output ,"_GO_Summary.csv", sep=""))
  }
}
