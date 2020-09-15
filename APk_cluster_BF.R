##############################
# Performs AP cluster analysis


cluster.apK <- function(zscore, nk = 2, c2u = 1:ncol(zscore), pc=1) {
  
    #############################################
    # Written by Ben Fadl, modified by Matthew J. Brooks
    # v0.1.0 tested on June 10th, 2019 by MJB
    #
    # zscore : data.frame of the z-score data to be clustered
    # nk : number of clusters to be found
    # c2u : columns to be used in clustering
    #############################################
    
  
    #Packages required
    require("apcluster")
    
    #PrincipleComponent sort function
    pc_sort <- function(df,
                        cols = 1:ncol(df),
                        pc = pc) {
      
      
      # #Testing variables
      # df <- zscore[k[[1]],]
      # cols = 1:ncol(df)
      # ##################
      
      #############################################
      # Written by Ben Fadl, modified by Matthew J. Brooks
      # v0.1.0 tested on June 10th, 2019 by MJB
      #
      # df : data.frame of the data to be sorted by principal component
      # cols : columns to be used in calculating the order
      # pc : pricipal component to be used in sorting
      #############################################
      
      #Perform PCA analysis on df
      pco <- prcomp(df[, cols])
      
      #sort the df by PC
      temp.df <- as.data.frame(pco$x)
      index <- sort(temp.df[, pc], index.return = T)$ix
      
      #Get the row.names of the PC ordered df
      ordered.rn <- rownames(temp.df[index, ])
    }
    
    #Make sure the data is a data.frame
    zscore <- data.frame(zscore)
    
    #Perform AP cluster with the number of clusters requested (k)
    k <- apclusterK(corSimMat(), zscore[, c2u], nk)@clusters
    
    #Assign rowname to clusters
    rn <- c()
    for (c in 1:length(k)) {
      rn <- c(rn, k[[c]])
    }
    
    #Sort genes by by cluster number and by PC1 within each cluster
    rn.sort <- c()
    for (i in 1:length(k)) {
      rn.sort <- c(rn.sort, pc_sort(zscore[k[[i]],], pc=pc))
    }
    zscore <- zscore[rn.sort,]
    
    #Make data.frame for export of gene's cluster assignment
    zscore.c <- data.frame(matrix(nrow = 0, ncol = 2))
    colnames(zscore.c) <- c("rn", "cluster")
    for (i in 1:length(k)) {
      zscore.c <- rbind(zscore.c, data.frame(rn = names(k[[i]]), cluster = i))
    }
    
    #Assign and return values
    assign("clust.df", zscore.c, envir = .GlobalEnv)
    assign("clust.list", lapply(k, names), envir = .GlobalEnv)
    assign("zscore.sort", zscore, envir = .GlobalEnv)
}


# #Testing variables
# zscore <- de.anova$aveSigZ[1:500,]
# nk = 8
# c2u = 1:ncol(zscore)
# pc = 2
# ##################
