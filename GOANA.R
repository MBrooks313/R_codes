GOANA <- function(gene.list, univ, species){
  ##This function performs Gene Ontology and KEGG enrichment analysis within the limma package. The user supplies a named 
  ##list (or list-of-lists) of ENTREZ gene IDs for clusters (or gene subsets) and a vector of ENTREZ gene IDs to be used as 
  ##a background list (expressed genes). Benjamini-Hochberg false discovery rate (FDR) is calculated. It exports a csv file of the complete enrichment analysis.
  ##This is the 1st of 3 programs for GO/KEGG enrichment.
  ##***WARNING*** This function takes several minutes to run depending on how many lists are being analyzed.
  ##This written by Matthew Brooks, last modified March 23rd, 2017
  
  ####USAGE####
  # gene.list = list (or list-of-list) of ENTREZ gene IDS.
  # univ = list of ENTREZ gene IDs for 
  # species = 'hsa', 'mmu', or 'rno'
  #############
  
  require(limma)
  
  #Load the species database
  if (species == "Hs") {require(org.Hs.eg.db); x<-as.list(org.Hs.egGO2ALLEGS); xx<-as.list(org.Hs.egSYMBOL2EG); k<-"hsa"}
  if (species == "Mm") {require(org.Mm.eg.db); x<-as.list(org.Mm.egGO2ALLEGS); xx<-as.list(org.Mm.egSYMBOL2EG); k<-"mmu"}
  if (species == "Rn") {require(org.Rn.eg.db); x<-as.list(org.Rn.egGO2ALLEGS); xx<-as.list(org.Rn.egSYMBOL2EG); k<-"rno"}
  #if (species == "Dm") require()
  #if (species == "Pt") require()
  
  name.list <- names(gene.list) #for list of IDs
  
  #Loop through each list given
  for (i in 1:length(name.list)){
    list <- gene.list[[i]] #for list of IDs
    
    #Perform Gene Ontology enrichment
    go <- goana(list, universe=univ, species=species)
    go <- go[order(go$P.DE),]
    rank.go <- 1:dim(go)[1]
    
    #Perform Kegg pathway enrichment
    kegg.g <- getGeneKEGGLinks(species.KEGG = k, convert = TRUE)
    kegg.p <- getKEGGPathwayNames(species.KEGG = k, remove.qualifier = FALSE)
    kg <- kegga(list, species=species, gene.pathway = kegg.g, pathway.names = kegg.p)
    kg <- kg[order(kg$P.DE),]
    rank.kg <- 1:dim(kg)[1]
    
    #Calculate FDR
    go$FDR <- go$P.DE*dim(go[go$DE!=0,])[1]/rank.go
    kg$FDR <- kg$P.DE*dim(kg[kg$DE!=0,])[1]/rank.kg
    
    #Get gene names of GO enrichment
    for (j in 1:dim(go)[1]){
      z <- gene.list[[i]][gene.list[[i]] %in% x[row.names(go[j,])][[1]]]
      go$Genes[j] <- paste(sort(names(xx[match(z,xx)])), collapse=",")
    }
    
    #Get gene names of Kegg enrichment
    for (j in 1:dim(kg)[1]){
      z <- gsub(".+(\\d....)$", "\\1",row.names(kg[j,]), perl = T)
      zz <- gene.list[[i]][gene.list[[i]] %in% kegg.g[grep(z, kegg.g$PathwayID),1]]
      kg$Genes[j] <- paste(sort(names(xx[match(zz,xx)])), collapse=",")
    }
    
    #Write files
    write.csv(go, file=paste(paste(name.list[i], "GOANA", sep="_"),".csv", sep=""))
    write.csv(kg, file=paste(paste(name.list[i], "KGANA", sep="_"),".csv", sep=""))
  }
  
  
}
