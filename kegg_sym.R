
kegg_sym <- funtion(keggID){
  
  ##############################################
  # Written by Matthew J Brooks on August 26th, 2024
  # Provide KEGG pathway ID and return gene symbols belonging to pathway
  #
  # keggID : Kegg pathway ID (character, ie. "mmu00010")
  ##############################################
  
  library(KEGGREST)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  
  orgdb <- org.Mm.eg.db::org.Mm.eg.db
  
  kegg.path <- keggList("pathway")
  kegg.names <- gsub("map", "mmu", names(kegg.path))
  kegg.path <- unlist(lapply(1:length(kegg.path), function(x){kegg.path[[x]]}))
  
  
  kegg.p <- keggID
  kegg <- KEGGREST::keggGet(kegg.p)[[1]]$GENE
  kegg.ez <- kegg[seq(1,length(kegg), by = 2)]
  kegg.sym <- AnnotationDbi::select(orgdb, keys
                                    = kegg.ez, keytype = "ENTREZID", columns = "SYMBOL")[,2]
  
  return(kegg.sym)
}

