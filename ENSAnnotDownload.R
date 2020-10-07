##############################
#This is the download of the gene and transcript level annotation from Ensembl.


ens_nnot <- function(species, release=NULL){
  
  
  ##################################
  # This was written on Oct 27th, 2018 by MJB, last modified on Oct 8th, 2020
  # ENSAnnotDownload)v0.3.1.R
  #
  # USAGE: annot <- ens_annot(species = "Mus musculus", release = "98")
  #
  # species <- species name for annotation, ex. "Homo sapiens"
  # release <- Ensembl release number, ex. "98", default is NULL which takes the current release
  ##################################

  
  
  library(biomaRt)
  
  # Get current annotation version if none stated
  if (is.null(release)){
    annot.ver = listEnsemblArchives()[which("*" == listEnsemblArchives()$current_release),]$`url`
  } else {
    annot.ver = listEnsemblArchives()[which(release == listEnsemblArchives()$version),]$`url`
  }
  
  #Get the database needed
  # listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host = annot.ver)) #Use to see which species annotation is available 
  spec <- strsplit(species, split = " ")
  spec <- paste0(tolower(gsub("(.).+", "\\1", spec[[1]][1])), 
                 spec[[1]][2], 
                 "_gene_ensembl")
  
  ens <- useMart("ENSEMBL_MART_ENSEMBL", dataset = spec, host = annot.ver)
  
  
  #Pull out transcript annotation
  ann.t <- data.frame(getBM(attributes = c("ensembl_transcript_id",
                                                 "ensembl_gene_id",
                                                 "chromosome_name",
                                                 "strand",
                                                 "transcript_start",
                                                 "transcript_end",
                                                 "transcript_length",
                                                 "external_gene_name",
                                                 "description",
                                                 "external_transcript_name",
                                                 "transcript_biotype",
                                                 "transcript_source"), mart = ens))
  print(paste0('Duplicated Transcripts: ', table(duplicated(ann.t$ensembl_transcript_id))[2]))
  ann.t <- data.frame(ann.t, row.names = 1)
  
  #Pull out gene annotation
  ann.g <- data.frame(getBM(attributes = c("ensembl_gene_id",
                                                "chromosome_name",
                                                "strand",
                                                "start_position",
                                                "end_position",
                                                "external_gene_name",
                                                "description",     
                                                "transcript_count",
                                                "gene_biotype",
                                                "source"), mart = ens))
  print(paste0('Duplicated Genes: ', table(duplicated(ann.g$ensembl_gene_id))[2]))
  ann.g <- data.frame(ann.g, row.names = 1)
  
  # Get version information for output
  ann.v <- listEnsemblArchives()[which(annot.ver == listEnsemblArchives()$url),]$`version`
    
  # Return
  ann <- list(trans = ann.t, gene = ann.g, version = ann.v)
  
}




# TESTING #########################
# listEnsemblArchives() #Use to get which ENSEMBL annotation version needed
# species <- "Homo sapiens" #testing
# annot.ver <- "94" # testing
#
# annot <- ens_annot(species = "Mus musculus", release = "94")
##################################
