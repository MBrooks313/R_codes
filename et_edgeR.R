et_edgeR <- function(dge, expt.nm, fc, fdr, res.path = "../data/processed/"){
  
  #############################################
  # Written by MJB on Sept 19th, 2019, last modified on 
  # dge = dge list object of data 
  # expt.nm = expt name to be used in return and exported tables
  # fc = fold change to filter for significance
  # fdr = false discovery rate to filter for significance
  # res.path = path for writing results tables
  #############################################
  
  #############################################
  # Change log 
  # 
  # ***Testing variables at bottom 
  #############################################
  
  library(edgeR)
  library(plyr); library(dplyr)
  library(tibble)
  
  # DGElist object for the samples needed
  y <- dge
  
  # Estimate dispersions needed for exact test
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  
  # Perform DE exact test
  et <- exactTest(y)
  et$table$FDR <- p.adjust(et$table$PValue, method = "fdr")
  
  # Filter for significance
  et.sig <- et$table %>%
    rownames_to_column('gene') %>% 
    filter(FDR < fdr, abs(logFC) >= fc) %>% 
    column_to_rownames('gene')
  
  # Merge results tables
  et_results <- list(all = data.frame(merge(et$genes, et$table, by = 0), row.names = 1),
                 sig = data.frame(merge(et$genes[row.names(et.sig),], et.sig, by = 0), row.names = 1))
  
  # Export results
  write.table(et_results$all, 
              paste0(res.path, "/DE.et.", expt.nm, ".all.tsv"),
              sep = "\t", quote = F, col.names = NA)
  write.table(et_results$sig, 
              paste0(res.path, "/DE.et.", expt.nm, ".sig.tsv"),
              sep = "\t", quote = F, col.names = NA)
  
  # Generate results summaries
  results <- decideTests(et, p.value=fdr, lfc = fc)
  d <- summary(results)
  
  # Print out results summaries
  # print(paste0("FC: ", fc, " , FDR: ", fdr))
  # print(d)
  
  #Return results
  return(assign(paste0("et_results_", expt.nm), et_results, envir = .GlobalEnv))
}


# Testing variables #################
# expt.nm <- "grp_test"
# dge <- gene.dge_in6
# fc <- 0.585
# fdr <- 0.1
# res.path = "../data/processed/"
####################################
