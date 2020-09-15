################################
# Runs functional gene enrichment using gProfileR

fge <- function(dat, gpout, catnm=10, paths = c("rea", "keg", "BP", "CC"), dat_col, sym_col, samp_ord, expt.nm, res.path = "../data/processed/"){
  
  #############################################
  # Written by MJB on Sept 20th, 2019, last modified on 
  # 
  # dat = expression data.frame to be used for heatmap
  # gpout = gProfileR output
  # catnm = number of categories to examine
  # paths = pathway enrichments to examine, c("rea", "keg", "BP", "CC"), 
  # dat_col = columns in expression data data.frame to be used in heatmap
  # sym_col = coulmn in expression data.frame with the gene symbol
  # samp_ord = sample order x-axis to be ued in heatmap
  # expt.nm = experiment name to be used in the subtitle and export file
  # res.path = base path for exports
  #############################################
  
  
  # Check if needed functions are loaded
  func = sum(existsFunction("GO_leaf"),
             existsFunction("dot_go"),
            existsFunction("bar_go"),
            existsFunction("heatmap_go"))
  
  if (func != 4){
    print("Please ensure needed functions are loaded: GO_leaf, dot_go, bar_go, and heatmap_go.")
    break
  }
  
  # Loop if the enrichment results are not empty 
  if (dim(gpout)[1] > 0){
    
    # Prep directories
    dir.gProf <- paste0(res.path, "/gProf_", expt.nm)
    dir.create(dir.gProf)
    dir.heat <- paste(dir.gProf, "Heatmaps", sep = "/")
    dir.create(dir.heat)

    # Set paths to examine
    path = paths
    
    # Loop for each type of enrichment
    for (k in path){
      dat_sub <- gpout %>%
        filter(
          domain == k
        )
      
      # GO redundancy reduction
      if (k %in% c("BP", "CC", "MF")){
        dat_sub <- GO_leaf(dat_sub, id.col = 9, ont = k)
      }
      
      # Make dot plot
      dot_go(dat = dat_sub, top = catnm, grp = "query.number", desc = "term.name", 
             p.val = "p.value", ratio = "recall", filenm = paste(dir.gProf, "/", expt.nm, "-", k, "-", catnm, sep = ""))
      
      # Loop for each cluster
      for (i in unique(dat_sub$query.number)) {
        dat_tmp <- dat_sub %>%
          filter(
            query.number == i
          )
        
        # Generate bar plot
        bar_go(dat = dat_tmp, top = catnm, desc = "term.name", 
               p.val = "p.value", genenum = "term.size", filenm = paste(paste(dir.gProf, "/", expt.nm, "-", k, sep = ""), "_", i, "-", catnm, sep = ""))
        
        # Generate heatmap for each pathway in bar plot
        print(i)
        heatmap_go(expdat = dat, godat = dat_tmp, sym_col = sym_col, dat_col = dat_col, samp_ord = samp_ord, top = catnm,
                   filenm = paste(dir.heat, "/", expt.nm, "-", k, "-", i, sep = ""))
      }
    }
  } else {
    print("The functional gene enrichment results supplied is empty!")
  }
}


# Testing variables #################
# dat <- data.frame(merge(gene_dge_e1$genes[row.names(et.sig),], gene_dge_e1$lcpm[row.names(et.sig),], by = 0), row.names = 1)
# dat <- dat[dat$chromosome_name %in% chrs,]
# gpout = gpout_all
# catnm = 20
# dat_col = 10:15
# samp_ord = c(1:5)
# sym_col = 5
# paths = c("rea", "keg", "BP")
# expt.nm = "Expt_test"
# res.path = "/Volumes/MB32/CutNRun/Mouse/Analysis/20200309_16-29/200604/data/processed/"
####################################
