################################
# Makes gene track plot from bigwig files


# Genetrack function
gene_track_MB <- function(goi, in_dir, dat_files, track_nm, samp_grp=NULL, 
                          out_dir, y_lim=NULL, seq_up=2000, seq_dwn=2000,
                          track_col="darkslategrey", ucsc_style = FALSE){
  
  #######################################
  # Written and tested by Matthew J Brooks on Sept 21st, 2020
  #
  # goi :  genes of interest, i.e gene symbols
  # in_dir : directory where big wigs are found
  # dat_files : list of bigwig files
  # track_nm :  vecotr of dat file short names for plot
  # samp_grp : factor of dat file groups or NULL (default), when !is.null() scales and colors are unified per group
  # out_dir : output directory for plot
  # y_lim : NULL (default) vector of min and max values for all samples in plot. when !is.null then each sample or group will define y_lim
  # seq_up :  distance upstream of gene to show in plot
  # seq_down :  distance downstream of gene to show in plot
  # track_col :  colors used for tracks samples (all the same) or based on group (group named vector)
  # ucsc_style : whether to have UCSC chrmosome name style, FALSE (default)
  ####################################### 
  
  
  
  # Load libraries
  library(AnnotationHub)
  library(ensembldb)
  library(wesanderson)
  library(Gviz)
  library(rtracklayer)
  library(trackViewer)
  library(org.Mm.eg.db)
  
  
  # Load annotation db if necessary
  if (!exists("ensDb")){
    
    ah <- AnnotationHub()
    
    ## Query for all available EnsDb databases
    ahdb <- query(ah, pattern = c("EnsDb", "Mus musculus", 99))
    
    ensDb <- ahdb[[1]]
    
  }
  
  
  # Make output directory if necessary
  dir.create(track_dir, recursive = T)
  
  #Sets GVIZ to 'ENSEMBL' or 'UCSC' chromosome naming sytle 
  options(ucscChromosomeNames = ucsc_style) 
  
  # Restrict results to main chromosomes
  chrs <- c(1:21, "X", "MT")
  
  # Track colors
  gene_col = wes_palette("BottleRocket1", 7, "discrete")[c(2)]
  name_col = wes_palette("BottleRocket1", 7, "discrete")[c(7)]
  
  
  
  #--------------------------------------------------------------------------#
  
  # Check if gene sysmbol exists
  if (exists(goi, org.Mm.egSYMBOL2EG)){
    
    entrez <- get(goi, org.Mm.egSYMBOL2EG)
    
    # Filter for only the protein coding genes
    goi_gr_ens <- getGeneRegionTrackForGviz(ensDb, 
                                            filter = list(TxBiotypeFilter("protein_coding"),
                                                          EntrezFilter(entrez),
                                                          SeqNameFilter(chrs)))
    if (length(goi_gr_ens) == 0){
      goi_gr_ens <- getGeneRegionTrackForGviz(ensDb, 
                                              filter = list(GeneNameFilter(goi),
                                                            SeqNameFilter(chrs)))
    }
    
    
    # Set upstream and downstream buffer region
    upstr <- seq_up
    dwnstr <- seq_dwn
    
    
    
    #------------------------#
    # Gene axis and region track
    
    # Gene axis
    gtrack <- GenomeAxisTrack()
    
    # Gene region 
    grtrack <- GeneRegionTrack(goi_gr_ens, name = goi, stacking = 'dense')
    displayPars(grtrack) <- list(fill = gene_col, size = 2,
                                 col = gene_col,
                                 col.title = name_col)
    
    
    
    #------------------------#
    # Get Data Tracks
    dTrack <- list()
    for (j in 1:length(dat_files)){
      
      # Get track data
      dTrack[j] <- DataTrack(range = paste0(in_dir, dat_files[j]),
                             genome = "mm10",
                             type = "histogram",
                             name = track_nm[j],
                             window = -1,
                             chromosome = seqlevels(goi_gr_ens))
    }
    
    
    
    #------------------------#
    # Get y-axis ranges for each data track
    y <- plotTracks(dTrack,
                    from = start(range(goi_gr_ens)) - upstr, 
                    to = end(range(goi_gr_ens)) + dwnstr)
    
    
    
    #------------------------#
    # Display parameters for dTracks
    
    #--------#
    # Samples not grouped
    if (is.null(samp_grp)){
      
      
      # y_lim not specified
      if (is.null(y_lim)){
        
        for (m in 1:length(dTrack)){
          displayPars(dTrack[[m]]) <- list(size = 2,
                                           fill = track_col, col.histogram = track_col,
                                           col.title = name_col)
        }
        
        
        # y_lim specified    
      } else {
        
        for (m in 1:length(dTrack)){
          displayPars(dTrack[[m]]) <- list(ylim = y_lim, size = 2,
                                           fill = track_col, col.histogram = track_col,
                                           col.title = name_col)
        }
        
      }
      
      
      #--------#  
      # Samples grouped
    } else {
      
      
      # y_lim not specified
      if (is.null(y_lim)){
        
        # Loop for each group
        for (k in levels(samp_grp)){
          
          # Get y_lim per group
          rng <- c()
          for (r in grep(k, samp_grp)){
            
            rng_tmp <- range(y[[r]]@data)
            rng <- c(rng, rng_tmp)
          }
          
          # Set new grouped y_lim
          y_lim <- c(0, max(rng))
          
          # Group defined display params
          for (m in grep(k, samp_grp)){
            grp_track_col <- track_col[k]
            displayPars(dTrack[[m]]) <- list(ylim=y_lim, size = 2,
                                             fill = grp_track_col, col.histogram = grp_track_col,
                                             col.title = name_col)
          }
          
        }
        
        
        
        # y_lim specified  
      } else {
        
        # Loop for each group
        for (k in levels(samp_grp)){
          
          # Group defined display params
          for (m in grep(k, samp_grp)){
            grp_track_col <- track_col[k]
            displayPars(dTrack[[m]]) <- list(ylim=y_lim, size = 2,
                                             fill = grp_track_col, col.histogram = grp_track_col,
                                             col.title = name_col)
          }
          
        }
        
      }
      
      
    }
    
    
    #---------------------------#
    # Generate the output plot
    
    pdf(paste0(out_dir, goi, "_peak_histo.pdf"), width=6, height=5)
    
    plotTracks(c(list(gtrack, grtrack), dTrack),
               from = start(range(goi_gr_ens)) - upstr, 
               to = end(range(goi_gr_ens)) + dwnstr,
               showID=T, shape = "arrow")
    
    dev.off()
    
    
  }
  
  
}





################
#Test variables
goi = "Rho"
in_dir = "/Volumes/MB32/CutNRun/Mouse/Analysis/20200309_16-29/200723/data/interim/bigwig/"
dat_files = list.files(path = in_dir, pattern = ".bw", recursive = TRUE)
dat_files <- dat_files[c(3,6,2,5,
                         9,12,8,11)]
out_dir = "03_GeneTrack_Plots/"
seq_up=2000; seq_dwn=2000
#y_lim = c(0,10)
y_lim = NULL
track_nm <- gsub("(.+)_IDR.+", "\\1", bw_idr)
samp_grp <- factor(gsub("(.+)_.+", "\\1", track_nm))
#samp_grp = NULL
#track_col="darkslategrey"
track_col <- wes_palette("Zissou1", 5, "discrete")[c(1,5)]
names(track_col) <- levels(samp_grp)
ucsc_style = FALSE
################