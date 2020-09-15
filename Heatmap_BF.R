hm.clust_BF <- function(dat, col.pal=NULL, clst=NA, sidebar=FALSE, sidebar.pal=NULL, sidebar.plc=0) {
  
  #############################################
  # Written by Ben Fadl, modified by Matthew J. Brooks
  # v0.1.0 tested on June 10th, 2019 by MJB
  #
  # zscore : data.frame of the Z-score data to be plotted
  # col.pal : color palette for plot
  # sidebar : add sidebar of colors to demarcate clusters (TRUE/FALSE)
  # sidebar.pal : color palette for sidebar
  # sidebar.plc : sidebar placement after column number, on left by default (0)
  # clst : cluster list for the zscore data.frame
  #############################################
  
  #Packages required
  require("ggplot2")
  
  #Set color scheme for heatmap
  if (is.null(col.pal)){
    col.hm <- colorRampPalette(c("#023858","#3690c0", "#e0f3f8","#fff7bc","#fec44f","#ec7014",
                                    "#fc4e2a","#bd0026","#800026"))(n=1000)
  } else {
    col.hm <- colorRampPalette(col.pal)(n=1000)
  }
  
  #Set color scheme for sidebar
  if (is.null(sidebar.pal)) {
    col.side <- rainbow(n = length(clst))
  } else {
    col.side <- colorRampPalette(sidebar.pal)(n = length(clst))
  }
  
  #Make sure the data is a data.frame
  dat <- data.frame(dat)
  
  #Transform wide-format data.frame to tall-format
  m <- data.table::melt(tibble::rownames_to_column(dat, var = "rn"), idvars = "rn")
  m$rn <- factor(m$rn, levels = rownames(dat))
  
  #Make the heatmap plot
  g <- ggplot(m, aes(x = variable, y = rev(rn))) +
    geom_tile(aes(fill = value), colour = "transparent") +
    scale_fill_gradientn(colors = col.hm,
                         name = "Z-score",
                         limits = c(-3, 3)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank()
    )
  
  if (sidebar == TRUE){
    
    ##Make column for demarcation of clusters
    #Initial settings
    cluster <- 1:length(clst)
    y.start <- dim(dat)[1] + .5
    
    #Placement along x-axis
    x <- sidebar.plc + 0.5
    
    #Add cluster colorbar to plot
    for (c in cluster) {
      c.size <- length(clst[[c]])
      y.end <- y.start - c.size
      g <-
        g + 
        #Add border around sidebar
        geom_segment(
          x = x,
          xend = x,
          y = y.start,
          yend = y.end,
          size = 7,
          color = "white"
        ) + 
        #Add sidebar
        geom_segment(
          x = x,
          xend = x,
          y = y.start,
          yend = y.end,
          size = 6,
          color = col.side[c]
        )
      y.start <- y.end
    }
  }
  
  #Return plot
  print(g)

}
