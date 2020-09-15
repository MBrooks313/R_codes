
heatmap_go <- function(expdat, godat, sym_col, dat_col, samp_ord=NA, top, filenm) {
  
  #############################################
  # Written by MJB on Oct 24th, 2018, last modified on Oct 25th, 2018
  # expdat = expression data.frame with annotation
  # godat = gProfiler output
  # sym_col = column number for the gene symbol from exprdat
  # dat_col = column number for the data to be plotted
  # samp_ord = sample order for factor of column names
  # top = top number of pathways to plot
  # filenm = prefix for file name generation
  #############################################
  
  
  library(tidyr)
  library(gtools)
  library(scales)
  
  mb.colors <- colorRampPalette(c("#023858","#3690c0", "#e0f3f8","#fff7bc","#fec44f","#ec7014",
                                  "#fc4e2a","#bd0026","#800026"))(n=1000)
  jet <- c("#7F0000","red","#FF7F00","yellow","white","cyan", "#007FFF", "blue","#00007F")
  col.HC <- colorRampPalette(rev(jet))(1000)
  
  id_col = "intersection"
  id_term <- "term.name"
  p_val <- "p.value"
  
  godat <- godat[order(godat[,`p_val`]),] #Order by p value
  
  if (is.na(top)){
    godat = godat
  } else if (top > dim(godat)[1]){
    godat = godat
  } else {
    godat = godat[1:top,]
  }
  
  
  #Loop through godata for each row
  for (i in 1:dim(godat)[1]){
    
    #Get gene IDs
    ids <- unlist(strsplit(godat[i,`id_col`], split = ","))
    
    #Get term.name
    term <- paste0(unlist(strsplit(godat[i,`id_term`], split = " ")), collapse = "-")
    term <- substring(term, 1, 30)
    term <- gsub("\\/", "", term)
    term <- gsub(",", "", term)
    
    #Subset expdata for gene ID
    #edat <- expdat[row.names(expdat) %in% ids,]
    edat <- expdat[toupper(expdat$external_gene_name) %in% ids,]
    edat <- edat[!duplicated(edat$external_gene_name),]
    
    #Tidy table
    dat.tidy <- edat[,c(sym_col, dat_col)] %>%
      gather(Sample, CPM, -external_gene_name)
    
    #Factor sort sample names and gene name
    if (is.na(samp_ord)){
      dat.tidy$Sample <- factor(dat.tidy$Sample, levels = mixedsort(levels(factor(dat.tidy$Sample))))
    } else {
      dat.tidy$Sample <- factor(dat.tidy$Sample, levels = mixedsort(levels(factor(dat.tidy$Sample)))[samp_ord])
    }
    
    dat.tidy$external_gene_name <- factor(dat.tidy$external_gene_name)
    
    #PCA method of ordering
    pca <- prcomp(edat[,dat_col], center = T, scale. = F)
    pca.order <- order(pca$x[,1])
  
    
    #Check to make ensure plot is decreasing from top to bottom
    dat.tmp <- dat.tidy
    dat.tmp$external_gene_name <- factor(dat.tmp$external_gene_name, levels = dat.tmp$external_gene_name[pca.order], ordered = T)
    
    first.mean <- dat.tmp %>%
      filter(external_gene_name == levels(external_gene_name)[1]) %>% summarize(mean(CPM))
    
    last.mean <- dat.tmp %>%
      filter(external_gene_name == levels(external_gene_name)[length(pca.order)]) %>% summarize(mean(CPM))
     
    if (first.mean > last.mean) {
      dat.tidy$external_gene_name <- factor(dat.tidy$external_gene_name, levels = dat.tidy$external_gene_name[rev(pca.order)], ordered = T)
    } else {
      dat.tidy$external_gene_name <- factor(dat.tidy$external_gene_name, levels = dat.tidy$external_gene_name[pca.order], ordered = T)
    }
    
    #Plot
    pdf(paste(filenm, "_", term, "_hm.pdf", sep=""), 
        width = dim(unique(dat.tidy[2]))[1] / 2 + 1, 
        height = dim(unique(dat.tidy[1]))[1] / 5 + 1, 
        useDingbats = F)
    
    print(ggplot(dat.tidy, aes(Sample, external_gene_name)) +
      geom_tile(aes(fill = CPM), colour = "white") + 
      #scale_fill_gradient(high = "white", low = "steelblue", name = "CPM (log2)") +
      scale_fill_gradientn(colours = mb.colors, name = "CPM (log2)", limits = c(0,16), oob=squish) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(), axis.title.y = element_blank()))
    dev.off()
    
  }
    
}


#Testing variables #################
# godat = dat_tmp
# expdat = gene.ave.lcpm.mstr
# sym_col = 5
# dat_col = 11:22
# top = 10
# filenm = "SC-kegg"
####################################
