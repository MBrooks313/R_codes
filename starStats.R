##############################
# Gets the stats from star log outputs

starStats <- function(star_results=getwd()) {
  #Modified and tested by Matthew J. Brooks on June 2nd, 2019
  #
  #star_results: directory containing the star results directories containing the Log.final.out file, defaults to current woring directory
  
  samps <- dir(path=star_results ,pattern=".+star$")
  files <- file.path(star_results, samps,"Log.final.out")
  #files <- (paste(path,"/Log.final.out",sep=''))
  
  #samps <- dir(path=kal_results ,pattern=".+")
  #files <- file.path(kal_results, samps,"run_info.json")
  
  samp.name <- c()
  fastq <- c()
  uniq <- c()
  mult <- c()
  len <- c()
  uniq.perc <- c()
  mult.perc <- c()
  unmap.short.perc <- c()
  unmap.other.perc <- c()
  
  for (i in 1:length(files)){
    file.data <- readLines(files[i])
    samp.name <- c(samp.name, sub("pf.+", "", files[i]))
    fastq <- c(fastq, as.numeric(sub(".+\t", "", file.data[6])))
    uniq <- c(uniq, as.numeric(sub(".+\t", "", file.data[9])))
    mult <- c(mult, as.numeric(sub(".+\t", "", file.data[24])))
    len <- c(len, as.numeric(sub(".+\t", "", file.data[7])))
    uniq.perc <- c(uniq.perc, sub(".+\t", "", file.data[10]))
    mult.perc <- c(mult.perc, sub(".+\t", "", file.data[25]))
    unmap.short.perc <- c(unmap.short.perc, sub(".+\t", "", file.data[30]))
    unmap.other.perc <- c(unmap.other.perc, sub(".+\t", "", file.data[31]))
  }
  stats <- data.frame(row.names=samp.name, Total=fastq, UniqueMap=uniq, MultiMap=mult)
  percs <- data.frame(row.names=samp.name, Ave.Frag.Len=len, UniqueMap=uniq.perc, MultiMap=mult.perc, 
                      UnMap.Short=unmap.short.perc, UnMap.Other=unmap.other.perc)
  assign("star.stats", stats ,envir=.GlobalEnv)
  assign("star.percs", percs ,envir=.GlobalEnv)
}