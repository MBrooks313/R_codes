##############################
# Gets the stats from kallisto outputs


kalStats <- function(kal_results=getwd()) {
  #Modified and tested by Matthew J. Brooks on June 2nd, 2019
  #
  #kal_results: directory containing the kallito results directories, defaults to current woring directory
  
  #Set results path and get all results files
  samps <- dir(path=kal_results ,pattern=".+")
  files <- file.path(kal_results, samps,"run_info.json")
  
  #Initialize variabels
  samp.name <- c()
  fastq <- c()
  align <- c()
  perc <- c()
  
  #Loop for getting all relavent parts of the run_info.json
  for (i in 1:length(files)){
    file.data <- readLines(files[i])
    samp.name <- c(samp.name, sub("pf.+", "", files[i]))
    fastq <- c(fastq, as.numeric(sub(".+\\s(\\d+)\\,", "\\1", file.data[4])))
    align <- c(align, as.numeric(sub(".+\\s(\\d+)\\,", "\\1", file.data[5])))
    perc <- c(perc, as.numeric(sub(".+\\s(.+)\\,", "\\1", file.data[7])))
  }
  
  #Make data.frame of the results
  stats <- data.frame(row.names=samp.name, Total=fastq, Align=align, Percentage = perc)
  
  #Export the results
  assign("kal.stats", stats ,envir=.GlobalEnv)
}
