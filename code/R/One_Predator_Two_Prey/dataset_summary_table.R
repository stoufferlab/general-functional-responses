rm(list = ls())
# set to FALSE if you want to match messages in real time 
# or TRUE to have them silently saved to file instead.
sinkMessages <- TRUE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
  Sys.getenv("LOGNAME"),
  stouffer = '../../../dropbox_data/Data',
  marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions
source('../lib/study_info.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# master list of datasets
datasets <- list.files('./Dataset_Code', full.names=TRUE, include.dirs=FALSE)

# remove template files which don't actually read data
datasets <- grep("template",datasets,invert=TRUE,value=TRUE)

# remove zzz files which are placeholders while a dataset is being cleaned/incorporated
datasets <- grep("zzz",datasets,invert=TRUE,value=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out <- dim(0)
for(i in 1:length(datasets)){
  
  # loads the data into data frame 'd' and specifies data-specific parameters
  source(datasets[i])
  
  # grab info from the google doc
  this.study <- study.info(datadir)
  
  # publication
  pub <- gsub('_',' ', this.study$dataname)
  
  # citation
  cite <- this.study$cite
  
  # dataset
  datasetsName <- sub('*.R$','', sub('*./Dataset_Code/','', datasets[i]))
  datasetsName <- gsub('_',' ', datasetsName)
  if(pub==datasetsName){
    datasetsName <- '-'
  }else{
    datasetsName <- gsub(paste0(pub,' '),'', datasetsName)
  }
  
  # used in analyses
  used <- ifelse(grepl("H|R", this.study$runswith),'Yes','No')
  
  # replacement
  repl <- ifelse(this.study$replacement,'Yes','No')
  
  # original data, or means and intervals (i.e. bootstrapped)
  orig <- ifelse(this.study$bootstrap,'No','Yes')
  
  # how we got the data
  datasource <- this.study$datasource
  
  # where the data came from (for extracted data)
  datafigtablesource <- this.study$datafigtablesource
  
  # pred/parasite
  pred <- ifelse(this.study$predator,'Predator','Parasitoid')
  
  # sample size
  if("Nconsumed1.mean" %in% colnames(d)){
    SS <- sum(d$n)
  }else{
    SS <- nrow(d)
  }
  
  # wrap it all up
  out <- rbind(out, 
         c(cite,
           datasetsName,
           used,
           orig,
           datasource,
           datafigtablesource,
           SS,
           repl,
           pred))
  print(paste(i," of ",length(datasets)))
  
}

tab <- data.frame(out)
colnames(tab) <- c('Study',
                   'Dataset',
                   'Used',
                   'Raw',
                   'Type',
                   'Source',
                   'Nobs',
                   'Replaced',
                   'Consumer')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Remove unused datasets at end of table
unused <- which(tab$Used=='No')
tab <- tab[-unused,]
tab$Used <- NULL

# reorder columns
tab <- tab[,c('Study','Dataset','Nobs','Replaced','Consumer',
              'Raw','Type','Source')]

# Export to LaTeX
wd <- getwd()
setwd('../../../results/R/OnePredTwoPrey_tables/')
  latex(
    tab,
    file='OnePredTwoPrey_datasets.tex',
    label='table:1pred2preydatasets', 
    rowname=NULL, 
    na.blank=TRUE, 
    longtable=TRUE,
    lines.page=100,
    caption="
      A summary of multi-species resource dependence datasets.
      ``Dataset'' refers to the specific experiment from the study, and `-' implies there was only one dataset available.
      ``Nobs'' indicates the sample size per resource consumed.
      ``Replacement'' refers to whether the consumed resources were replaced during the study, which dictated our use of a binomial versus a Poisson likelihood.
      ``Consumer'' refers to whether the consumer was a predator or a parasitoid.
      ``Raw'' refers to whether we were able to use the raw data at the level of each treatment replicate, or whether we instead used means and associated uncertainty intervals to produce bootstrapped datasets.
      ``Type'' refers to whether the data was provided to us by the author, was obtained from an online repository, or was extracted from the publication.
      ``Source'' refers to the figures and tables from which the data where extracted.
    "
  )
setwd(wd)


# Table of data sets by confidence interval type
# load('../../../../results/R/OnePredOnePrey_fits_profiled/ffr.fits.prof.AA.Rdata')
# labels <- unlist(lapply(ffr.fits, function(x) x$study.info$datasetName))
# method <- unlist(lapply(ffr.fits, function(x) x$profile$method))
# table(method)
# CI.method <- data.frame(DataSet=labels, Method=method)

