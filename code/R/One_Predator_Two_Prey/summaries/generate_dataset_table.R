
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we need bundle_fits and study.info
source('../../lib/plot_coefs.R')
source('../../lib/study_info.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredTwoPrey_fits')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out <- dim(0)
for(i in 1:length(ffr.fits)){

  # grab info from the google doc
  this.study <- study.info(ffr.fits[[i]]$study.info$datadir)
  
  # publication
  pub <- gsub('_',' ', this.study$dataname)
  
  # citation
  cite <- this.study$cite
  
  # dataset
  datasetsName <- ffr.fits[[i]]$study.info$datasetName
  datasetsName <- gsub('_',' ', datasetsName)
  if(pub==datasetsName){
    datasetsName <- '-'
  }else{
    datasetsName <- gsub(paste0(pub,' '),'', datasetsName)
  }
  
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
  if("data.Nconsumed1.mean" %in% names(ffr.fits[[i]]$study.info)){
    SS <- sum(ffr.fits[[i]]$study.info$data.n)
  }else{
    SS <- length(ffr.fits[[i]]$study.info$data.Nconsumed1)
  }
  
  # wrap it all up
  out <- rbind(out, 
         c(cite,
           datasetsName,
           orig,
           datasource,
           datafigtablesource,
           SS,
           repl,
           pred))
  print(paste(i," of ",length(ffr.fits)))
  
}

tab <- data.frame(out)
colnames(tab) <- c('Study',
                   'Dataset',
                   'Raw',
                   'Type',
                   'Source',
                   'Nobs',
                   'Replaced',
                   'Consumer')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reorder columns
tab <- tab[,c('Study','Dataset','Nobs','Replaced','Consumer',
              'Raw','Type','Source')]

# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredTwoPrey_tables/')
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
