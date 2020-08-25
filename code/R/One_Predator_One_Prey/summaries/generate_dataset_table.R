
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions
source('../../lib/plot_coefs.R')
source('../../lib/study_info.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in the table of dataset details
dataset_details <- read.csv(
  '../../../../data/dataset_details.csv'
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in the dataset-specific fits into a mega container
ffr.fits <- bundle_fits('../../../../results/R/OnePredOnePrey_fits')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out <- dim(0)
for(i in 1:length(ffr.fits)){
  
  # grab info about experimental design, etc
  this.study <- study.info(
    dataset_details,
    ffr.fits[[i]]$study.info$datadir,
    "One_Predator_One_Prey"
  )

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
  
  # how to cite the data
  datacitation <- this.study$datacitation
  
  # where the data came from (for extracted data)
  datafigtablesource <- this.study$datafigtablesource
  
  # pred/parasite
  pred <- ifelse(this.study$predator,'Predator','Parasitoid')
  
  # sample size
  if("data.Nconsumed.mean" %in% names(ffr.fits[[i]]$study.info)){
    SS <- sum(ffr.fits[[i]]$study.info$data.n)
  }else{
    SS <- length(ffr.fits[[i]]$study.info$data.Nconsumed)
  }
  
  # wrap it all up
  out <- rbind(out, 
         c(cite, 
           datasetsName, 
           orig,
           datasource,
           datafigtablesource,
           datacitation,
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
                   'Citation',
                   'Nobs',
                   'Replaced',
                   'Consumer')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For SampleSizeBias paper
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')
  latex(tab,file='OnePredOnePrey_datasets.tex',
        label='table:datasets', 
        rowname=NULL, 
        na.blank=TRUE,
        longtable=TRUE,
        lines.page=100,
        caption="A summary of used datasets. 
        ``Dataset'' refers to the specific experiment from the study, and ‘-’ implies there was only one dataset available.
        ``Raw'' refers to whether we were able to use the raw data at the level of each treatment replicate, or whether we instead used means and associated uncertainty intervals to produce bootstrapped datasets. 
          ``Type'' refers to whether the data was provided to us by the author, was obtained from an online repository, or was extracted from the publication.
          ``Source'' refers to the figures and tables from which the data where extracted.
        ``Nobs'' indicates the sample size.
        ``Replaced'' refers to the whether consumed prey were replaced during the study (or whether the parasitoid was considered discriminatory or not), which dictated our use of a binomial versus a Poisson likelihood. 
        ")
setwd(wd)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tweak things for HiddenLayers paper
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reorder columns
tab <- tab[,c('Study','Dataset','Nobs','Replaced',
              'Consumer','Raw','Type','Source','Citation')]

# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredOnePrey_tables/')
  latex(
    tab,
    file='OnePredOnePrey_datasets_DBS.tex',
    label='table:1pred1preydatasets', 
    rowname=NULL, 
    na.blank=TRUE,
    longtable=TRUE,
    lines.page=100,
    caption="
      A summary of single-resource consumer dependence datasets.
      ``Dataset'' refers to the specific experiment from the study, and `-' implies there was only one dataset available.
      ``Nobs'' indicates the sample size.
      ``Replaced'' refers to whether the consumed resource was replaced during the study, which dictated our use of a binomial versus a Poisson likelihood.
      ``Consumer'' refers to the whether the consumer was a predator or a parasitoid.
      ``Raw'' refers to whether we were able to use the raw data at the level of each treatment replicate, or whether we instead used means and associated uncertainty intervals to produce bootstrapped datasets.
      ``Type'' refers to whether the data was provided to us by the author, was obtained from an online repository, or was extracted from the publication.
      ``Source'' refers to the figures and tables from which the data where extracted.
    "
  )
setwd(wd)
