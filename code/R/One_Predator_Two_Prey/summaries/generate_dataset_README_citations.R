
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
  
  # has this dataset been okayed to post to the repo
  ok2post <- this.study$ok2post
  if(is.na(ok2post)){ok2post<-FALSE}
  
  if(ok2post){
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
    
    # where the data came from (for extracted data)
    datafigtablesource <- this.study$datafigtablesource
    
    # how to cite the data
    datacitation <- this.study$datacitation
  
    # wrap it all up
    out <- rbind(out, 
                 c(cite,
                   datasetsName,
                   datacitation,
                   datafigtablesource))
    }
  print(paste(i," of ",length(ffr.fits)))
  
}

tab <- data.frame(out)
colnames(tab) <- c('Study',
                   'Dataset',
                   'Citation',
                   'Source')

# Export to LaTeX
wd <- getwd()
setwd('../../../../results/R/OnePredTwoPrey_tables/')
latex(
  tab,
  file='OnePredTwoPrey_dataset_citations4README.tex',
  rowname=NULL, 
  na.blank=TRUE, 
  longtable=TRUE,
  lines.page=100,
  caption="
  Multi-resource datasets
  "
)
setwd(wd)
