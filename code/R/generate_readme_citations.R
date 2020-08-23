
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Hmisc) # for LaTeX table export
options(xdvicmd='open')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we need bundle_fits and study.info
source('lib/plot_coefs.R')
source('lib/study_info.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in the dataset-specific fits and combine
ffr.fits_OnePrey <- bundle_fits('../../results/R/OnePredOnePrey_fits')
ffr.fits_TwoPrey <- bundle_fits('../../results/R/OnePredTwoPrey_fits')
ffr.fits <- c(ffr.fits_OnePrey, ffr.fits_TwoPrey)

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
    
    # paper citation
    cite <- this.study$cite
    
    # how to cite the data
    datacitation <- this.study$datacitation
    
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
    
    # wrap it all up
    out <- rbind(out, 
                 c(cite,
                   datacitation))
    }
  print(paste(i," of ",length(ffr.fits)))
  
}

tab <- data.frame(out)
colnames(tab) <- c('Study citation',
                   'Data source (repository) citation')


# Alphabetical order
tab <- tab[order(tab$Study),]
# Reduce down to unique studies (when datasetname isn't provided)
tab <- unique(tab)

# Export to LaTeX
wd <- getwd()
setwd('../../data/readme_citations/')
latex(
  tab,
  file='readme_citations_tab.tex',
  rowname=NULL, 
  na.blank=TRUE, 
  longtable=FALSE,
  lines.page=100,
  caption="Please cite appropriately! "
)
setwd(wd)