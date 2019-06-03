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
  
  # bootstrapped
  boot <- ifelse(this.study$bootstrap,'No','Yes')
  
  # pred/parasite
  pred <- ifelse(this.study$predator,'Predator','Parasitoid')
  
  # Sample size
  if("Nconsumed.mean" %in% colnames(d)){
    SS <- sum(d$n)
  }else{
    SS <- nrow(d)
  }
  
  # wrap it all up
  out <- rbind(out, 
         c(pub, datasetsName, used, repl, boot, pred, SS))
  
}

tab <- data.frame(out)
colnames(tab) <- c('Publication','Dataset','Used','Replacement','Raw data','Consumer','Sample size')

# export to LaTeX
# too long for single LaTeX table, so split up
tab1 <- tab[1:30,]
tab2 <- tab[31:60,]
tab3 <- tab[61:nrow(tab),]

wd <- getwd()
setwd('../../../results/R/OnePredOnePrey_figs/')

  latex(tab1,file='OnePredOnePrey_datasets_1.tex',label='table:datasets', rowname=NULL, na.blank=TRUE, caption='Summary of considered datasets.')
  latex(tab2,file='OnePredOnePrey_datasets_2.tex',label='table:datasets', rowname=NULL, na.blank=TRUE, caption='Summary of considered datasets continued.')
  latex(tab3,file='OnePredOnePrey_datasets_3.tex',label='table:datasets', rowname=NULL, na.blank=TRUE, caption='Summary of considered datasets continued.')

setwd(wd)
