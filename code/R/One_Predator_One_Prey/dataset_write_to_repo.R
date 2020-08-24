rm(list = ls())
# set to FALSE if you want to match messages in real time 
# or TRUE to have them silently saved to file instead.
sinkMessages <- TRUE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify where the data files are located
dropboxdir <- switch(
  Sys.getenv("LOGNAME"),
  stouffer = '~/Dropbox/Projects/GenFuncResp/Data',
  marknovak = '~/Dropbox/Research/Projects/GenFuncResp/Data'
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# utility functions
source('../lib/study_info.R')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# folder in which to place all datasets
mainDir <- '../../../data/One_Predator_One_Prey/'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in the table of dataset details
dataset_details <- read.table(
  '../../../data/dataset_details.csv'
)

# master list of datasets
datasets <- list.files('./Dataset_Code', pattern=".R$", full.names=TRUE, include.dirs=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(i in 1:length(datasets)){
  
  # loads the data into data frame 'd' and specifies data-specific parameters
  source(datasets[i])
  
  # grab info about experimental design, etc
  this.study <- study.info(
    dataset_details,
    datadir,
    "One_Predator_One_Prey"
  )

  # has this dataset been okayed to post to the repo
  ok2post <- this.study$ok2post
  if(is.na(ok2post)){ok2post<-FALSE}
  
  if(ok2post==TRUE){
    
    # study folder
    subDir <- this.study$dataname
    
    # Create study folder (if it doesn't already exist)
    dir.create(file.path(mainDir, subDir))
    
    datasetsName <- sub('*.R$','', sub('*./Dataset_Code/','', datasets[i]))
    
    write.csv(d,paste0(mainDir, subDir,'/', datasetsName,'.csv'))
    
    print(paste(datasetsName, "exported"))
    
  }
  
  print(paste(i," of ",length(datasets)))
}


