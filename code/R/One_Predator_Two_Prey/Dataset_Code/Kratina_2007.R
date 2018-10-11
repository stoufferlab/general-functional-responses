
# read in the simplify the raw data
datadir <- 'Kratina_2007'
filename <- 'Kratina_2007.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Pred", "Prey1","Prey2", "Prey1Nconsumed.Total","Prey2Nconsumed.Total",'Time')]
colnames(d) <- c('Npredator',"Nprey1","Nprey2","Nconsumed1","Nconsumed2",'Time')
