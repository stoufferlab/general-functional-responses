
# read in the simplify the raw data
datadir <- 'Wong_2005'
filename <- 'Wong_2005_ss.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey1.m','Prey1Eaten','Prey2.s','Prey2Eaten','Time')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')


