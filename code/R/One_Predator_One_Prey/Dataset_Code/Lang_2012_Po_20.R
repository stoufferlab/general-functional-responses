
# read in the raw data
datadir <- 'Lang_2012'
filename <- 'Lang_2012_Poe_20C.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("PredNo", "Initial", "Killed",'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')