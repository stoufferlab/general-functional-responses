
# read in the simplify the raw data
datadir <- 'Lang_2012'
filename <- 'Lang_2012_Pter_10C.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("PredNo", "Initial", "Killed",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
