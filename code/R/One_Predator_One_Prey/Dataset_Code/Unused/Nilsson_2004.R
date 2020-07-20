
# read in the simplify the raw data
datadir <- 'Nilsson_2004'
filename <- 'Nilsson_2004.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Pred", "Prey", "TotalFeedingRate")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")