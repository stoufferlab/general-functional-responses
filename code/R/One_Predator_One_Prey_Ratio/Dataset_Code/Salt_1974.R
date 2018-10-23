
# read in the simplify the raw data
datadir <- 'Salt_1974'
filename <- 'Salt_1974.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','FeedingRate.Total','FeedingRate.Total.SE','n','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')
