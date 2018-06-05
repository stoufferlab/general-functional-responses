
# read in the simplify the raw data
datadir <- 'Crowley_1989'
filename <- 'Crowley_1989.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Pred','Prey','FeedingRate.Total.Mean','FeedingRate.Total.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
