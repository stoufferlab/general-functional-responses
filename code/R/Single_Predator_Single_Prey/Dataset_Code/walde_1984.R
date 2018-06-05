
# read in the simplify the raw data
datadir <- 'Walde_1984'
filename <- 'Walde_1984.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c("Pred","Prey",'Captures.Total.mean','Captures.Total.se','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
