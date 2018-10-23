
# read in the raw data
datadir <- 'vonWesternhagen_1976'
filename <- 'vonWesternhagen_1976_Fig2_4hrs.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Predators','Prey','Total.Eaten.mean','SE','n','Hours')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n", "Time")
