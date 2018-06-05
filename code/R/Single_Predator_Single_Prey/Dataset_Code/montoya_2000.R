
# read in the raw data
datadir <- 'Montoya_2000'
filename <- 'Montoya_2000.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Parasitoids','Hosts','Parasitized.Total.mean','Parasitized.Total.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
