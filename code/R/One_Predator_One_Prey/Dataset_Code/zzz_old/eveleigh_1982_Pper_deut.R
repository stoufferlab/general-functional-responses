
# read in the raw data
datadir <- 'Eveleigh_1982'
filename <- 'Eveleigh_1982_Pper_deut.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Pred','Prey','Total.killed.mean','Total.killed.se','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
