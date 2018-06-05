
# read in the raw data
datadir <- 'Kfir_1983'
filename <- 'Kfir_1983.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Pred','Prey','Parasitized.Total.mean','Parasitized.Total.se','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
