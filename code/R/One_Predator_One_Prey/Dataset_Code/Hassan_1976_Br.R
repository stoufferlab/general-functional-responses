
# read in the raw data
datadir <- 'Hassan_1976'
filename <- 'Hassan_1976_Br.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# necessary data for analysis
d <- rawdata[,c('Parasites','Hosts','Parasiized.Total.mean','Parasitized.Total.se','n','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')
