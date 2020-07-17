
# read in the simplify the raw data
datadir <- 'Huffaker_1982'
filename <- 'Huffaker_1982.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','Attacks.mean','Attack.SE','n','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')
