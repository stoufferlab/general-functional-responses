
# read in the simplify the raw data
datadir <- 'Katz_1985'
filename <- 'Katz_1985.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Predators','Prey','Eaten.Total.mean','Eaten.Total.se','n','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')
