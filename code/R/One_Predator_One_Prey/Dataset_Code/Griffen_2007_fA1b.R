
# read in the simplify the raw data
datadir <- 'Griffen_2007'
filename <- 'Griffen_2007_FigA1b.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','Total.Feeding.mean','Total.Feeding.se','n','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')
