
# read in the simplify the raw data
datadir <- 'Elliot_2003'
filename <- 'Elliot_2003_Iso.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
# d <- rawdata[,c('Preds','Prey','PreyEatenTotal.Mean','PreyEatenTotal.SE','n')]
# colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")

# No SEs given, so treat means as raw data
d <- rawdata[,c('Preds','Prey','PreyEatenTotal.Mean','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')