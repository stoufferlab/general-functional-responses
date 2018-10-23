
# read in the simplify the raw data
datadir <- 'Elliot_2003'
filename <- 'Elliot_2003_Din.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','PreyEatenTotal.Mean','PreyEatenTotal.SE','n', 'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", 'n', 'Time')
