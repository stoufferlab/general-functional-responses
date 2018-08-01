
# read in the raw data
datadir <- 'Uttley_1980'
filename <- 'Uttley_1980_NotoIns1.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Pred','Prey','PreyEatenTotal.Mean','PreyEatenTotal.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")