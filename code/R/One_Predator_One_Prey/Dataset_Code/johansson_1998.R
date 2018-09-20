
# read in the raw data
datadir <- 'Johansson_1998'
filename <- 'Johansson_1998.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','Prey.Eaten.Total.Mean','Prey.Eaten.Total.se','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")