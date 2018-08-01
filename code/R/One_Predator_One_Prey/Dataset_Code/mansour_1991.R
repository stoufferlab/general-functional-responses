
# read in the raw data
datadir <- 'Mansour_1991'
filename <- 'Mansour_1991.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','TotalEaten.Mean','TotalEaten.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")