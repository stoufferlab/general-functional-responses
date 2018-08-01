
# read in the raw data
datadir <- 'Kumar_1985'
filename <- 'Kumar_1985_Dl.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Parasitoids','Hosts','Parasitized.Mean','Parasitized.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")