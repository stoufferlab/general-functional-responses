
# read in the raw data
datadir <- 'Jones_1988'
filename <- 'Jones_1988_Fig1.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Parasitoids','Hosts','Parasitized.Mean','Parasitized.se','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
