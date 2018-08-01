
# read in the simplify the raw data
datadir <- 'Mills_2004'
filename <- 'Mills_2004.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Parasitoids", "Hosts", "TotalParasitized")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")