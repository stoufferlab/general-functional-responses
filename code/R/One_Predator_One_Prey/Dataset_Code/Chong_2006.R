
# read in the simplify the raw data
datadir <- 'Chong_2006'
filename <- 'Chong_2006.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Parasitoids","Hosts","Parasitized","Time")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
