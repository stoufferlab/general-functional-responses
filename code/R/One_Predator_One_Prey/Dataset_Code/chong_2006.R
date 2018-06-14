
# read in the simplify the raw data
datadir <- 'Chong_2006'
filename <- 'Chong_2006_expmt2.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Parasitoids","Hosts","Parasitized")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
