
# read in the simplify the raw data
datadir <- 'Edwards_1961'
filename <- 'Edwards_1961_ts2.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Parasites','Hosts','Parasitized','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')