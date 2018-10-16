
# read in the simplify the raw data
datadir <- 'Mertz_1968'
filename <- 'Mertz_1968.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Preds", "Prey", "Eaten",'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
