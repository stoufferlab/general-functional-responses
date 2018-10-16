
# read in the simplify the raw data
datadir <- 'Creswell_1998'
filename <- 'Creswell_1998.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey','Nconsumed', "Time")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
