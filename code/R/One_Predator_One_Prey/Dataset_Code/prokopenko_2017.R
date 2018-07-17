
# read in the simplify the raw data
datadir <- 'Prokopenko_2017'
filename <- 'Prokopenko_2017.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("P.well", "N.well", "Killed.well", "avgtime")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
