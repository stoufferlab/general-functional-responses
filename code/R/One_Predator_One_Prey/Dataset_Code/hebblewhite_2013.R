
# read in the simplify the raw data
datadir <- 'Hebblewhite_2013'
filename <- 'Hebblewhite_2013.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey','KillRate')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")

