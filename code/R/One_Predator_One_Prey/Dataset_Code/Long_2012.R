
# read in the simplify the raw data
datadir <- 'Long_2012'
filename <- 'Long_2012a.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Predators", "Prey", "Eaten",'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
