
# read in the raw data
datadir <- 'Hossie_2016'
filename <- 'Hossie_2016_clumped.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Predators", "Prey", "Killed")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
