
# read in the simplify the raw data
datadir <- 'Schenk_2005'
filename <- 'Schenk_2005.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Preds", "Prey", "PreyEaten")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")