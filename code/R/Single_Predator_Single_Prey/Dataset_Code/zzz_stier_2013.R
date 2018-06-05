
# read in the simplify the raw data
datadir <- 'Stier_2013'
filename <- 'Stier_2013.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("pred", "total.settlers", "dead")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
