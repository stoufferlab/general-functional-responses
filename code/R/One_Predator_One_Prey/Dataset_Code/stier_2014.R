
# read in the simplify the raw data
datadir <- 'Stier_2014'
filename <- 'Stier_2014.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("pred", "total.settlers", "dead",'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
