
# read in the simplify the raw data
datadir <- 'Medoc_2015'
filename <- 'Medoc_2015_be.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Preds", "Prey", "Total.Eaten",'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')