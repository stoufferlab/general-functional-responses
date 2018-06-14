
# read in the simplify the raw data
datadir <- 'Chant_1966'
filename <- 'Chant_1966.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey','Eaten')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
