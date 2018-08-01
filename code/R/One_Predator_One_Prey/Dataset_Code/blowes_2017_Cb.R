
# read in the simplify the raw data
datadir <- 'Blowes_2017'
filename <- 'Blowes_2017_Cb.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey','Bites')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
