
# read in the simplify the raw data
datadir <- 'Mattila_1998'
filename <- 'Mattila_1998.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey1.B','Prey1Eaten.Mean','Prey1Eaten.SE','Prey2.M','Prey2Eaten.Mean','Prey2Eaten.SE','n')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
