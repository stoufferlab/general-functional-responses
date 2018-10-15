
# read in the simplify the raw data
datadir <- 'Buckel_2000'
filename <- 'Buckel_2000_small.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey1','Prey1.TotalEaten.Mean','Prey1.TotalEaten.SE','Prey2','Prey2.TotalEaten.Mean','Prey2.TotalEaten.SE','n')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
