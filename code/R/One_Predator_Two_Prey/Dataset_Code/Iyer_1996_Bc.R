
# read in the simplify the raw data
datadir <- 'Iyer_1996'
filename <- 'Iyer_1996_Bc.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey1.Bc','Prey1Eaten.Mean','Prey1Eaten.SE','Prey2.Hm','Prey2Eaten.Mean','Prey2Eaten.SE','n','Time')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n", "Time")
