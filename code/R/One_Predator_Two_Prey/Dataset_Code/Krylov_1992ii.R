
# read in the simplify the raw data
datadir <- 'Krylov_1992ii'
filename <- 'Krylov_1992ii.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey1.S','TotalPrey1Eaten.Mean','TotalPrey1Eaten.SE','Prey2.L','TotalPrey2Eaten.Mean','TotalPrey2Eaten.SE','n')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
