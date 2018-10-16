
# read in the simplify the raw data
datadir <- 'Akre_1979'
filename <- 'Akre_1979_Satiated.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey1.Sim','Prey1.Eaten.Mean','Prey1.Eaten.SE','Prey2.Daph','Prey2.Eaten.Mean','Prey2.Eaten.SE','n')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
