
# read in the simplify the raw data
datadir <- 'Ranta_1985'
filename <- 'Ranta_1985_18.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred',
                'Prey1.L','Prey1.Eaten.Mean','Prey1.Eaten.SE',
                'Prey2.M','Prey2.Eaten.Mean','Prey2.Eaten.SE','n')]
colnames(d) <- c("Npredator", 
                 "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", 
                 "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
