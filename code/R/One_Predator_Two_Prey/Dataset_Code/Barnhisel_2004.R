
# read in the simplify the raw data
datadir <- 'Barnhisel_2004'
filename <- 'Barnhisel_2004.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey1.Daphnia','Prey1.Eaten','Prey1.Eaten.SE','Prey2.Bythotrephes','Prey2.Eaten','Prey2.Eaten.se','n')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
