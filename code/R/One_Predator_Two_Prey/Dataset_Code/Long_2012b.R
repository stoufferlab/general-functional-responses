
# read in the simplify the raw data
datadir <- 'Long_2012'
filename <- 'Long_2012b.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey1.Cr','Prey1.Eaten','Prey2.Cl','Prey2.Eaten','Time')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')
