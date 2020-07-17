
# read in the simplify the raw data
datadir <- 'Lester_2002'
filename <- 'Lester_2002_Ty_d.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey1.TSSM','Prey1.Eaten','Prey2.ERM','Prey2.Eaten','Time')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')
