
# read in the simplify the raw data
datadir <- 'Nachappa_2006'
filename <- 'Nachappa_2006.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey1.TLS','Prey1.Eaten','Prey2.FAW','Prey2.Eaten','Time')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')
