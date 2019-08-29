
# read in the simplify the raw data
datadir <- 'Elliot_2005'
filename <- 'Elliot_2005_i5.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey','TotalEaten','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
