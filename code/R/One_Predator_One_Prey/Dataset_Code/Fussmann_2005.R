
# read in the simplify the raw data
datadir <- 'Fussmann_2005'
filename <- 'Fussmann_2005.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey','Nconsumed','Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')

