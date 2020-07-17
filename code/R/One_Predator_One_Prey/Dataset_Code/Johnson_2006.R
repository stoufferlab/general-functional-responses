
# read in the simplify the raw data
datadir <- 'Johnson_2006'
filename <- 'Johnson_2006.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey','PreyEaten','Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
