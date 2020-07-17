
# read in the simplify the raw data
datadir <- 'Omkar_2004'
filename <- 'Omkar_2004.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Preds", "Prey", "Eaten", "Time.hrs")]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
