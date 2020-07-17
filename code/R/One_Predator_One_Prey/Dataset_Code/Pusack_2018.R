
# read in the simplify the raw data
datadir <- 'Pusack_2018'
filename <- 'Pusack_2018.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("drill.abundance", "oyster.abundance", "total.no.oysters.consumed", "days")]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
