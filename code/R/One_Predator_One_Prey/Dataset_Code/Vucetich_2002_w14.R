
# read in the simplify the raw data
datadir <- 'Vucetich_2002'
filename <- 'IsleRoyale_Whole2014.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','Prey','TotalEatenPerObsTime','ObsTime')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
