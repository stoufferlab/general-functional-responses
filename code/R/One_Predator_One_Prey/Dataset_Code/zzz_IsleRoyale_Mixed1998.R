
# read in the simplify the raw data
datadir <- 'Vucetich_2002'
filename <- 'IsleRoyale_Mixed1998.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Pred", "Prey", "TotalKillsPerDay")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")
