
# read in the simplify the raw data
datadir <- 'Reeve_1997'
filename <- 'Reeve_1997.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Prey", "Pred", "Eaten")]
colnames(d) <- c("Nprey", "Npredator", "Nconsumed")

# WARNING non-integers and NA values muck things up
d <- d[!is.na(d$Nconsumed),]
d$Nconsumed <- round(d$Nconsumed)
