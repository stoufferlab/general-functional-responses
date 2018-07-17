
# read in the simplify the raw data
datadir <- 'Wasserman_2016'
filename <- 'Wasserman_2016_mb.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# give columns the common name format
d <- rawdata[,c("Pred","Prey","Prey.Eaten")]
colnames(d) <- c("Npredator","Nprey","Nconsumed")
