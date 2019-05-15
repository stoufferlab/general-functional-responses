
# read in the simplify the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_ls.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Lynx.p.100km", "Squirrels.p.100km", "Squirrel.kills", "Lynx.Time.km")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
