
# read in the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_coyote_squirrel.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Coyotes.p.100km", "Squirrels.p.100km", "Squirrel.kills", "Coyote.Time.km")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
