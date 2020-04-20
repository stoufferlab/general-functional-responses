
# read in the simplify the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_c.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Coyotes.p.100km','Hares.p.100km','HareKills','Squirrels.p.100km','Squirrel.kills','Coyote.Time.km')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')
