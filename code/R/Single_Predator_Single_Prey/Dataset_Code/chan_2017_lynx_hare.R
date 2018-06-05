
# read in the simplify the raw data
datadir <- 'Chan_2017'
filename <- 'Chan2017_Lynx_annual.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("Lynx.p.100km", "Hares.p.100km", "HareKills", "Lynx.Time.km")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
