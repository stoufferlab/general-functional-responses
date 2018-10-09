
# read in the simplify the raw data
datadir <- 'Kratina_2009'
filename <- 'Kratina_2009.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("pred", "prey", "eaten",'Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')

# WARNING
d <- d[d$Nconsumed>0,]
d <- d[round(d$Nconsumed)==d$Nconsumed,]
