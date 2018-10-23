
# read in the simplify the raw data
datadir <- 'Cronin_1993'
filename <- 'Cronin_1993.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Parasitoids.Mean','Host.Mean','TotalParasitized.Mean','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed","Time")
