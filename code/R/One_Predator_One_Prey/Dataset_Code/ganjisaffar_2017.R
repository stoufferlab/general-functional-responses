
# read in the simplify the raw data
datadir <- 'Ganjisaffar_2017'
filename <- 'Ganjisaffar_2017.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Predator.Density','Prey.Density','Eaten')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")