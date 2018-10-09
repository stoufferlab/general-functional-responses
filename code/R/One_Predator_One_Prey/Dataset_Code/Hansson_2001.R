
# read in the simplify the raw data
datadir <- 'Hansson_2001'
filename <- 'Hansson_2001.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey','TotalPreyEaten','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')

