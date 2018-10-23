
# read in the raw data
datadir <- 'Jones_1988'
filename <- 'Jones_1986_Exp5.1.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Parasitoids','Hosts','Attacked','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
