
# read in the raw data
datadir <- 'Krylov_1992'
filename <- 'Krylov_1992i.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','TotalPreyEatenPer24hrs','TotalPreyEatenPer24hrs.SE','n','Time')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')