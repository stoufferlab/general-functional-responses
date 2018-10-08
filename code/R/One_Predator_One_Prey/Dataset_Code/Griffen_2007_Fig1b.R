
# read in the simplify the raw data
datadir <- 'Griffen_2007'
filename <- 'Griffen_2007_Fig1b.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Pred','Prey','TotalEaten.Mean','TotalEaten.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")
