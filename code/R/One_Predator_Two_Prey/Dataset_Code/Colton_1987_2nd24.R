
# read in the simplify the raw data
datadir <- 'Colton_1987'
filename <- 'Colton_1987_2nd24.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Pred','Prey1.SIN','Prey1.Eaten','Prey2.DIN','Prey2.Eaten')]
colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2")
