
# read in the simplify the raw data
datadir <- 'Kalinkat_2011'
filename <- 'Kalinkat_2011_Pard.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c('Preds','large_prey_N0','large_prey_NE','small_prey_N0','small_prey_NE','Time')]
	colnames(d) <- c('Npredator',"Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')
