
# read in the simplify the raw data
datadir <- 'Lang_2012'
filename <- 'Lang_2012_Pter_10C.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("PredNo", "Initial", "Killed",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
