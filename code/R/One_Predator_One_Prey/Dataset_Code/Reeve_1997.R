
# read in the simplify the raw data
datadir <- 'Reeve_1997'
filename <- 'Reeve_1997.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Prey", "Pred", "Eaten",'Time')]
	colnames(d) <- c("Nprey", "Npredator", "Nconsumed", 'Time')
}
