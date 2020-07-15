
# read in the simplify the raw data
datadir <- 'Long_2012'
filename <- 'Long_2012.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Predators", "Prey", "Eaten",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
