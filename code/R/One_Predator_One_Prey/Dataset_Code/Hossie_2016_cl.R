
# read in the raw data
datadir <- 'Hossie_2016'
filename <- 'Hossie_2016_clumped.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Predators", "Prey", "Killed",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
