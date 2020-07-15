
# read in the simplify the raw data
datadir <- 'Mills_2004'
filename <- 'Mills_2004.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Parasitoids", "Hosts", "TotalParasitized",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
