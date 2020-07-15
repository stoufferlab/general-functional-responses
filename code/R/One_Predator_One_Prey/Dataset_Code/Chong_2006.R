
# read in the simplify the raw data
datadir <- 'Chong_2006'
filename <- 'Chong_2006.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Parasitoids","Hosts","Parasitized","Time")]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
}
