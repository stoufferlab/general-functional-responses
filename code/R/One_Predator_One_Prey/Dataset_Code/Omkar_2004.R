
# read in the simplify the raw data
datadir <- 'Omkar_2004'
filename <- 'Omkar_2004.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Preds", "Prey", "Eaten", "Time.hrs")]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
}
