
# read in the simplify the raw data
datadir <- 'Vahl_2005'
filename <- 'Vahl_2005_Knot.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Preds", "Prey", "Eaten",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
