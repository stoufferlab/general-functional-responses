
# read in the simplify the raw data
datadir <- 'Fussmann_2005'
filename <- 'Fussmann_2005.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c('Preds','Prey','Nconsumed','Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
