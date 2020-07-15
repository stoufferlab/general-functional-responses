
# read in the simplify the raw data
datadir <- 'Elliot_2005'
filename <- 'Elliot_2005_i2.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c('Preds','Prey','TotalEaten','Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
