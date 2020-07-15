
# read in the simplify the raw data
datadir <- 'Johnson_2006'
filename <- 'Johnson_2006.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c('Preds','Prey','PreyEaten','Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
