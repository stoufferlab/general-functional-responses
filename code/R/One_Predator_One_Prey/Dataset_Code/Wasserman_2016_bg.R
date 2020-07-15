
# read in the simplify the raw data
datadir <- 'Wasserman_2016'
filename <- 'Wasserman_2016_bg.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Pred", "Prey", "Prey.Eaten",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
