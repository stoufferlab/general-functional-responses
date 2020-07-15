
# read in the simplify the raw data
datadir <- 'Kratina_2009'
filename <- 'Kratina_2009.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("pred", "prey", "eaten",'Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed",'Time')
}
