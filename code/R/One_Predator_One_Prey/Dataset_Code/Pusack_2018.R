
# read in the simplify the raw data
datadir <- 'Pusack_2018'
filename <- 'Pusack_2018.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("drill.abundance", "oyster.abundance", "total.no.oysters.consumed", "days")]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
}
