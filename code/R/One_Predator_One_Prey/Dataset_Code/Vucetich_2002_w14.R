
# read in the simplify the raw data
datadir <- 'Vucetich_2002'
filename <- 'IsleRoyale_Whole2014.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c('Preds','Prey','TotalEatenPerObsTime','ObsTime')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")
}
