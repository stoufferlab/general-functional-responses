
# read in the raw data
datadir <- 'Montoya_2000'
filename <- 'Montoya_2000.csv'
d <- read.data(datadir, filename, "One_Predator_One_Prey", dropboxdir)

# turn into a standard dataframe with standard column names
if(!is.null(d)){
	d <- d[,c('Parasitoids','Hosts','Parasitized.Total.mean','Parasitized.Total.SE','n','Time')]
	colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n",'Time')
}
