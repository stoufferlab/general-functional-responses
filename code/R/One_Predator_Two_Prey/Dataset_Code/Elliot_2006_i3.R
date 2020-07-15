
# read in the simplify the raw data
datadir <- 'Elliot_2006'
filename <- 'Elliot_2006_Instar3.csv'
d <- read.data(datadir, filename, "One_Predator_Two_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c('Pred','Prey1.S','Prey1.Eaten.Mean','Prey1.Eaten.SE','Prey2.L','Prey2.Eaten.Mean','Prey2.Eaten.SE','n')]
	colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1.mean", "Nconsumed1.se", "Nprey2", "Nconsumed2.mean", "Nconsumed2.se", "n")
}
