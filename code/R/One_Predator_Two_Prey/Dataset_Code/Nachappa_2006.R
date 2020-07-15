
# read in the simplify the raw data
datadir <- 'Nachappa_2006'
filename <- 'Nachappa_2006.csv'
d <- read.data(datadir, filename, "One_Predator_Two_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c('Pred','Prey1.TLS','Prey1.Eaten','Prey2.FAW','Prey2.Eaten','Time')]
	colnames(d) <- c("Npredator", "Nprey1", "Nconsumed1", "Nprey2", "Nconsumed2",'Time')
}
