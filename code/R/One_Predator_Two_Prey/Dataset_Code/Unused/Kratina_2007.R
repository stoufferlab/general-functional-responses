
# read in the simplify the raw data
datadir <- 'Kratina_2007'
filename <- 'Kratina_2007.csv'
d <- read.data(datadir, filename, "One_Predator_Two_Prey", dropboxdir)

# rename to standard column names used in fitting code
if(!is.null(d)){
	d <- d[,c("Pred", "Prey1","Prey2", "Prey1Nconsumed.Total","Prey2Nconsumed.Total",'Time')]
	colnames(d) <- c('Npredator',"Nprey1","Nprey2","Nconsumed1","Nconsumed2",'Time')
}
