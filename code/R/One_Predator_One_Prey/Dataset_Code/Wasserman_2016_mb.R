
# read in the simplify the raw data
datadir <- 'Wasserman_2016'
filename <- 'Wasserman_2016_mb.csv'
columns <- rbind(
	c('Npredator', 'Pred'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Prey.Eaten'),
	c('Time',      'Time')
)
