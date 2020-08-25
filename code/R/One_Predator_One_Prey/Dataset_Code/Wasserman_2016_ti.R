
# read in the raw data
datadir <- 'Wasserman_2016'
filename <- 'Wasserman_2016_ti.csv'
columns <- rbind(
	c('Npredator', 'Pred'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Prey.Eaten'),
	c('Time',      'Time')
)
