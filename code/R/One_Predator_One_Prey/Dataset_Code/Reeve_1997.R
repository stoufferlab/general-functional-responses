
# read in the raw data
datadir <- 'Reeve_1997'
filename <- 'Reeve_1997.csv'
columns <- rbind(
	c('Npredator', 'Pred'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Eaten'),
	c('Time',      'Time')
)
