
# read in the raw data
datadir <- 'Kratina_2009'
filename <- 'Kratina_2009.csv'
columns <- rbind(
	c('Npredator', 'pred'),
	c('Nprey',     'prey'),
	c('Nconsumed', 'eaten'),
	c('Time',      'time')
)
