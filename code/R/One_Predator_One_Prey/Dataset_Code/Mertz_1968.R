
# read in the simplify the raw data
datadir <- 'Mertz_1968'
filename <- 'Mertz_1968.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Eaten'),
	c('Time',      'Time')
)
