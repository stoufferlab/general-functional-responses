
# read in the simplify the raw data
datadir <- 'Medoc_2015'
filename <- 'Medoc_2015_be.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Total.Eaten'),
	c('Time',      'Time')
)
