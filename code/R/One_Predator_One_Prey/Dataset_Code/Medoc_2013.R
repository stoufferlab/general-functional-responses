
# read in the simplify the raw data
datadir <- 'Medoc_2013'
filename <- 'Medoc_2013.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Total.Eaten'),
	c('Time',      'Time')
)
