
# read in the raw data
datadir <- 'Johnson_2006'
filename <- 'Johnson_2006.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'PreyEaten'),
	c('Time',      'Time')
)
