
# read in the raw data
datadir <- 'Elliot_2005'
filename <- 'Elliot_2005_i5.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'TotalEaten'),
	c('Time',      'Time')
)
