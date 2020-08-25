
# read in the raw data
datadir <- 'Vahl_2005'
filename <- 'Vahl_2005_Turn.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Eaten'),
	c('Time',      'Time')
)
