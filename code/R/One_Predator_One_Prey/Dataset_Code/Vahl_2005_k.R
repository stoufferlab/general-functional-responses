
# read in the simplify the raw data
datadir <- 'Vahl_2005'
filename <- 'Vahl_2005_Knot.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Eaten'),
	c('Time',      'Time')
)
