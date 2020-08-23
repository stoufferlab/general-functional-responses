
# read in the raw data
datadir <- 'Mansour_1991'
filename <- 'Mansour_1991.csv'
columns <- rbind(
	c('Npredator',      'Preds'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'TotalEaten.Mean'),
	c('Nconsumed.se',   'TotalEaten.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
