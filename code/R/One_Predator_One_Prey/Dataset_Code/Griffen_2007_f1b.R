
# read in the raw data
datadir <- 'Griffen_2007'
filename <- 'Griffen_2007_f1b.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'TotalEaten.Mean'),
	c('Nconsumed.se',   'TotalEaten.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
