
# read in the raw data
datadir <- 'Griffen_2007'
filename <- 'Griffen_2007_f1a.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'TotalEaten.Mean'),
	c('Nconsumed.se',   'TotalEaten.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
