
# read in the raw data
datadir <- 'Crowley_1989'
filename <- 'Crowley_1989.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Nconsumed.Total.Mean'),
	c('Nconsumed.se',   'Nconsumed.Total.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
