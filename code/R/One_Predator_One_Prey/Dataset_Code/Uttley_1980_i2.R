
# read in the raw data
datadir <- 'Uttley_1980'
filename <- 'Uttley_1980_i2.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'PreyEatenTotal.Mean'),
	c('Nconsumed.se',   'PreyEatenTotal.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
