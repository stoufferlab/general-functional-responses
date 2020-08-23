
# read in the simplify the raw data
datadir <- 'Walde_1984'
filename <- 'Walde_1984.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Captures.Total.mean'),
	c('Nconsumed.se',   'Captures.Total.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
