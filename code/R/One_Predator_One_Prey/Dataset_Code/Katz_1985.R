
# read in the raw data
datadir <- 'Katz_1985'
filename <- 'Katz_1985.csv'
columns <- rbind(
	c('Npredator',      'Predators'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Eaten.Total.mean'),
	c('Nconsumed.se',   'Eaten.Total.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
