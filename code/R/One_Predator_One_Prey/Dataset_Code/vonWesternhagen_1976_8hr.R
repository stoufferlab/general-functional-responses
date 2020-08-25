
# read in the data
datadir <- 'vonWesternhagen_1976'
filename <- 'vonWesternhagen_1976_8hr.csv'
columns <- rbind(
	c('Npredator',      'Predators'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Total.Eaten.mean'),
	c('Nconsumed.se',   'SE'),
	c('n',              'n'),
	c('Time',           'Hours')
)
