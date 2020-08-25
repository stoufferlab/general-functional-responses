
# read in the raw data
datadir <- 'Kumar_1985'
filename <- 'Kumar_1985_Dl.csv'
columns <- rbind(
	c('Npredator',      'Parasitoids'),
	c('Nprey',          'Hosts'),
	c('Nconsumed.mean', 'Parasitized.Mean'),
	c('Nconsumed.se',   'Parasitized.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
