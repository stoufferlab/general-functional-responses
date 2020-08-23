
# read in the raw data
datadir <- 'Kumar_1985'
filename <- 'Kumar_1985_Cc.csv'
columns <- rbind(
	c('Npredator',      'Parasitoids'),
	c('Nprey',          'Hosts'),
	c('Nconsumed.mean', 'Parasitized.Mean'),
	c('Nconsumed.se',   'Parasitized.Mean'),
	c('n',              'n'),
	c('Time',           'Time')
)
