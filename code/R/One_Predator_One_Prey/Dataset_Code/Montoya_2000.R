
# read in the raw data
datadir <- 'Montoya_2000'
filename <- 'Montoya_2000.csv'
columns <- rbind(
	c('Npredator',      'Parasitoids'),
	c('Nprey',          'Hosts'),
	c('Nconsumed.mean', 'Parasitized.Total.mean'),
	c('Nconsumed.se',   'Parasitized.Total.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
