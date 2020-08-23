
# read in the raw data
datadir <- 'Hassan_1976'
filename <- 'Hassan_1976_Ag.csv'
columns <- rbind(
	c('Npredator',      'Parasites'),
	c('Nprey',          'Hosts'),
	c('Nconsumed.mean', 'Parasiized.Total.mean'),
	c('Nconsumed.se',   'Parasitized.Total.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
