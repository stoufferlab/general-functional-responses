
# read in the raw data
datadir <- 'Kfir_1983'
filename <- 'Kfir_1983.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Parasitized.Total.mean'),
	c('Nconsumed.se',   'Parasitized.Total.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
