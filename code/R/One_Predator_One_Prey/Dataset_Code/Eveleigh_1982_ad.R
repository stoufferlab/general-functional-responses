
# read in the raw data
datadir <- 'Eveleigh_1982'
filename <- 'Eveleigh_1982_ad.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Total.killed.mean'),
	c('Nconsumed.se',   'Total.killed.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
