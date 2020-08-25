
# read in the raw data
datadir <- 'Eveleigh_1982'
filename <- 'Eveleigh_1982_ap.csv'
columns <- rbind(
	c('Npredator',      'Pred'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Nconsumed.mean'),
	c('Nconsumed.se',   'Nconsumed.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
