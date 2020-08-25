
# read in the raw data
datadir <- 'Edwards_1961'
filename <- 'Edwards_1961_ts2.csv'
columns <- rbind(
	c('Npredator', 'Parasites'),
	c('Nprey',     'Hosts'),
	c('Nconsumed', 'Parasitized'),
	c('Time',      'Time')
)
