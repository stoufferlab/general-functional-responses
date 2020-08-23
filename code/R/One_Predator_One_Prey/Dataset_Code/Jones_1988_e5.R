
# read in the raw data
datadir <- 'Jones_1988'
filename <- 'Jones_1988_e5.csv'
columns <- rbind(
	c('Npredator', 'Parasitoids'),
	c('Nprey',     'Hosts'),
	c('Nconsumed', 'Attacked'),
	c('Time',      'Time')
)
