
# read in the raw data
datadir <- 'Mills_2004'
filename <- 'Mills_2004.csv'
columns <- rbind(
	c('Npredator', 'Parasitoids'),
	c('Nprey',     'Hosts'),
	c('Nconsumed', 'TotalParasitized'),
	c('Time',      'Time')
)
