
# read in the simplify the raw data
datadir <- 'Edwards_1961'
filename <- 'Edwards_1961_ts1.csv'
columns <- rbind(
	c('Npredator', 'Parasite'),
	c('Nprey',     'Host'),
	c('Nconsumed', 'Parasitized'),
	c('Time',      'Time')
)
