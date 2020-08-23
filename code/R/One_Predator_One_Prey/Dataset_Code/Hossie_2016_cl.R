
# read in the raw data
datadir <- 'Hossie_2016'
filename <- 'Hossie_2016_clumped.csv'
columns <- rbind(
	c('Npredator', 'Predators'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'Killed'),
	c('Time',      'Time')
)
