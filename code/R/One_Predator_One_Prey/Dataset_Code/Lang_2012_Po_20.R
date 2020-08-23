
# read in the raw data
datadir <- 'Lang_2012'
filename <- 'Lang_2012_Poe_20C.csv'
columns <- rbind(
	c('Npredator', 'PredNo'),
	c('Nprey',     'Initial'),
	c('Nconsumed', 'Killed'),
	c('Time',      'Time')
)
