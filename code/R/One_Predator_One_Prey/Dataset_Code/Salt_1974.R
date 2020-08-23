
# read in the simplify the raw data
datadir <- 'Salt_1974'
filename <- 'Salt_1974.csv'
columns <- rbind(
	c('Npredator',      'Preds'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'FeedingRate.Total'),
	c('Nconsumed.se',   'FeedingRate.Total.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
