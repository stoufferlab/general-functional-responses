
# read in the raw data
datadir <- 'Pusack_2018'
filename <- 'Pusack_2018.csv'
columns <- rbind(
	c('Npredator', 'drill.abundance'),
	c('Nprey',     'oyster.abundance'),
	c('Nconsumed', 'total.no.oysters.consumed'),
	c('Time',      'days')
)
