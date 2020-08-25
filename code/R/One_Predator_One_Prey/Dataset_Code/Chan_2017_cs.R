
# read in the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_cs.csv'
columns <- rbind(
	c('Npredator', 'Coyotes.p.100km'),
	c('Nprey',     'Squirrels.p.100km'),
	c('Nconsumed', 'Squirrel.kills'),
	c('Time',      'Coyote.Time.km')
)
