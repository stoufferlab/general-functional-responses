
# read in the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_ls.csv'
columns <- rbind(
	c('Npredator', 'Lynx.p.100km'),
	c('Nprey',     'Squirrels.p.100km'),
	c('Nconsumed', 'Squirrel.kills'),
	c('Time',      'Lynx.Time.km')
)
