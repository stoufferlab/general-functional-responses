
# read in the simplify the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_lh.csv'
columns <- rbind(
	c('Npredator', 'Lynx.p.100km'),
	c('Nprey',     'Hares.p.100km'),
	c('Nconsumed', 'HareKills'),
	c('Time',      'Lynx.Time.km')
)
