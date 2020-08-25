
# read in the raw data
datadir <- 'Chan_2017'
filename <- 'Chan_2017_ch.csv'
columns <- rbind(
	c('Npredator', 'Coyotes.p.100km'),
	c('Nprey',     'Hares.p.100km'),
	c('Nconsumed', 'HareKills'),
	c('Time',      'Coyote.Time.km')
)
