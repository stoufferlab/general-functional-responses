
# read in the data
datadir <- 'Nachappa_2006'
filename <- 'Nachappa_2006.csv'
columns <- rbind(
	c('Npredator',  'Pred'),
	c('Nprey1',     'Prey1.TLS'),
	c('Nconsumed1', 'Prey1.Eaten'),
	c('Nprey2',     'Prey2.FAW'),
	c('Nconsumed2', 'Prey2.Eaten'),
	c('Time',       'Time')
)
