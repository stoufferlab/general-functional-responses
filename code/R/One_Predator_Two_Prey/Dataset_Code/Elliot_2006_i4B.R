
# read in the simplify the raw data
datadir <- 'Elliot_2006'
filename <- 'Elliot_2006_Instar4Baet.csv'
columns <- rbind(
	c('Npredator',       'Pred'),
	c('Nprey1',          'Prey1.B'),
	c('Nconsumed1.mean', 'Prey1.Eaten.Mean'),
	c('Nconsumed1.se',   'Prey1.Eaten.SE'),
	c('Nprey2',          'Prey2.L'),
	c('Nconsumed2.mean', 'Prey2.Eaten.Mean'),
	c('Nconsumed2.se',   'Prey2.Eaten.SE'),
	c('n',               'n',),
	c('Time',            'Time')
)
