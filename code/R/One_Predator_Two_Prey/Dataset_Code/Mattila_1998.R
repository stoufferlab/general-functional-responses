
# read in the simplify the raw data
datadir <- 'Mattila_1998'
filename <- 'Mattila_1998.csv'
columns <- rbind(
	c('Npredator',       'Pred'),
	c('Nprey1',          'Prey1.B'),
	c('Nconsumed1.mean', 'Prey1Eaten.Mean'),
	c('Nconsumed1.se',   'Prey1Eaten.SE'),
	c('Nprey2',          'Prey2.M'),
	c('Nconsumed2.mean', 'Prey2Eaten.Mean'),
	c('Nconsumed2.se',   'Prey2Eaten.SE'),
	c('n',               'n',),
	c('Time',            'Time')
)
