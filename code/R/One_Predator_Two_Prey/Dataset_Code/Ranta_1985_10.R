
# read in the data
datadir <- 'Ranta_1985'
filename <- 'Ranta_1985_10.csv'
columns <- rbind(
	c('Npredator',       'Pred'),
	c('Nprey1',          'Prey1.L'),
	c('Nconsumed1.mean', 'Prey1.Eaten.Mean'),
	c('Nconsumed1.se',   'Prey1.Eaten.SE'),
	c('Nprey2',          'Prey2.M'),
	c('Nconsumed2.mean', 'Prey2.Eaten.Mean'),
	c('Nconsumed2.se',   'Prey2.Eaten.SE'),
	c('n',               'n'),
	c('Time',            'Time')
)
