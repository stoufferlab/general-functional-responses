
# read in the simplify the raw data
datadir <- 'Krylov_1992'
filename <- 'Krylov_1992_ii.csv'
columns <- rbind(
	c('Npredator',       'Preds'),
	c('Nprey1',          'Prey1.S'),
	c('Nconsumed1.mean', 'TotalPrey1Eaten.Mean'),
	c('Nconsumed1.se',   'TotalPrey1Eaten.SE'),
	c('Nprey2',          'Prey2.L'),
	c('Nconsumed2.mean', 'TotalPrey2Eaten.Mean'),
	c('Nconsumed2.se',   'TotalPrey2Eaten.SE'),
	c('n',               'n',),
	c('Time',            'Time')
)
