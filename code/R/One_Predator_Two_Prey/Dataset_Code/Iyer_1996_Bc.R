
# read in the simplify the raw data
datadir <- 'Iyer_1996'
filename <- 'Iyer_1996_Bc.csv'
columns <- rbind(
	c('Npredator',       'Preds'),
	c('Nprey1',          'Prey1.Bc'),
	c('Nconsumed1.mean', 'Prey1Eaten.Mean'),
	c('Nconsumed1.se',   'Prey1Eaten.SE'),
	c('Nprey2',          'Prey2.Hm'),
	c('Nconsumed2.mean', 'Prey2Eaten.Mean'),
	c('Nconsumed2.se',   'Prey2Eaten.SE'),
	c('n',               'n',),
	c('Time',            'Time')
)
