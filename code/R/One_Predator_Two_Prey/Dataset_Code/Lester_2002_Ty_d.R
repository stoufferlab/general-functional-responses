
# read in the simplify the raw data
datadir <- 'Lester_2002'
filename <- 'Lester_2002_Ty_d.csv'
columns <- rbind(
	c('Npredator',  'Pred'),
	c('Nprey1',     'Prey1.TSSM'),
	c('Nconsumed1', 'Prey1.Eaten'),
	c('Nprey2',     'Prey2.ERM'),
	c('Nconsumed2', 'Prey2.Eaten'),
	c('Time',       'Time')
)
