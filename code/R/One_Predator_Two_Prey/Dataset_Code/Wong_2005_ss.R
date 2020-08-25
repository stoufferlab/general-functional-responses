
# read in the simplify the raw data
datadir <- 'Wong_2005'
filename <- 'Wong_2005_ss.csv'
columns <- rbind(
	c('Npredator',  'Preds'),
	c('Nprey1',     'Prey1.m'),
	c('Nconsumed1', 'Prey1Eaten'),
	c('Nprey2',     'Prey2.s'),
	c('Nconsumed2', 'Prey2Eaten'),
	c('Time',       'Time')
)
