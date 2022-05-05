
# read in the data
datadir <- 'Long_2012_2p'
filename <- 'Long_2012_2p.csv'
columns <- rbind(
	c('Npredator',  'Pred'),
	c('Nprey1',     'Prey1.Cr'),
	c('Nconsumed1', 'Prey1.Eaten'),
	c('Nprey2',     'Prey2.Cl'),
	c('Nconsumed2', 'Prey2.Eaten'),
	c('Time',       'Time')
)
