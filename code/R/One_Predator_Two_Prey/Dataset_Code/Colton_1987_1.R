
# read in the data
datadir <- 'Colton_1987'
filename <- 'Colton_1987_1.csv'
columns <- rbind(
	c('Npredator',  'Pred'),
	c('Nprey1',     'Prey1.SIN'),
	c('Nconsumed1', 'Prey1.Eaten'),
	c('Nprey2',     'Prey2.DIN'),
	c('Nconsumed2', 'Prey2.Eaten'),
	c('Time',       'Time')
)
