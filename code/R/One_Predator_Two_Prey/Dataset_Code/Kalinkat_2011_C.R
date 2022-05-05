
# read in the data
datadir <- 'Kalinkat_2011'
filename <- 'Kalinkat_2011_Cal.csv'
columns <- rbind(
	c('Npredator',  'Preds'),
	c('Nprey1',     'large_prey_N0'),
	c('Nconsumed1', 'large_prey_NE'),
	c('Nprey2',     'small_prey_N0'),
	c('Nconsumed2', 'small_prey_NE'),
	c('Time',       'Time')
)
