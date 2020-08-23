
# read in the simplify the raw data
datadir <- 'Griffen_2007'
filename <- 'Griffen_2007_fA1b.csv'
columns <- rbind(
	c('Npredator',      'Preds'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Total.Feeding.mean'),
	c('Nconsumed.se',   'Total.Feeding.se'),
	c('n',              'n'),
	c('Time',           'Time')
)
