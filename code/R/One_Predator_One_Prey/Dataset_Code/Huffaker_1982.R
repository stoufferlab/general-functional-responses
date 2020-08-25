
# read in the raw data
datadir <- 'Huffaker_1982'
filename <- 'Huffaker_1982.csv'
columns <- rbind(
	c('Npredator',      'Preds'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'Attacks.mean'),
	c('Nconsumed.se',   'Attack.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
