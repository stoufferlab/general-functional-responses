
# read in the raw data
datadir <- 'Krylov_1992'
filename <- 'Krylov_1992_i.csv'
columns <- rbind(
	c('Npredator',      'Preds'),
	c('Nprey',          'Prey'),
	c('Nconsumed.mean', 'TotalPreyEatenPer24hrs'),
	c('Nconsumed.se',   'TotalPreyEatenPer24hrs.SE'),
	c('n',              'n'),
	c('Time',           'Time')
)
