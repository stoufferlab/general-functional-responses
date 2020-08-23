
# read in the simplify the raw data
datadir <- 'Vucetich_2002'
filename <- 'IsleRoyale_Whole1998.csv'
columns <- rbind(
	c('Npredator', 'Preds'),
	c('Nprey',     'Prey'),
	c('Nconsumed', 'TotalEatenPerObsTime'),
	c('Time',      'ObsTime')
)
