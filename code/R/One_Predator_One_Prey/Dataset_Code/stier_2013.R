
# read in the simplify the raw data
datadir <- 'Stier_2013'
filename <- 'Stier_2013.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# rename to standard column names used in fitting code
d <- rawdata[,c("predator_density_per_patch", "initial_prey_density_per_patch", "prey_eaten")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")