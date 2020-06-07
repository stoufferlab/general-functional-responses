
# read in the raw data
datadir <- 'Prokopenko_2017'
filename <- 'Prokopenko_2017.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# Prokopenko et al. fit continuous-time models to data recorded at ~3 min intervals over the course of ~1hr.  They did not replace prey.  
# Since time-points are not independent, this script aggregates the interval observations to count the total number of prey eaten over the course of the hour.

# label unique replicates
rawdata$rep <- paste(rawdata$well.id, rawdata$date, sep='-')

# add up killed prey over time period for each replicate
Nconsumed <- aggregate(list(Nconsumed=rawdata$Killed.well), 
                 by=list(rep=rawdata$rep), sum)

# add up average intervals to get total (average) time for each replicate
# (averages between min and max possible interval periods were used in Prokopenko)
Time <- aggregate(list(Time=rawdata$avgtime), 
                    by=list(rep=rawdata$rep), sum)

# identify predator count of each replicate
Npredator <- aggregate(list(Npredator=rawdata$P.well), 
                  by=list(rep=rawdata$rep), mean)

# identify starting prey count of each replicate
Nprey <- aggregate(list(Nprey=rawdata$N.well), 
                       by=list(rep=rawdata$rep), max)

dat <- merge(merge(merge(Npredator, Nprey), Nconsumed), Time)

# drop replicate identifier
d <- dat[,-1]

