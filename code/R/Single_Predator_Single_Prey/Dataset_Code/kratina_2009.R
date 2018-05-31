
# read in the simplify the raw data
datadir <- 'Kratina_2009'
filename <- 'Kratina_2009_data.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# rename to standard column names used in fitting code
d <- rawdata[,c("pred", "prey", "eaten")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")

# DEBUG DEBUG DEBUG
d <- d[d$Nconsumed>0,]
d <- d[round(d$Nconsumed)==d$Nconsumed,]

# DEBUG
# grab units from paper
d$Time <- 4

# DEBUG initial estimate of attack rate in type-I FR
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# fit everything else
x0.rd <- c(-3,0)
