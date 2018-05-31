
# read in the simplify the raw data
datadir <- 'Lang_2012'
filename <- 'Lang_Interference_FR_Pter_10C.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# rename to standard column names used in fitting code
d <- rawdata[,c("PredNo", "Initial", "Killed")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")

# DEBUG DEBUG DEBUG 
# non-integer numbers consumed
d <- d[round(d$Nconsumed)==d$Nconsumed,]

# DEBUG
# grab units from paper
d$Time <- 1

# DEBUG initial estimate of attack rate in Holling Type I
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# DEBUG initial estimate of attack rate in Hassell-Varley
# NOTE: optimization is on log-transformed values
x0.rd <- c(-3,0)
