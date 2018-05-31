
# read in the simplify the raw data
datadir <- 'Chant_1966'
filename <- 'Chant_1966.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# rename to standard column names used in fitting code
d <- rawdata
colnames(d) <- c("Npredator", "Nprey", "Nconsumed")

# DEBUG
# grab units from paper
d$Time <- 1

# DEBUG initial estimate of attack rate in type-I FR
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# fit everything else
x0.rd <- c(-3,0)
