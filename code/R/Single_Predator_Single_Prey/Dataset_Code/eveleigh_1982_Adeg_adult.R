
# read in the raw data
datadir <- 'Eveleigh_1982'
filename <- 'Eveleigh_1982_Adeg_adult.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# turn into a standard dataframe with standard column names
d <- rawdata[!is.na(rawdata$n),]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")

# DEBUG
# predation rates here are per predator
# must scale them and SEs to be absolute measures
d$Nconsumed.mean <- d$Nconsumed.mean * d$Npredator
d$Nconsumed.se <- d$Nconsumed.se * d$Npredator

# WARNING: this should probably be bootstrapped properly!
# bootstrap the experiment
d <- bootstrap.data(d, expttype)

# DEBUG
# grab units from paper
d$Time <- 1

# DEBUG initial estimate of attack rate in Holling Type I
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# DEBUG initial estimate of attack rate in Hassell-Varley
# NOTE: optimization is on log-transformed values
x0.rd <- c(-3,0)
