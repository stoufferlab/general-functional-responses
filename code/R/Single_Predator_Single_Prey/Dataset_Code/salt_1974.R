
# read in the simplify the raw data
datadir <- 'Salt_1974'
filename <- 'Salt_1974.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# turn into a standard dataframe with standard column names
d <- rawdata[,c('Preds','Prey','FeedingRate.Total','FeedingRate.Total.SE','n')]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed.mean", "Nconsumed.se", "n")

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
