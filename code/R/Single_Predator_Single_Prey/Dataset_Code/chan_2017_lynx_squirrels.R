
# read in the simplify the raw data
datadir <- 'Chan_2017'
filename <- 'Chan2017_Lynx_annual.csv'
rawdata <- read.csv(paste(dropboxdir,datadir,filename,sep="/"))

# grab some info from the google doc
this.study <- study.info(datadir)
expttype <- this.study$expttype
Pminus1 <- this.study$Pminus1

# rename to standard column names used in fitting code
d <- rawdata[,c("Lynx.p.100km", "Squirrels.p.100km", "Squirrel.kills", "Lynx.Time.km")]
colnames(d) <- c("Npredator", "Nprey", "Nconsumed", "Time")

# DEBUG initial estimate of attack rate in Holling Type I
# NOTE: optimization is on log-transformed values
x0.hl <- c(-3)

# DEBUG initial estimate of attack rate in Hassell-Varley
# NOTE: optimization is on log-transformed values
x0.rd <- c(-3,0)

