dat <- dsum

dat$log.Nconsumed.mean <- log(dat$Nconsumed.mean)
dat$log.Nconsumed.var <- log(dat$Nconsumed.var)
dat$log.Npredator <- log(dat$Npredator)

# remove treatments with zero mean or variance or NA
rem <- which(is.infinite(dat$log.Nconsumed.mean) | 
               is.infinite(dat$log.Nconsumed.var)|
               is.na(dat$log.Nconsumed.var))
if(length(rem) > 0){
  dat <- dat[-rem,]
  cat(paste("          ***Removed", length(rem), "of", nrow(dsum),
            "zero-mean or zero-variance treatment(s) from", datasetsName,"\n"))
}

# Does dataset have sufficient replication remaining? 
# (nMin and lMin specified in parent script)
cleared <- sum(dat$n >= nMin) >= lMin
# note that additional treatments with lower replication are kept too!
