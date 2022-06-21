dat <- dsum

dat$log.Nconsumed1.mean <- log(dat$Nconsumed1.mean)
dat$log.Nconsumed2.mean <- log(dat$Nconsumed2.mean)
dat$log.Nconsumed1.var <- log(dat$Nconsumed1.var)
dat$log.Nconsumed2.var <- log(dat$Nconsumed2.var)
dat$log.Nprey1 <- log(dat$Nprey1)
dat$log.Nprey2 <- log(dat$Nprey2)

# create prey-specific datasets to enable removal of 
# treatments with zero mean or variance or NAs
dat1 <- dat2 <- dat

rem1 <- which(is.infinite(dat1$log.Nconsumed1.mean) | 
                is.infinite(dat1$log.Nconsumed1.var)|
                is.na(dat1$log.Nconsumed1.var))
rem2 <- which(is.infinite(dat2$log.Nconsumed2.mean) | 
                is.infinite(dat2$log.Nconsumed2.var)|
                is.na(dat2$log.Nconsumed2.var))

if(length(rem1) > 0){ dat1 <- dat1[-rem1,] }
if(length(rem2) > 0){ dat2 <- dat2[-rem2,] }
if(length(rem1) > 0 | length(rem2) > 0){
  cat(paste("          ***Removed", length(rem1), "of", nrow(dsum), 
            "and",
            length(rem2), "of", nrow(dsum),
            "zero-mean or zero-variance treatment(s) from", datasetsName,
            "\nfor prey 1 and prey 2, respectively.\n"))
}


# Make log(NPrey=0) = -Inf -> 0 to enable models for the effect of prey1 on prey2's 
# mean-variance relationship (and vice versa)
mod1 <- which(is.infinite(dat1$log.Nprey2))
mod2 <- which(is.infinite(dat2$log.Nprey1))

if(length(mod1) > 0){ dat1$log.Nprey2[mod1] <- 0 }
if(length(mod2) > 0){ dat2$log.Nprey1[mod2] <- 0 }
if(length(mod1) > 0 | length(mod2) > 0){
  cat(paste("          ***Replaced", length(mod1), "of", nrow(dsum), 
            "and",
            length(mod2), "of", nrow(dsum),
            "log(Nprey) = -inf values with zeros in", datasetsName,
            "\nfor prey 1 and prey 2, respectively", 
            "to enable fitting of alternative prey effects.\n"))
}

# Do _both_ datasets have sufficient replication remaining? 
# (nMin and lMin are both specified in parent script)
# Note: additional treatments with lower replication are kept too!
cleared <- sum(dat1$n >= nMin) >= lMin & sum(dat2$n >= nMin) >= lMin

# For aggregating datasets, relabel columns
colnames(dat1) <- sub('1', '.MainPrey', colnames(dat1))
colnames(dat1) <- sub('2', '.AltPrey', colnames(dat1))
colnames(dat2) <- sub('2', '.MainPrey', colnames(dat2))
colnames(dat2) <- sub('1', '.AltPrey', colnames(dat2))

dat1 <- data.frame(MainPrey = 'Prey1', dat1)
dat2 <- data.frame(MainPrey = 'Prey2', dat2)

dat <- rbind(dat1, dat2)


