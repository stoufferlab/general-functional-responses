# variables to include
# max:min N
# max:min P
# max P:N
# min P:N
# var P:N
# levels P
# levels N
# parasitoid/pred
# replacement
library(bbmle)

# load in the actual fits (or rather a list of fits and other stuff)
load("functional_response_fits.Rdata")

# a container for all the different attributes of the data being fit
rd <- matrix(NA, length(ffr.fits), 5)
rownames(rd) <- names(ffr.fits)
colnames(rd) <- names(ffr.fits[[1]]$study.info)

# grab some info from the study.info beast
for(ff in 1:length(ffr.fits)){
	rd[ff,] <- unlist(ffr.fits[[ff]]$study.info)
}
rd <- as.data.frame(rd)
rd$expttype <- as.factor(rd$expttype)
for(cc in c('Pminus1','delong','predator')){
	rd[,cc] <- as.logical(rd[,cc])
}

#########################################
# start adding further columns manually
#########################################

# ratio dependent m
rd$exponent <- unlist(lapply(ffr.fits, function(x) (coef(x[["Arditi-Akcakaya"]])["exponent"])))
rd$exponent.se <- unlist(lapply(ffr.fits, function(x) coef(summary(x[["Arditi-Akcakaya"]]))["exponent","Std. Error"]))

# ratio dependent ah
rd$ah <- unlist(lapply(ffr.fits, function(x){
	exp(coef(x[["Arditi-Akcakaya"]])["attack"]) * exp(coef(x[["Arditi-Akcakaya"]])["handling"])
}))

# beddington-deangelis c
rd$bda.interference <- unlist(lapply(ffr.fits, function(x) (coef(x[["Beddington-DeAngelis"]])["interference"])))
rd$bda.interference.se <- unlist(lapply(ffr.fits, function(x) coef(summary(x[["Beddington-DeAngelis"]]))["interference","Std. Error"]))

# crowley-martin c
rd$cm.interference <- unlist(lapply(ffr.fits, function(x) (coef(x[["Crowley-Martin"]])["interference"])))
rd$cm.interference.se <- unlist(lapply(ffr.fits, function(x) coef(summary(x[["Crowley-Martin"]]))["interference","Std. Error"]))

# sample size
rd$sample.size <- unlist(lapply(ffr.fits, function(x) x$sample.size))

# ratio of max to min numbers of prey
rd$prey.ratio <- unlist(lapply(ffr.fits, function(x){max(x$data$Nprey) / min(x$data$Nprey)}))

# ratio of max to min numbers of predators
rd$pred.ratio <- unlist(lapply(ffr.fits, function(x){max(x$data$Npredator) / min(x$data$Npredator)}))

# max of prey:predators ratio
rd$max.prey.pred.ratio <- unlist(lapply(ffr.fits, function(x){max(x$data$Nprey / x$data$Npredator)}))

# min of prey:predators ratio
rd$min.prey.pred.ratio <- unlist(lapply(ffr.fits, function(x){min(x$data$Nprey / x$data$Npredator)}))

# variance of prey:predators ratio
rd$var.prey.pred.ratio <- unlist(lapply(ffr.fits, function(x){var(x$data$Nprey / x$data$Npredator)}))

# number of prey abundance treatments
rd$prey.levels <- unlist(lapply(ffr.fits, function(x){length(unique(x$data$Nprey))}))

# number of predator abundance treatments
rd$pred.levels <- unlist(lapply(ffr.fits, function(x){length(unique(x$data$Npredator))}))

# # if we wish to log some things some things
# to.log <- c("prey.ratio","pred.ratio","max.prey.pred.ratio","min.prey.pred.ratio","var.prey.pred.ratio")
# rd[,to.log] <- log(rd[,to.log])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figures
# ~~~~~~~
par(mfrow=c(2,2))

plot(exp(exponent) ~ sample.size, rd, xlab='sample size', ylab='m')
for(i in 1:nrow(rd)){
	n <- rd[i,"sample.size"]
	m <- rd[i,"exponent"]
	m.se <- rd[i,"exponent.se"]
	segments(n, exp(m-m.se), n, exp(m+m.se))
}

plot(exp(bda.interference) ~ sample.size, rd, xlab='sample size', ylab='beddington-deangelis', ylim=c(0,4))
for(i in 1:nrow(rd)){
	n <- rd[i,"sample.size"]
	c <- rd[i,"bda.interference"]
	c.se <- rd[i,"bda.interference.se"]
	segments(n, exp(c-c.se), n, exp(c+c.se))
}

plot(exp(cm.interference) ~ sample.size, rd, xlab='sample size', ylab='crowley-martin', ylim=c(0,4))
for(i in 1:nrow(rd)){
	n <- rd[i,"sample.size"]
	c <- rd[i,"cm.interference"]
	c.se <- rd[i,"cm.interference.se"]
	segments(n, exp(c-c.se), n, exp(c+c.se))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regressions
# ~~~~~~~~~~~
# Our hypothesis is that the source of sample-size dependence in the estimator for m ('m-hat') is due to (1) predator individuals varying in their interference effects, and (2) interference effects being less than one (i.e. a predator is not a full predator).  That is, letting 'g' stand for how much a given predator individual actually represents and g_r the average g in a given experimental replicate and g-bar the average g_r across replicates, the assumed ratio-dependent model P^m actually equals (g-bar * P)^m-hat, hence the estimator 'm-hat' will overestimate 'm' according to m-hat = m * log(P) / log(g-bar * P) = m * log(P) / (log(P) + log(g-bar)), which will result in m-hat > m whenever (1) the true mean interference effect 'g-bar' is less than 1 (i.e. log(g-bar) < 0) or (2) there is a distribution of g values across individuals (i.e. g has non-zero variance, i.e. indivdual variation) such that even when the true g-bar = 1, log(est-g-bar) will be less than 0 because the number of experimental replicates is low enough that convergence on a normal distribution of g-bar (in accord with the central limit theorem) has not been achieved (the geometric mean of g's is < arithmetic mean of g's).

# Given m-hat = m * log(P) / (log(P) + log(g)), we predict that 'm-hat' will coverge on 'm' with increasing predator abundance (regardless of the true g-bar) and with increasing replication (sample size n), and that the effect of increasing replication will diminish as predator density increases.  In other words, we expect there to be a statistical interaction between 'n' and 'P'.

# Practically speaking, since the lowest P used in experiments is almost always P=1, we should be able to use the max:min ratio of predator densities as a proxy for the maximum P used in the experiment.  Hence we use pred.ratio. (Perhaps the mean P is more relevant depending on how P increases across treatment levels -- 1, 2, 3, 4 vs. 1, 2, 4, 8 ?)

# In the data.frame, 'sample.size' refers to the total number of replicates.  However, in terms of the CTL, we expect that number of replicates per treatment level to be the relavant variable.  We haven't kept track of this information for each treatment level (and in some studies it does indeed vary), so for now just assume the average replication and that the experiments are all fully factorial (which they are not).

rd$sample.size.per.rep<-rd$sample.size/(rd$prey.levels*rd$pred.levels)

# According to the Berry-Esseen theorem, the (maximum) rate of convergence on the normal distribution is bounded by n^{-1/2} = 1/sqrt(n), where n is the sample size.  To keep the interpretation easy, use sqrt(n).

rd$sample.size.per.rep.sqrt <- sqrt(rd$sample.size.per.rep)

plot(exp(exponent) ~ sample.size.per.rep, rd)
plot(exp(exponent) ~ sample.size.per.rep.sqrt, rd)

fit<-lm(exp(exponent) ~ sample.size.per.rep.sqrt, rd); summary(fit)

# Because they might contribute, we might want to throw in predator type ('parasitoid' vs. 'predator') and whether or not eaten prey were replaced during the experiment ('expttype').  However, 'expttype' is confounded with 'sample size', so will cause problems. 


library(ggplot2)
p <- ggplot(rd,aes(sample.size.per.rep.sqrt,exp(exponent)))

p + geom_point(aes(shape=predator, color=expttype))

p + geom_point(aes(shape=expttype, color=predator))

fit<-lm(exp(exponent) ~ predator + expttype + sample.size.per.rep.sqrt * log(pred.ratio), rd); summary(fit)

fit<-lm(exp(exponent) ~ predator  + sample.size.per.rep.sqrt * log(pred.ratio), rd); summary(fit)


# We also have a problem distinguishing the effects of pred.ratio and sample size because they are negatively correlated across studies.
plot(rd$sample.size.per.rep.sqrt,log(rd$pred.ratio))
cor.test(rd$sample.size.per.rep.sqrt,log(rd$pred.ratio))

# Try again after removing the three 'atypically' high pred.ratio studies
kp <- which(log(rd$pred.ratio)<3.5)
plot(rd$sample.size.per.rep.sqrt[kp],log(rd$pred.ratio[kp]))
cor.test(rd$sample.size.per.rep.sqrt[kp],log(rd$pred.ratio[kp]))

fit<-lm(exp(exponent) ~ predator + expttype + sample.size.per.rep.sqrt * log(pred.ratio), rd, subset=kp); summary(fit)
