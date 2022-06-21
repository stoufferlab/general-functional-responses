
# remove all prior fits from the workspace to avoid fit of prior dataset
# being used when a fit on the current dataset fails
rm(list=ls(pattern='^fit.'))

source('lib/ttest.R') 

# Use no weights (i.e. equal weights) or weight treatments
# (i.e. mean-variance estimates) by each treatment's number of replicates
# (only relevant when replication varies across treatments)
# 'weightbyrepl' is specified is sourcing parent script
dat$weights <- if(weightbyrepl){ dat$n }else{ rep(1, length(dat$n)) }

dat1 <- subset(dat, MainPrey=='Prey1')
dat2 <- subset(dat, MainPrey=='Prey2')

# Almost all datasets converge with glm's default maxit = 25,
# a few require maxit = 35, and one (Elliot_2006_i5) requires 100,
# so just let them all have up to 100 iterations!
maxit <- 100 
  
# linear models
fit.lm.1 <- lm(log.Nconsumed.MainPrey.var ~
                 log.Nconsumed.MainPrey.mean,
             data = dat1,
             weights = dat1$weights)
fit.lm.2 <- lm(log.Nconsumed.MainPrey.var ~
                 log.Nconsumed.MainPrey.mean,
               data = dat2,
               weights = dat2$weights)

# generalized linear model
fit.glm.1 <- glm(Nconsumed.MainPrey.var ~ 
                   log.Nconsumed.MainPrey.mean,
                 data = dat1,
                 weights = dat1$weights,
                 family = gaussian(link = log),
                 control = list(maxit = maxit))
fit.glm.2 <- glm(Nconsumed.MainPrey.var ~ 
                   log.Nconsumed.MainPrey.mean,
                 data = dat2,
                 weights = dat2$weights,
                 family = gaussian(link = log),
                 control = list(maxit = maxit))

# w/ other-prey effect
# linear model
fit.lm.main.1 <- lm(log.Nconsumed.MainPrey.var ~
                         log.Nconsumed.MainPrey.mean + log.Nprey.AltPrey,
                 data = dat1,
                 weights = dat1$weights)
fit.lm.main.2 <- lm(log.Nconsumed.MainPrey.var ~
                         log.Nconsumed.MainPrey.mean + log.Nprey.AltPrey,
                 data = dat2,
                 weights = dat2$weights)

# generalized linear model
fit.glm.main.1 <- glm(Nconsumed.MainPrey.var ~ 
                      log.Nconsumed.MainPrey.mean + log.Nprey.AltPrey,
                   data = dat1,
                   weights = dat1$weights,
                   family = gaussian(link = log),
                   control = list(maxit = maxit))
fit.glm.main.2 <- glm(Nconsumed.MainPrey.var ~ 
                           log.Nconsumed.MainPrey.mean + log.Nprey.AltPrey,
                   data = dat2,
                   weights = dat2$weights,
                   family = gaussian(link = log),
                   control = list(maxit = maxit))

#...which is equivalent to...
fit.glm.int.1 <- glm(Nconsumed.MainPrey.var ~ 
                                 log.Nconsumed.MainPrey.mean + 
                                 I(log.Nprey.AltPrey*log.Nconsumed.MainPrey.mean),
                         data = dat1,
                         weights = dat1$weights,
                         family = gaussian(link = log),
                         control = list(maxit = maxit))
fit.glm.int.2 <- glm(Nconsumed.MainPrey.var ~ 
                                 log.Nconsumed.MainPrey.mean + 
                                 I(log.Nprey.AltPrey*log.Nconsumed.MainPrey.mean),
                         data = dat2,
                         weights = dat2$weights,
                         family = gaussian(link = log),
                         control = list(maxit = maxit))

# ...for which the corresponding linear model is...
fit.lm.int.1 <- lm(log.Nconsumed.MainPrey.var ~
                               log.Nconsumed.MainPrey.mean + 
                               log.Nprey.AltPrey:log.Nconsumed.MainPrey.mean,
                       data = dat1,
                       weights = dat1$weights)
fit.lm.int.2 <- lm(log.Nconsumed.MainPrey.var ~
                               log.Nconsumed.MainPrey.mean + 
                               log.Nprey.AltPrey:log.Nconsumed.MainPrey.mean,
                       data = dat2,
                       weights = dat2$weights)

# Rename the coefficients (causes problems for glm.confint() )
names(fit.lm.1$coefficients)[1] <-
  names(fit.lm.2$coefficients)[1] <-
  names(fit.lm.main.1$coefficients)[1] <-
  names(fit.lm.main.2$coefficients)[1] <-
  names(fit.lm.int.1$coefficients)[1] <-
  names(fit.lm.int.2$coefficients)[1] <-
  'log.a'
names(fit.glm.1$coefficients)[1] <-
  names(fit.glm.2$coefficients)[1] <-
  names(fit.glm.main.1$coefficients)[1] <-
  names(fit.glm.main.2$coefficients)[1] <-
  names(fit.glm.int.1$coefficients)[1] <-
  names(fit.glm.int.2$coefficients)[1] <-
  'a'
names(fit.lm.1$coefficients)[2] <-
  names(fit.lm.2$coefficients)[2] <-
  names(fit.lm.main.1$coefficients)[2] <-
  names(fit.lm.main.2$coefficients)[2] <-
  names(fit.lm.int.1$coefficients)[2] <-
  names(fit.lm.int.2$coefficients)[2] <-
  names(fit.glm.1$coefficients)[2] <-
  names(fit.glm.2$coefficients)[2] <-
  names(fit.glm.main.1$coefficients)[2] <-
  names(fit.glm.main.2$coefficients)[2] <-
  names(fit.glm.int.1$coefficients)[2] <-
  names(fit.glm.int.2$coefficients)[2] <-
  'b'
names(fit.lm.main.1$coefficients)[3] <-
  names(fit.lm.main.2$coefficients)[3] <-
  names(fit.glm.main.1$coefficients)[3] <-
  names(fit.glm.main.2$coefficients)[3] <-
  'c'
names(fit.lm.int.1$coefficients)[3] <-
  names(fit.lm.int.2$coefficients)[3] <-
  names(fit.glm.int.1$coefficients)[3] <-
  names(fit.glm.int.2$coefficients)[3] <-
  'd'


# aggregate the fits
fits <- list(lm = list(prey1 = fit.lm.1,
                       prey2 = fit.lm.2),
             lm.main = list(prey1 = fit.lm.main.1,
                            prey2 = fit.lm.main.2),
             lm.int = list(prey1 = fit.lm.int.1,
                           prey2 = fit.lm.int.2),
             glm = list(prey1 = fit.glm.1,
                        prey2 = fit.glm.2),
             glm.main = list(prey1 = fit.glm.main.1,
                             prey2 = fit.glm.main.2),
             glm.int = list(prey1 = fit.glm.int.1,
                            prey2 = fit.glm.int.2)
             )


#  test null hypothesis of slope b (2nd parameter) being equal to 1 (function default)
ttest <- list(lm = list(prey1 = ttest(fit.lm.1, 2),
                        prey2 = ttest(fit.lm.2, 2)),
                lm.main = list(prey1 = ttest(fit.lm.main.1, 2),
                               prey2 = ttest(fit.lm.main.2, 2)),
                lm.int = list(prey1 = ttest(fit.lm.int.1, 2),
                              prey2 = ttest(fit.lm.int.2, 2)),
                glm = list(prey1 = ttest(fit.glm.1, 2),
                           prey2 = ttest(fit.glm.2, 2)),
                glm.main = list(prey1 = ttest(fit.glm.main.1, 2),
                                prey2 = ttest(fit.glm.main.2, 2)),
                glm.int = list(prey1 = ttest(fit.glm.int.1, 2),
                               prey2 = ttest(fit.glm.int.2, 2))
)
