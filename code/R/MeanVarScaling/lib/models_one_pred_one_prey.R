
# remove all prior fits from the workspace to avoid fit of prior dataset
# being used when a fit on the current dataset fails
rm(list=ls(pattern='^fit.'))

source('lib/ttest.R') 

# Use no weights (i.e. equal weights) or weight treatments
# (i.e. mean-variance estimates) by each treatment's number of replicates
# (only relevant when replication varies across treatments)
# 'weightbyrepl' is specified is sourcing parent script
weights <- if(weightbyrepl){ dat$n }else{ rep(1, length(dat$n)) }

# nls included to confirm equivalence to glm

# Almost all datasets converge with default maxit = 25,
# a few require maxit = 35, and some require 100,
# so just let them all have up to 100 iterations!
maxit <- 200

#~~~~~~~~~~~~~~~~~
# Power-law models
#~~~~~~~~~~~~~~~~~

# linear model
fit.lm <- lm(log.Nconsumed.var ~
               log.Nconsumed.mean,
             data = dat,
             weights = weights)

# generalized linear model
fit.glm <- glm(Nconsumed.var ~ 
                 log.Nconsumed.mean,
               data = dat,
               weights = weights,
               family = gaussian(link = log))

# nonlinear least squares (gives same estimates as glm, as expected)
# fit.nls <- nls(Nconsumed.var ~
#                  a * Nconsumed.mean ^ b,
#                data = dat,
#                weights = weights,
#                start = list(a = exp(coef(fit.glm)[1]),
#                             b = coef(fit.glm)[2]),
#                control = list(maxiter = maxit))

# quadratic model
fit.lm.q <- lm(log.Nconsumed.var ~
                 log.Nconsumed.mean + I(log.Nconsumed.mean^2),
               data = dat,
               weights = weights)

# w/ Predator effect
# linear model
fit.lm.main <- lm(log.Nconsumed.var ~
                    log.Nconsumed.mean + log.Npredator,
                  data = dat,
                  weights = weights)

# generalized linear model
fit.glm.main <- glm(Nconsumed.var ~ 
                      log.Nconsumed.mean + log.Npredator,
                    data = dat,
                    weights = weights,
                    family = gaussian(link = log))

# nls (again gives the same estimates as glm)
# fit.nls.main <- nls(Nconsumed.var ~
#                       a * Nconsumed.mean ^ b * Npredator ^ c,
#                     data = dat,
#                     weights = weights,
#                     start = list(a = coef(fit.nls)[1],
#                                  b = coef(fit.nls)[2],
#                                  c = 0),
#                  control = list(maxiter = maxit))

# # alternative predator effect formulation...
# # fit.nls.int <- nls(Nconsumed.var ~
# #                          a * Nconsumed.mean ^ (b + d * Npredator),
# #                        data = dat,
# #                        weights = weights,
# #                        start = list(a = coef(fit.nls)[1],
# #                                     b = coef(fit.nls)[2],
# #                                     d = 0),
# #                        control = list(maxiter = maxit))
# 
# #...which is equivalent to...
# fit.glm.int <- glm(Nconsumed.var ~ 
#                             log.Nconsumed.mean + I(Npredator*log.Nconsumed.mean),
#                           data = dat,
#                           weights = weights,
#                           family = gaussian(link = log))
# 
# # ...for which the corresponding linear model is...
# fit.lm.int <- lm(log.Nconsumed.var ~
#                           log.Nconsumed.mean + Npredator:log.Nconsumed.mean,
#                      data = dat,
#                      weights = weights)




# a different alternative predator effect formulation...
# fit.nls.int <- nls(Nconsumed.var ~
#                          a * Nconsumed.mean ^ (b + d * log.Npredator),
#                        data = dat,
#                        weights = weights,
#                        start = list(a = coef(fit.nls)[1],
#                                     b = coef(fit.nls)[2],
#                                     d = 0),
#                        control = list(maxiter = maxit))

#...which is equivalent to...
fit.glm.int <- glm(Nconsumed.var ~ 
                         log.Nconsumed.mean + I(log.Npredator*log.Nconsumed.mean),
                       data = dat,
                       weights = weights,
                       family = gaussian(link = log))

# ...for which the corresponding linear model is...
fit.lm.int <- lm(log.Nconsumed.var ~
                       log.Nconsumed.mean + log.Npredator:log.Nconsumed.mean,
                     data = dat,
                     weights = weights)


#~~~~~~~~~~~~~~~~~~
# Polynomial models
#~~~~~~~~~~~~~~~~~~

fit.mean <- lm(Nconsumed.var ~ 1,
               data = dat,
               weights = weights)

fit.quad <- lm(Nconsumed.var ~
                 -1 +
                 poly(Nconsumed.mean, 2, raw = TRUE),
               data = dat,
               weights = weights)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rename the coefficients (causes problems for glm.confint() )
names(fit.lm$coefficients)[1] <-
  names(fit.lm.q$coefficients)[1] <-
  names(fit.lm.main$coefficients)[1] <-
  names(fit.lm.int$coefficients)[1] <-
  'log.a'
names(fit.glm$coefficients)[1] <-
  names(fit.glm.main$coefficients)[1] <-
  names(fit.glm.int$coefficients)[1] <-
  'a'
names(fit.lm$coefficients)[2] <-
  names(fit.lm.q$coefficients)[2] <-
  names(fit.lm.main$coefficients)[2] <-
  names(fit.lm.int$coefficients)[2] <-
  names(fit.glm$coefficients)[2] <-
  names(fit.glm.main$coefficients)[2] <-
  names(fit.glm.int$coefficients)[2] <-
  'b'
names(fit.lm.main$coefficients)[3] <-
  names(fit.glm.main$coefficients)[3] <-
  'c'
names(fit.lm.int$coefficients)[3] <-
  names(fit.glm.int$coefficients)[3] <-
  'd'

names(fit.mean$coefficients) <- 'b0'
names(fit.quad$coefficients) <- c('b1', 'b2')

# aggregate the fits
fits <- list(lm = fit.lm,
             lm.q = fit.lm.q,
             lm.main = fit.lm.main,
             lm.int = fit.lm.int,
             glm = fit.glm,
             glm.main = fit.glm.main,
             glm.int = fit.glm.int,
             mean = fit.mean,
             quad = fit.quad)

# test null hypothesis of slope b (2nd parameter) being equal to 1 (function default)
ttest <- list(lm = ttest(fit.lm, 2),
              lm.q = ttest(fit.lm.q, 2),
              lm.main = ttest(fit.lm.main, 2),
              lm.int = ttest(fit.lm.int, 2),
              glm = ttest(fit.glm, 2),
              glm.main = ttest(fit.glm.main, 2),
              glm.int = ttest(fit.glm.int, 2))


