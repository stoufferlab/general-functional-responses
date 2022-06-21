

library(devtools)
# devtools::install_github("thomas-fung/mpcmp")
library(mpcmp) # For Conway-Maxwell-Poisson rcomp()

reps.n <- 4
Nmax <- 50
Ns <- seq(1, Nmax, 10)
Ps <- seq(1, 2, 1)

dat <- expand.grid(Nprey = rep(Ns, each = reps.n),
                    Npredator = Ps)

fr <- function(N, P, 
               a = 2, h = 0.05, c = 0.1){
  (a*N*P)/(1 + a*h*N + c*(P-1))
}

fr.stoch <- function(N, P){
  rcomp(n = 1,
        mu = fr(N = N, P = P), # mean
        nu = 1 # dispersion
)}

fr.stoch.rep <- function(dat){

    dat$Nconsumed <- mapply(fr.stoch, N = dat$Nprey, P = dat$Npredator)
    
    dsum <- as.data.frame(as.matrix(
      aggregate(list(Nconsumed = dat$Nconsumed), 
                      by = list(Npredator = dat$Npredator,
                                Nprey = dat$Nprey),
                      FUN = function(x){ c(
                        mean = mean(x),
                        var = var(x),
                        se = sd(x)/length(x),
                        n = length(x)
                      )})))
    
    fit <- lm(log(dsum$Nconsumed.var) ~ log(dsum$Nconsumed.mean))
    
    # par(mfrow = c(2, 1))
    
    # plot(dat$Nprey, dat$Nconsumed)
    # for(p in seq_along(sort(unique(dat$Npredator)))){
    #   y <- function(x){fr(x, P = p)}
    #   curve(y, 0, Nmax, add = TRUE)
    # }
    
    # plot(log(dsum$Nconsumed.var) ~ log(dsum$Nconsumed.mean))
    # abline(fit)
    
    out <- coef(fit)[2]
    names(out) <- 'b'
    
    return(out)
}

y <- function(n, x){ replicate(n, fr.stoch.rep(dat))}

y(10, dat)
