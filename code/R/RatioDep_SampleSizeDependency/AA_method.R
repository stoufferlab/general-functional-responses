# Implement 'Method II' of Arditi & Akcakaya 1990

# Assumes values of P represent discrete levels, with variation in N for each P 
# (i.e. a 'manipulative experiment'), rather than a set of continuous densities 

# source('../LogLikelihoods/holling_like_single_predator_single_prey_unbounded.R')


# I've made a change a few more changes now.

x<-model.matrix(~0+as.factor(NP$P))

a<-as.vector(1:ncol(x))

x%*%a



AA.method <- function(N,P,Neaten)
  Plevels <- sort(unique(P))



