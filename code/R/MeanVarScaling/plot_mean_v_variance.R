library(ggplot2)
library(ggforce)

type <- 'One_Predator_One_Prey'
dat <- read.csv(paste0('../../../results/R/MeanVarScaling/AggregatedData/AggData_',
                       type, '.csv'))

pdf(paste0('../../../results/R/MeanVarScaling/DataPlot_',type,'.pdf'))
for( p in 1:3){
  print(ggplot(data = dat,
       mapping = aes(Nconsumed.mean,
                     Nconsumed.var,
                     colour = factor(Npredator))) +
         geom_abline(linetype="dashed", size=0.1) +
         geom_point() + 
         theme(legend.position="none",
               text = element_text(size=5),
               strip.text = element_text(size = 5)) +
         facet_wrap_paginate(~dataname, 
                             scales="free",
                             ncol = 4, nrow = 5, page = p))
}
dev.off() 


pdf(paste0('../../../results/R/MeanVarScaling/DataPlot_',type,'_log.pdf'))
for( p in 1:3){
  print(ggplot(data = dat,
               mapping = aes(Nconsumed.mean,
                             Nconsumed.var,
                             colour = factor(Npredator))) +
          coord_trans(x ='log10', y='log10') +
          geom_abline(linetype="dashed", size=0.1) +
          geom_point() + 
          theme(legend.position="none",
                text = element_text(size=5),
                strip.text = element_text(size = 5)) +
          facet_wrap_paginate(~dataname, 
                              scales="free",
                              ncol = 4, nrow = 5, page = p))
}
dev.off() 



type <- 'One_Predator_Two_Prey'
dat <- read.csv(paste0('../../../results/R/MeanVarScaling/AggregatedData/AggData_',
                       type, '.csv'))

pdf(paste0('../../../results/R/MeanVarScaling/DataPlot_',type,'.pdf'))
for( p in 1:3){
  print(ggplot(data = dat,
               mapping = aes(Nconsumed.MainPrey.mean, 
                             Nconsumed.MainPrey.var,
                             colour = factor(Nprey.AltPrey))) +
          geom_abline(linetype="dashed", size=0.1) +
          geom_point() + 
          theme(legend.position="none",
                text = element_text(size=5),
                strip.text = element_text(size = 5)) +
          facet_wrap_paginate(~dataname+MainPrey, 
                              scales="free",
                              ncol = 4, nrow = 5, page = p))
}
dev.off() 

pdf(paste0('../../../results/R/MeanVarScaling/DataPlot_',type,'_log.pdf'))
for( p in 1:3){
  print(ggplot(data = dat,
               mapping = aes(Nconsumed.MainPrey.mean, 
                             Nconsumed.MainPrey.var,
                             colour = factor(Nprey.AltPrey))) +
          coord_trans(x ='log10', y='log10') +
          geom_abline(linetype="dashed", size=0.1) +
          geom_point() + 
          theme(legend.position="none",
                text = element_text(size=5),
                strip.text = element_text(size = 5)) +
          facet_wrap_paginate(~dataname+MainPrey, 
                              scales="free",
                              ncol = 4, nrow = 5, page = p))
}
dev.off() 



