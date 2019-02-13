source('../../lib/plot_coefs.R') # for plot_coefs() and order.of.fits()
source('../../lib/holling_method_one_predator_one_prey.R')
source('../../lib/ratio_method_one_predator_one_prey.R')

# load('../../../../results/R/ffr.fits_OnePredOnePrey.Rdata')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
fit.order <- order.of.fits(ffr.fits, order=TRUE, model="Arditi.Akcakaya", order.parm="Sample size")
ffr.fits <- ffr.fits[fit.order]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
SS <- unlist(lapply(ffr.fits, function(x) x$study.info$sample.size))

# Summarize mean of AIC estimates across bootstrapped fits 
AIC.H1 <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Holling.Type.I']],       function(y){AIC(logLik(y))}  )))   }))
AIC.H2 <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Holling.Type.II']],      function(y){AIC(logLik(y))}  )))   }))
AIC.R <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Ratio']],                function(y){AIC(logLik(y))}  )))   }))
AIC.AG <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Arditi.Ginzburg']],      function(y){AIC(logLik(y))}  )))   }))
AIC.HV <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Hassell.Varley']],       function(y){AIC(logLik(y))}  )))   }))
AIC.AA <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Arditi.Akcakaya']],      function(y){AIC(logLik(y))}  )))   }))
AIC.BD <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Beddington.DeAngelis']], function(y){AIC(logLik(y))}  )))   }))
AIC.CM <- unlist(lapply(ffr.fits, function(x){ 
  mean(unlist(lapply(  x$fits[['Crowley.Martin']],       function(y){AIC(logLik(y))}  )))   }))
  

AICs <- data.frame(AIC.H1, AIC.H2, AIC.R, AIC.AG, AIC.HV, AIC.AA, AIC.BD, AIC.CM)

rownames(AICs)<-sub('./Dataset_Code/','', rownames(AICs))
rownames(AICs)<-sub('.R','', rownames(AICs))
rownames(AICs) <- paste0(rownames(AICs), ' (',SS,')')
colnames(AICs) <- sub('AIC.', '', colnames(AICs))

minAICs <- apply(AICs, 1, min)

dAICs <- AICs - minAICs

rnkAICs <- t(apply(AICs, 1, order, decreasing=T))
colnames(rnkAICs) <- colnames(AICs)

#~~~~~~~~~~~
# Rank order
#~~~~~~~~~~~
pdf('../../../../results/R/OnePredOnePrey_AIC_ranks.pdf',height=6,width=5)
par(mar=c(3,10,1,1), mgp=c(1.5,0.1,0), tcl=-0.1, las=1, cex=0.7, yaxs='i')
  plot(1:nrow(rnkAICs), 1:nrow(rnkAICs),
       type='n', yaxt='n',
       xlim=c(1,ncol(rnkAICs)),
       ylim=c(0,nrow(rnkAICs)+1),
       xlab='Model rank by AIC',
       ylab='')
  axis(side=2, at=1:nrow(rnkAICs), labels=rownames(rnkAICs), cex.axis=0.5, las=2)
  
  Mcols <- rainbow(ncol(rnkAICs))
  for(m in 1:ncol(rnkAICs)){
    points(rnkAICs[,m], 1:nrow(rnkAICs), type='p', pch=21, col=Mcols[m], bg=Mcols[m])
  }  
  legend('bottomright',legend=colnames(rnkAICs),pch=21, pt.bg=Mcols, col=Mcols, bg='white')
dev.off()

#~~~~~
# dAIC
#~~~~~
plot(1:nrow(dAICs), 1:nrow(dAICs),
     type='n', yaxt='n',
     # xlim=c(0,max(dAICs)),
     xlim=c(0,4),
     xlab='Model delta-AIC',
     ylab=''
)
axis(side=2, at=1:nrow(dAICs), labels=rownames(dAICs), cex.axis=0.5, las=2)

Mcols <- rainbow(ncol(dAICs))
for(m in 1:ncol(dAICs)){
  points(dAICs[,m], 1:nrow(dAICs), pch=21, col=Mcols[m], bg=Mcols[m])
}  

legend('bottomright',legend=colnames(dAICs),pch=21, pt.bg=Mcols, col=Mcols)
