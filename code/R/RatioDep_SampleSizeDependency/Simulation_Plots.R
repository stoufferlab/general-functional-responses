# ~~~~~~~~~~~~
# Plot results of simulation
# ~~~~~~~~~~~~
load(file='output_Sim/Simulation_dat.Rdata')

# Rearrange models to enable use of facets
dat.AA <- dat[,c(1:14,15:18)]
colnames(dat.AA)[15:18] <- sub('AA.','',colnames(dat.AA)[15:18])
dat.AA$model <- 'AA'
dat.qAA <- dat[,c(1:14,19:22)]    
colnames(dat.qAA)[15:18] <- sub('qAA.','',colnames(dat.qAA)[15:18])
dat.qAA$model <- 'qAA'
dat2 <- rbind(dat.AA,dat.qAA)


ggplot(dat2,aes(R,exponent, colour=factor(mxP))) + geom_point() + geom_smooth(se=FALSE) + facet_grid(model~caseName) + coord_trans(x="log")
ggplot(dat2,aes(R,attack, colour=factor(mxP))) + geom_point() + geom_smooth(se=FALSE) + facet_grid(model~caseName) + coord_trans(x="log")
ggplot(dat2,aes(R,handling, colour=factor(mxP))) + geom_point() + geom_smooth(se=FALSE) + facet_grid(model~caseName) + coord_trans(x="log")

ggplot(dat2,aes(attack,exponent)) + geom_point(aes(colour=factor(mxP))) + facet_grid(model~caseName)
ggplot(dat2,aes(handling,exponent)) + geom_point(aes(colour=factor(mxP))) + facet_grid(model~caseName)
ggplot(dat2,aes(attack,handling)) + geom_point(aes(colour=factor(mxP))) + facet_grid(model~caseName)

