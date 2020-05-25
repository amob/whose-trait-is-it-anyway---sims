

# load(file=paste(Sys.getenv("SCRATCH"),'/popsetIBD_smmu_confl_randflo.RData',sep=""))

load(file=paste(Sys.getenv("SCRATCH"),'/expset_abo.RData',sep=""))

expdat <- expsetabo$expdat

plantmeans <- tapply(expdat$trait, expdat$popP,mean)
micrmeans <- tapply(expdat$trait, expdat$popM,mean)

cor(plantmeans,micrmeans)

pdf(file=paste(Sys.getenv("SCRATCH"),'/ExperimentPhenores.pdf',sep=""),width=4,height=4)
plot(plantmeans~micrmeans,ylab="Population mean plant value",xlab="Population mean microbe value",main=round(cor(plantmeans,micrmeans),digits=3))
dev.off()

load(file=paste(Sys.getenv("SCRATCH"),'/expset_aba.RData',sep=""))

expdat <- expsetaba$expdat

plantmeans <- tapply(expdat$trait, expdat$popP,mean)
micrmeans <- tapply(expdat$trait, expdat$popM,mean)

pdf(file=paste(Sys.getenv("SCRATCH"),'/ExperimentPhenores_aba.pdf',sep=""),width=4,height=4)
plot(plantmeans~micrmeans,ylab="Population mean plant value",xlab="Population mean microbe value",main=round(cor(plantmeans,micrmeans),digits=3))
dev.off()