#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to run simulations across specific scenarios for illustrative purposes
##

source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

gens <- 500
#test.4a <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=5,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)#at about .1 for fiterr relative to sdfit of 1 is when error in fitness starts to obscure fitness differences, not calculated. just apprx guess.	
#test.4b <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=5,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#test.4d <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=1,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#test.4e <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=1,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)

#fourdemosims <- list(test.4a,test.4b,test.4d,test.4e)
#save(fourdemosims,file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_PFF_expDFE.R",sep=""))
fourdmosims <- load(paste(Sys.getenv("SCRATCH"),"/Simulation Results_PFF_expDFE.R",sep=""))
test.4a <- fourdemosims[[1]]
test.4b <- fourdemosims[[2]]
test.4d <- fourdemosims[[3]]
test.4e <- fourdemosims[[4]]

titles <- c("~0 link to microbe fitness","~0 link to plant fitness","microbe link > plant link","plant link > microbe link")

pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_fourdemos.pdf",sep=""),height=6,width=6)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(oma=c(1,0,1,0))
windowplot(1,gens+1,5,test.4a,ylim=c(-2,6),titles[1])
	abline(h=3,col=rgb(0,0,1),lty=2)
	abline(h=0)
windowplot(1,gens+1, 5,test.4b,ylim=c(-2,6),titles[2])
	abline(h=3,col=rgb(1,0,0),lty=2)
	abline(h=0)
	legend(gens*-0.05,y=6,c("Plant","Microbe","Holobiont"), fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")
	legend(gens*0.35,y=6,c("Plant optima","Microbe optima","Midpoint"), lty = 2, col = c(rgb(0,0,1),rgb(1,0,0),rgb(0,0,0,alpha=.5)), bty="n")
windowplot(1,gens+1,5,test.4d,ylim=c(-2,6),titles[3])
	abline(h=2,col=rgb(0,0,1),lty=2); abline(h=3,col=rgb(1,0,0),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot(1,gens+1,5,test.4e,ylim=c(-2,6),titles[4])
	abline(h=3,col=rgb(0,0,1),lty=2); abline(h=2,col=rgb(1,0,0),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
dev.off()

#demo figure for total variance "G", variance G proportions of M and P, distance from optimum, local dynamics?
#take two simulations from above figure, and plot each below
FC4b <- getfitcon(10, gens+1, 1, test.4b,zoP=2,zoM=3, wP=5, wM=0.75,pfP=0.6,pfM=0.6)   
VmVp4b <- extractVmVp(test.4b, 1,gens+1,1)
win4b  <- extractwinning(test.4b,first=1,last=gens+1,1,zoP=2,zoM= 3)
dyn4b <- extractDyn(test.4b,first=1,last=gens,10)
FC4e <- getfitcon(10, gens+1, 1, test.4e,zoP=3,zoM=2, wP=0.75, wM=1,pfP=0.6,pfM=0.6)   
VmVp4e <- extractVmVp(test.4e, 1,gens+1,1)
win4e  <- extractwinning(test.4e,first=1,last=gens+1,1,zoP=3,zoM= 2)
dyn4e <- extractDyn(test.4e,first=1,last=gens,10)
pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_fourdemos_statsovertime.pdf",sep=""),height=12,width=6)
	par(mfrow=c(4,2))
	par(mar=c(3,5,1,1))
	par(oma=c(5,0,1,0))
 	plot(FC4b$fitnesscorrelation~c(10:(gens+1)),main="~0 link to plant fitness",ylab="Fitness correlation",xlab="",ylim=c(-1,1),pch=NA)
	       lines(FC4b$fitnesscorrelation~c(10:(gens+1)),lty=2)
		abline(h=0)
	plot(FC4e$fitnesscorrelation~c(10:(gens+1)),main="plant link > microbe link",ylab="",xlab="",ylim=c(-1,1),pch=NA)
		lines(FC4e$fitnesscorrelation~c(10:(gens+1)),lty=2)
		abline(h=0)
 	plot(VmVp4b$Vp ~ c(1:length(VmVp4b$Vp)),pch=NA,ylim=c(0,1),ylab="Proportion of variance",xlab="")
		lines( 1:length(VmVp4b$PVp),VmVp4b$PVp,col=rgb(0,0,1)) 
		lines(1:length(VmVp4b$PVm),VmVp4b$PVm , col=rgb(1,0,0)) 
	plot(VmVp4e$Vp ~ c(1:length(VmVp4e$Vp)),pch=NA,ylim=c(0,1),ylab="",xlab="")
		lines( 1:length(VmVp4e$PVp),VmVp4e$PVp,col=rgb(0,0,1))
		lines(1:length(VmVp4e$PVm),VmVp4e$PVm , col=rgb(1,0,0))
	plot(win4b$dP ~ c(1:length(win4b$dP)),pch=NA,ylim=c(0,max(c(win4b$dP,win4b$dM))),ylab="Distance from optima",xlab="")
		lines(win4b$dP~c(1:(gens+1)),col=rgb(0,0,1))
		lines(win4b$dM~c(1:(gens+1)),col=rgb(1,0,0))
	plot(win4e$dP ~ c(1:length(win4e$dP)),pch=NA,ylim=c(0,max(c(win4e$dP,win4e$dM))),ylab="",xlab="")
		lines(win4e$dP~c(1:(gens+1)),col=rgb(0,0,1))
		lines(win4e$dM~c(1:(gens+1)),col=rgb(1,0,0))
	plot(dyn4b$tcoefP ~ c(1:length(dyn4b$tcoefP)),pch=NA,ylim=c(min(c(dyn4b$tcoefP,dyn4b$tcoefM)),max(c(dyn4b$tcoefP,dyn4b$tcoefM))),ylab="Local trait change",xlab="generations")
		abline(h=0)
		lines(dyn4b$tcoefP~c(1:length(dyn4b$tcoefP)),col=rgb(0,0,1))
		lines(dyn4b$tcoefM~c(1:length(dyn4b$tcoefM)),col=rgb(1,0,0))	
	plot(dyn4e$tcoefP ~ c(1:length(dyn4e$tcoefP)),pch=NA,ylim=c(min(c(dyn4e$tcoefP,dyn4e$tcoefM)),max(c(dyn4e$tcoefP,dyn4e$tcoefM))),ylab="Local trait change",xlab="generations")
 		abline(h=0)
		lines(dyn4e$tcoefP~c(1:length(dyn4e$tcoefP)),col=rgb(0,0,1))
		lines(dyn4e$tcoefM~c(1:length(dyn4e$tcoefM)),col=rgb(1,0,0))
dev.off()

#also run at least one shifting optima simulation

zoptP <- c(rep(0,times=100),rep(3,times=100),rep(2,times=100),rep(3,times=200)) 
zoptM <- c(rep(0,times=100),rep(2,times=200),rep(2,times=100),rep(3,times=100))
zvecs <- list(PlantZ = zoptP, MicrZ = zoptM)
testvarE <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,
		zoP=3,zoM=2,wP=0.75,wM=1,timesteps=gens,Lambda=30,mutprb=0.0005,
		fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,
		FLFC=0.1,zoptvects =zvecs)

save(testvarE,file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_testvarE.R",sep=""))

#load(paste(Sys.getenv("SCRATCH"),"/Simulation Results_testvarE.R",sep=""))

#figure, summarizing testvarE

FCVE <- getfitcon(10, gens+1, 1, testvarE,zoP=3,zoM=2, wP=0.75, wM=1,pfP=0.6,pfM=0.6)
VmVpVE <- extractVmVp(testvarE, 1,gens+1,1)
winVE  <- extractwinning(testvarE,first=1,last=gens+1,1,zoP=3,zoM= 2,zoptvects=zvecs) #have zoptvects, so zoP and zoM don't matter in this case
dynVE <- extractDyn(testvarE,first=1,last=gens,10)

pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_varE_statsovertime.pdf",sep=""),height=9,width=6)
        par(mfrow=c(3,2))
        par(mar=c(3,5,1,1))
        par(oma=c(5,0,1,0))
	windowplot(1,gens+1,5,testvarE,ylim=c(-2,6),"variable optima, stronger link to plant fitness")
        	abline(h=0)
        plot(FCVE$fitnesscorrelation~c(10:(gens+1)),main="",ylab="Fitness correlation",xlab="",ylim=c(-1,1),pch=NA)
		lines(FCVE$fitnesscorrelation~c(10:(gens+1)),lty=2)
		abline(h=0)
        plot(VmVpVE$Vp ~ c(1:length(VmVpVE$Vp)),pch=NA,ylim=c(0,1),ylab="Proportion of genetic variance",xlab="")
		lines( 1:length(VmVpVE$PVp),VmVpVE$PVp,col=rgb(0,0,1))
		lines(1:length(VmVpVE$PVm),VmVpVE$PVm , col=rgb(1,0,0))
        plot(winVE$dP ~ c(1:length(winVE$dP)),pch=NA,ylim=c(0,max(c(winVE$dP,winVE$dM))),ylab="Distance from optima",xlab="")
		lines(winVE$dP~c(1:(gens+1)),col=rgb(0,0,1))
		lines(winVE$dM~c(1:(gens+1)),col=rgb(1,0,0))
        plot(dynVE$tcoefP ~ c(1:length(dynVE$tcoefP)),pch=NA,ylim=c(min(c(dynVE$tcoefP,dynVE$tcoefM)),max(c(dynVE$tcoefP,dynVE$tcoefM))),ylab="Local trait change",xlab="generations")
		abline(h=0)
		lines(dynVE$tcoefP~c(1:length(dynVE$tcoefP)),col=rgb(0,0,1))
		lines(dynVE$tcoefM~c(1:length(dynVE$tcoefM)),col=rgb(1,0,0))
	plot(zvecs$PlantZ~c(1:length(zvecs$PlantZ)),ylim=c(-2,6),ylab="trait optima",xlab="generations",pch=NA)
		lines(zvecs$PlantZ~c(1:length(zvecs$PlantZ)),col=rgb(0,0,1))
               lines(zvecs$MicrZ~c(1:length(zvecs$MicrZ)),col=rgb(1,0,0))
		abline(h=0)
dev.off()

#microbes > plants. quick test with low cost of free living and high cost of free living
#fitness cost comes before error
testmanyMa <- sim.cotrait(NP=100,NM=200,nlP=100,nlM=200,nlnP=3,nlnM=3,
                zoP=3,zoM=2,wP=0.75,wM=1,timesteps=gens,Lambda=30,mutprb=0.0005,
                fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,
                FLFC=0.1)#at time of running march 10, this is 10% less than the MINIMUM of symbiotic fitness. e.g. still a decent cost
testmanyMb <- sim.cotrait(NP=100,NM=200,nlP=100,nlM=200,nlnP=3,nlnM=3,
                zoP=3,zoM=2,wP=0.75,wM=1,timesteps=gens,Lambda=30,mutprb=0.0005,
                fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,
                FLFC=0.5)
manyM <- list(testmanyMa,testmanyMb)
save(manyM,file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_manymicrobes.R",sep=""))
#load(file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_manymicrobes.R",sep=""))
#testmanyMa <- manyM[[1]]
#testmanyMb <- manyM[[2]]

FCmMa <- getfitcon(10, gens+1, 1, testmanyMa,zoP=3,zoM=2, wP=0.75, wM=1,pfP=0.6,pfM=0.6)
VmVpmMa <- extractVmVp(testmanyMa, 1,gens+1,1)
FCmMb <- getfitcon(10, gens+1, 1, testmanyMb,zoP=3,zoM=2, wP=0.75, wM=1,pfP=0.6,pfM=0.6)
VmVpmMb <- extractVmVp(testmanyMb, 1,gens+1,1)


pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_manymicrobes_statsovertime.pdf",sep=""),height=9,width=6)
        par(mfrow=c(3,2))
        par(mar=c(3,5,1,1))
        par(oma=c(5,0,1,0))
        windowplot(1,gens+1,5,testmanyMa,ylim=c(-2,6),"low cost of free-living, stronger link to plant fitness")
                abline(h=3,col=rgb(0,0,1),lty=2)
                abline(h=0)
        windowplot(1,gens+1,5,testmanyMb,ylim=c(-2,6),"high cost of free-living, stronger link to plant fitness")
                abline(h=3,col=rgb(0,0,1),lty=2)
                abline(h=0)
	 plot(FCmMa$fitnesscorrelation~c(10:(gens+1)),main="",ylab="Fitness correlation",xlab="",ylim=c(-1,1),pch=NA)
                lines(FCmMa$fitnesscorrelation~c(10:(gens+1)),lty=2)
                abline(h=0)
        plot(FCmMb$fitnesscorrelation~c(10:(gens+1)),main="",ylab="",xlab="",ylim=c(-1,1),pch=NA)
                lines(FCmMb$fitnesscorrelation~c(10:(gens+1)),lty=2)
                abline(h=0)
       plot(VmVpmMa$Vp ~ c(1:length(VmVpmMa$Vp)),pch=NA,ylim=c(0,1),ylab="Proportion of genetic variance",xlab="generations")
                lines( 1:length(VmVpmMa$PVp),VmVpmMa$PVp,col=rgb(0,0,1))
                lines(1:length(VmVpmMa$PVm),VmVpmMa$PVm , col=rgb(1,0,0))
       plot(VmVpmMb$Vp ~ c(1:length(VmVpmMb$Vp)),pch=NA,ylim=c(0,1),ylab="",xlab="generations")
                lines( 1:length(VmVpmMb$PVp),VmVpmMb$PVp,col=rgb(0,0,1))
                lines(1:length(VmVpmMb$PVm),VmVpmMb$PVm , col=rgb(1,0,0))
dev.off()


