

#also run at least one shifting optima simulation

zoptP <- c(rep(0,times=100),rep(3,times=100),rep(2,times=100),rep(3,times=200)) 
zoptM <- c(rep(0,times=100),rep(2,times=200),rep(2,times=100),rep(3,times=100))
zvecs <- list(PlantZ = zoptP, MicrZ = zoptM)
testvarE <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,
		zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=length(zoptP),Lambda=25,mutprb=0.0005,
		prbHorz = 0.2,pfP = 0.6, pfM=0.6,
		FLFC=0.1,zoptvects =zvecs)


#save(testvarE,file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_testvarE.R",sep=""))

#load(paste(Sys.getenv("SCRATCH"),"/Simulation Results_testvarE.R",sep=""))

#figure, summarizing testvarE


VmVpVE <- extractVmVp(testvarE, 1,length(zoptP)+1,1)

#pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_varE_statsovertime.pdf",sep=""),height=9,width=6)
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims//Simulation_Results_varE_statsovertime.pdf",height=6,width=3)
        par(mfrow=c(2,1))
        par(mar=c(2,3,1,1))
        par(oma=c(1.5,0,0,0))
	windowplot(1,length(zoptP)+1,5,testvarE,ylim=c(-2,6),"",xlab="",ylab="")
        	abline(h=0)
           	lines(zoptP~c(2:(length(zoptP)+1)),lty=2,col=rgb(0,0.5,0))
           	lines(zoptM~c(2:(length(zoptP)+1)),lty=2,col=rgb(0.5,0,0.5))
#         	mtext("variable optima",side=3,line=0.5)  
        	mtext("breeding value",side=2,line=2)
		legend(0,y=6.5,c("Plant","Microbe","Interacting"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
		legend(225,y=6.5,c("Optima"), lty=2,  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
        plot(VmVpVE$Vp ~ c(1:length(VmVpVE$Vp)),pch=NA,ylim=c(0,0.12),ylab="Genetic variance",xlab="")
		  	mtext("genetic variance",side=2,line=2)
	       	mtext("generations",side=1,line=2)
	       	abline(v=c(101,201,301,401),lty=3,col=rgb(0,0,0,alpha=0.5))
			lines( 1:length(VmVpVE$PVp),VmVpVE$Vp,col=rgb(0,0.5,0))
			lines(1:length(VmVpVE$PVm),VmVpVE$Vm , col=rgb(0.5,0,0.5))
			lines( 1:length(VmVpVE$PVp),VmVpVE$Vb,col=rgb(0.5,0.5,0.5))
dev.off()

