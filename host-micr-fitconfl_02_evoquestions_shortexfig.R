#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to produce illustrative graphics
	#run simulations across specific scenarios for illustrative purposes
##

# source('~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/host-micr-fitconfl_01_simfunction.R') 
 ## possible local directories are included, but commented out, as running long simulations on a personal laptop is not advisable
source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 
# 
# ###THESE SIMULATIONS ARE THE ONES THAT GET RUN
# # #Temporary commented out, because they take a very long time to run on a personal laptop and this is not recommended.
gens <- 300
# #4 scenarios
# #micr direct very weak 
# test.4b <- 		sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=3,wP=0.75,wM=10,  timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
# # micr direct very weak, indirect stronger
# test.4bff <- 	sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=3,wP=0.75,wM=10,  timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
# # conflict, links =,  direct only
# test.4f <- 		sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
# #conflict, links =,  direct and indirect links
# test.4fff <- 	sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1) 
# #matched optima, no ff
# test.4a <-		sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)#
# #mattched optima, ff
# test.4aff <- 	sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
# # sixdemosims <- list(test.4b,test.4bff,test.4a,test.4aff,test.4f,test.4fff)
# # save(sixdemosims,file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_sixdemos.RData',sep=""))
# 
# 
# load(file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_sixdemos.RData',sep=""))
# test.4b <- sixdemosims[[1]]
# test.4bff <- sixdemosims[[2]]
# test.4a <- sixdemosims[[3]]
# test.4aff <- sixdemosims[[4]]
# test.4f <- sixdemosims[[5]]
# test.4fff <- sixdemosims[[6]]

#additional simulation for supplement
#host direct very weak
# test.4c <- 		sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=3,wP=10,wM=0.75,  timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
# host direct very weak, indirect stronger
# save(test.4c,file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4c.RData',sep=""))
# test.4cff <- 	sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=3,zoM=3,wP=10,wM=0.75,  timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
# save(test.4cff,file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4cff.RData',sep=""))
# load(file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4c.RData',sep=""))
# load(file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4cff.RData',sep=""))

# # conflict, links =,  direct only, swap zopt
# test.4e <- 		sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=2,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
# save(test.4e,file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4e.RData',sep=""))
# #conflict, links =,  direct and indirect links
test.4eff <- 	sim.cotrait(NP=2000,NM=2000,nlP=20,nlM=40,nlnP=400,nlnM=800,zoP=2,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1) 
# save(test.4eff,file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4eff.RData',sep=""))
load(file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4e.RData',sep=""))
# load(file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_hostweaksims4eff.RData',sep=""))

# 
# pdf(paste(Sys.getenv("SCRATCH"), "/Simulation_Results_sixdemos_n.pdf",sep=""),height=5,width=7)
# # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/Simulation_Results_fourdemos_n.pdf",height=5,width=5)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# windowplot(1,gens+1, 3,test.4b,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[3])
# 	abline(h=3,col=rgb(0,0.5,0),lty=2,lwd=3); # abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
# 	abline(h=0)
# 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# 	mtext("breeding values",side=2,line=2,adj=-1.25)
# windowplot(1,gens+1, 3,test.4bff,ylim=c(-3,6),"",xlabs="",ylabs="")#"titles[4]")
# 	abline(h=3,col=rgb(0,0.5,0),lty=2,lwd=3); #abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	legend(gens*-0.05,y=0,c("Host","Microbe","Expressed"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0.5,0.5,0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
# windowplot(1,gens+1, 3,test.4a,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[3])
# 	abline(h=3,col=rgb(0,0,0),lty=2,lwd=3); # abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
# 	abline(h=0)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# 	mtext("Same optima",side=3,line=2)
# windowplot(1,gens+1, 3,test.4aff,ylim=c(-3,6),"",ylabs="")#"titles[4]")
# 	abline(h=3,col=rgb(0,0,0),lty=2,lwd=3); #abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
# 	abline(h=0)
# windowplot(1,gens+1,3,test.4f,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[9])
# 	abline(h=3,col=rgb(0,0.5,0),lty=2,lwd=3); abline(h=2,col=rgb(0.5,0,0.5),lty=2,lwd=3); abline(h=2.5,col=rgb(0,0,0),lty=2,lwd=3) 
# 	abline(h=0)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# 	mtext("Optima differ",side=3,line=2)
# windowplot(1,gens+1,3,test.4fff,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[10])
# 	abline(h=3,col=rgb(0,0.5,0),lty=2,lwd=3); abline(h=2,col=rgb(0.5,0,0.5),lty=2,lwd=3); abline(h=2.5,col=rgb(0,0,0),lty=3,lwd=3) 
# 	abline(h=0)
# 	legend(gens*-0.05,y=0,c("Host optimum","Microbe optimum","Both optima","Midpoint"), lty = c(2,2,2,3),lwd=2, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0),rgb(0,0,0)), bty="n")
# dev.off()
# 
# # #EXRA FIGURE
# pdf(paste(Sys.getenv("SCRATCH"), "/Simulation_Results_extrademos_n.pdf",sep=""),height=5.75,width=3.1)
# layout(matrix(1:2,ncol=1,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# windowplot(1,gens+1, 3,test.4c,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[3])
# 	abline(h=3,col=rgb(0.5,0,0.5),lty=2,lwd=2)
# 	abline(h=0)
# 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("breeding values",side=2,line=2,adj=-1)
# 	legend(gens*-0.05,y=0,c("Microbe optimum"), lty = 2,lwd=2, col = rgb(0.5,0,0.5), bty="n")
# windowplot(1,gens+1, 3,test.4cff,ylim=c(-3,6),"",xlabs="",ylabs="")#"titles[4]")
# 	abline(h=3,col=rgb(0.5,0,0.5),lty=2,lwd=2)
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	legend(gens*-0.05,y=0,c("Host","Microbe","Expressed"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
# 	mtext("generations",side=1,line=2)
# dev.off()
# 

# 
# 
# VmVp4b <- extractVmVp(test.4b, 1,gens+1,1)
# VmVp4bff <- extractVmVp(test.4bff, 1,gens+1,1)
# VmVp4a <- extractVmVp(test.4a, 1,gens+1,1)
# VmVp4aff <- extractVmVp(test.4aff, 1,gens+1,1)
# VmVp4f <- extractVmVp(test.4f, 1,gens+1,1)
# VmVp4fff <- extractVmVp(test.4fff, 1,gens+1,1)
# 
# win4b  <- extractwinning(test.4b,first=1,last=gens+1,1,zoP=3,zoM= 2) # some are sort of obvious..., but the rest matter
# win4bff  <- extractwinning(test.4bff,first=1,last=gens+1,1,zoP=3,zoM= 2)
# win4a  <- extractwinning(test.4a,first=1,last=gens+1,1,zoP=3,zoM= 3) # 
# win4aff  <- extractwinning(test.4aff,first=1,last=gens+1,1,zoP=3,zoM= 3)
# win4f  <- extractwinning(test.4f,first=1,last=gens+1,1,zoP=3,zoM= 2)
# win4fff  <- extractwinning(test.4fff,first=1,last=gens+1,1,zoP=3,zoM= 2)
#  
#  
# pdf(paste(Sys.getenv("SCRATCH"), "/Variance_Results_sixdemos_n.pdf",sep=""),height=5,width=7)
# # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/Variance_Results_fourdemos_n.pdf",height=4.5,width=5.5)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# 	plot(VmVp4b$Vp ~ c(1:length(VmVp4b$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4b$Vp),VmVp4b$Vp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4b$Vm),VmVp4b$Vm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4b$Vb),VmVp4b$Vb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4b$dP==min(win4b$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does host mean reach its optima or beyond
# 	 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 		mtext("No direct link to microbe fitness",side=3,line=0.5)
# 		mtext("genetic variance",side=2,line=2,adj=-1.5)
# # 		text(x=100,y=0.11, labels =expression(Z[opt[P]]~reached))
# 	plot(VmVp4bff$Vp ~ c(1:length(VmVp4bff$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4bff$Vp),VmVp4bff$Vp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4bff$Vm),VmVp4bff$Vm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4bff$Vb),VmVp4bff$Vb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4bff$dP==min(win4bff$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does host mean reach its optima or beyond
# 		mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# # 		text(x=95,y=0.11, labels =expression(Z[opt[P]]~reached))
# 	plot(VmVp4a$Vp ~ c(1:length(VmVp4a$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4a$Vp),VmVp4a$Vp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4a$Vm),VmVp4a$Vm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4a$Vb),VmVp4a$Vb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4a$dP==min(win4a$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does host mean reach its optima or beyond
# 		mtext("Equal links to fitness",side=3,line=0.5)
# 		mtext("Same optima",side=3,line=2)
# 	plot(VmVp4aff$Vp ~ c(1:length(VmVp4aff$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4aff$Vp),VmVp4aff$Vp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4aff$Vm),VmVp4aff$Vm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4aff$Vb),VmVp4aff$Vb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4aff$dP==min(win4aff$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does host mean reach its optima or beyond
# # 		text(x=95,y=0.11, labels =expression(Z[opt[P]]~reached))
# 		mtext("generations",side=1,line=2)
# 	plot(VmVp4f$Vp ~ c(1:length(VmVp4f$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4f$Vp),VmVp4f$Vp,col=rgb(0,0.5,0))
# 		lines(1:length(VmVp4f$Vm),VmVp4f$Vm , col=rgb(0.5,0,0.5))
# 		lines(1:length(VmVp4f$Vb),VmVp4f$Vb , col=rgb(0.5,0.5,0.5))
# 		abline(v=which(win4f$dM==min(win4f$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# # 		text(x=55,y=0.11, labels =expression(Z[opt[M]]~reached))
# 		mtext("Equal links to fitness",side=3,line=0.5)
# 		mtext("Different optima",side=3,line=2)
# 		legend(gens*0.4,y=0.10,c(expression(V[H]),expression(V[M]),expression(V[A])),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
# 	plot(VmVp4fff$Vp ~ c(1:length(VmVp4fff$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4fff$Vp),VmVp4fff$Vp,col=rgb(0,0.5,0))
# 		lines(1:length(VmVp4fff$Vm),VmVp4fff$Vm , col=rgb(0.5,0,0.5))
# 		lines(1:length(VmVp4fff$Vb),VmVp4fff$Vb , col=rgb(0.5,0.5,0.5))
# 		abline(v=which(win4fff$dM==min(win4fff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# # 		text(x=60,y=0.11, labels =expression(Z[opt[M]]~reached))
# dev.off()

##EXTRA FIGURE
# VmVp4c <- extractVmVp(test.4c, 1,gens+1,1)
# VmVp4cff <- extractVmVp(test.4cff, 1,gens+1,1)
# win4c  <- extractwinning(test.4c,first=1,last=gens+1,1,zoP=3,zoM= 3)
# win4cff  <- extractwinning(test.4cff,first=1,last=gens+1,1,zoP=3,zoM= 3)
# pdf(paste(Sys.getenv("SCRATCH"), "/Variance_Results_extrademos_n.pdf",sep=""),height=5,width=3)
# layout(matrix(1:2,ncol=1,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# 	plot(VmVp4c$Vp ~ c(1:length(VmVp4c$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4c$Vp),VmVp4c$Vp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4c$Vm),VmVp4c$Vm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4c$Vb),VmVp4c$Vb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4c$dM==min(win4c$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# 	 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 		mtext("No direct link to host fitness",side=3,line=0.5)
# 		mtext("genetic variance",side=2,line=2,adj=-1.75)
# 	plot(VmVp4cff$Vp ~ c(1:length(VmVp4cff$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 		lines( 1:length(VmVp4cff$Vp),VmVp4cff$Vp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4cff$Vm),VmVp4cff$Vm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4cff$Vb),VmVp4cff$Vb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4cff$dM==min(win4cff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# 		mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 		mtext("generations",side=1,line=2)
# 		legend(gens*0.4,y=0.10,c(expression(V[H]),expression(V[M]),expression(V[A])),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
# dev.off()

# rangep <- range(c(as.vector(test.4b$Plant),as.vector(test.4bff$Plant),as.vector(test.4a$Plant), as.vector(test.4aff$Plant), as.vector(test.4f$Plant),as.vector(test.4fff$Plant)))
# rangem <- range(c(as.vector(test.4b$Microbe),as.vector(test.4bff$Microbe),as.vector(test.4a$Microbe),as.vector(test.4aff$Microbe),as.vector(test.4f$Microbe),as.vector(test.4fff$Microbe)))
# maxabsp <- max(abs(rangep))
# maxabsm <- max(abs(rangem))

 

##trajectories
# getting allele specific information

# pdf(paste(Sys.getenv("SCRATCH"), "/AllTraject_plant_six.pdf",sep=""),height=5,width=7)
# # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/AllTraject_plant.pdf",width=8,height=5)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# plottraj(test.4b$Plant,maxpos=rangep[2],maxneg=rangep[1])
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# 	mtext("Frequency",side=2,line=2,adj=-0.6)
# 	abline(v=which(win4b$dP==min(win4b$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4bff$Plant,maxpos=rangep[2],maxneg=rangep[1])
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	abline(v=which(win4bff$dP==min(win4bff$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4a$Plant,maxpos=rangep[2],maxneg=rangep[1])
# 	mtext("Same optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# 	abline(v=which(win4a$dP==min(win4a$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4aff$Plant,maxpos=rangep[2],maxneg=rangep[1])
# 	mtext("generations",side=1,line=2)
# 	abline(v=which(win4aff$dP==min(win4aff$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4f$Plant,maxpos=rangep[2],maxneg=rangep[1])
# 	mtext("Different optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# 	abline(v=which(win4f$dM==min(win4f$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# plottraj(test.4fff$Plant,maxpos=rangep[2],maxneg=rangep[1])
# 	abline(v=which(win4fff$dM==min(win4fff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# dev.off()

## EXTRA FIGURE
# rangepc <- range(c(as.vector(test.4c$Plant),as.vector(test.4cff$Plant)))
# rangemc <- range(c(as.vector(test.4c$Microbe),as.vector(test.4cff$Microbe)))
# maxabspc <- max(abs(rangepc))
# maxabsmc <- max(abs(rangemc))
# pdf(paste(Sys.getenv("SCRATCH"), "/AllTraject_plant_extra.pdf",sep=""),height=5,width=3)
# layout(matrix(1:2,ncol=1,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# plottraj(test.4c$Plant,maxpos=rangepc[2],maxneg=rangepc[1])
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("Frequency",side=2,line=2,adj=-0.6)
# 	abline(v=which(win4c$dM==min(win4c$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# plottraj(test.4cff$Plant,maxpos=rangepc[2],maxneg=rangepc[1])
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	abline(v=which(win4cff$dP==min(win4cff$dP)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# 	mtext("generations",side=1,line=2)
# dev.off()
# pdf(paste(Sys.getenv("SCRATCH"), "/AllTraject_micr_extra.pdf",sep=""),height=5,width=3)
# layout(matrix(1:2,ncol=1,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# plottraj(test.4c$Microbe, type="micr",maxpos=rangemc[2],maxneg=rangemc[1])
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("Frequency",side=2,line=2,adj=-0.6)
# 	abline(v=which(win4c$dM==min(win4c$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4cff$Microbe, type="micr",maxpos=rangemc[2],maxneg=rangemc[1])
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	abline(v=which(win4cff$dM==min(win4cff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# 	mtext("generations",side=1,line=2)
# dev.off()

# pdf(paste(Sys.getenv("SCRATCH"), "/AllTraject_micr_six.pdf",sep=""),height=5,width=7)
# # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/AllTraject_micr.pdf",width=8,height=5)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# plottraj(test.4b$Microbe, type="micr",maxpos=rangem[2],maxneg=rangem[1])
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# 	mtext("Frequency",side=2,line=2,adj=-0.6)
# 	abline(v=which(win4b$dP==min(win4b$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4bff$Microbe, type="micr",maxpos=rangem[2],maxneg=rangem[1])
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	abline(v=which(win4bff$dP==min(win4bff$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4a$Microbe, type="micr",maxpos=rangem[2],maxneg=rangem[1])
# 	mtext("Same optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# 	abline(v=which(win4a$dP==min(win4a$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4aff$Microbe, type="micr",maxpos=rangem[2],maxneg=rangem[1])
# 	mtext("generations",side=1,line=2)
# 	abline(v=which(win4aff$dP==min(win4aff$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# plottraj(test.4f$Microbe, type="micr",maxpos=rangem[2],maxneg=rangem[1])
# 	mtext("Different optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# 	abline(v=which(win4f$dM==min(win4f$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# plottraj(test.4fff$Microbe, type="micr",maxpos=rangem[2],maxneg=rangem[1])
# 	abline(v=which(win4fff$dM==min(win4fff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# dev.off()

# soj4bM <- sojT(test.4b$Micr,type="micr")
# soj4bffM <- sojT(test.4bff$Micr,type="micr")
# soj4aM <- sojT(test.4a$Micr,type="micr")
# soj4affM <- sojT(test.4aff$Micr,type="micr")
# soj4fM <- sojT(test.4f$Micr,type="micr")
# soj4fffM <- sojT(test.4fff$Micr,type="micr")
# 
# 
# pdf(paste(Sys.getenv("SCRATCH"), "/sojsegT_micr_six.pdf",sep=""),height=5,width=7)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# plot( unlist(soj4bM$effs)[unlist(soj4bM$effs)!=0 & is.na(unlist(soj4bM$soj))] ~
# 	unlist(soj4bM$finalfrq) [ unlist(soj4bM$effs)!=0 & is.na(unlist(soj4bM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4bM$origingen))[unlist(soj4bM$effs)!=0 & is.na(unlist(soj4bM$soj))]/gens,0,0),
# 		ylim=c(-0.6,0.4),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# 	mtext("effect size",side=2,line=2,adj=-0.5)
# plot( unlist(soj4bffM$effs)[unlist(soj4bffM$effs)!=0 & is.na(unlist(soj4bffM$soj))] ~
# 	unlist(soj4bffM$finalfrq) [ unlist(soj4bffM$effs)!=0 & is.na(unlist(soj4bffM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4bffM$origingen))[unlist(soj4bffM$effs)!=0 & is.na(unlist(soj4bffM$soj))]/gens,0,0),
# 		ylim=c(-0.6,0.4),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# plot(unlist(soj4aM$effs)[unlist(soj4aM$effs)!=0 & is.na(unlist(soj4aM$soj))] ~
# 	 unlist(soj4aM$finalfrq) [ unlist(soj4aM$effs)!=0 & is.na(unlist(soj4aM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4aM$origingen))[unlist(soj4aM$effs)!=0 & is.na(unlist(soj4aM$soj))]/gens,0,0),
# 		ylim=c(-0.6,0.4),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("Same optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# plot( unlist(soj4affM$effs)[unlist(soj4affM$effs)!=0 & is.na(unlist(soj4affM$soj))] ~
# 	unlist(soj4affM$finalfrq) [ unlist(soj4affM$effs)!=0 & is.na(unlist(soj4affM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4affM$origingen))[unlist(soj4affM$effs)!=0 & is.na(unlist(soj4affM$soj))]/gens,0,0),
# 		ylim=c(-0.6,0.4),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
#  	mtext("derived allele frequency",side=1,line=2)
# plot( unlist(soj4fM$effs)[unlist(soj4fM$effs)!=0 & is.na(unlist(soj4fM$soj))]~
# 	unlist(soj4fM$finalfrq) [ unlist(soj4fM$effs)!=0 & is.na(unlist(soj4fM$soj))],  
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4fM$origingen))[unlist(soj4fM$effs)!=0 & is.na(unlist(soj4fM$soj))]/gens,0,0),
# 		ylim=c(-0.6,0.4),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("Different optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# plot(unlist(soj4fffM$effs)[unlist(soj4fffM$effs)!=0 & is.na(unlist(soj4fffM$soj))]~
# 	unlist(soj4fffM$finalfrq) [ unlist(soj4fffM$effs)!=0 & is.na(unlist(soj4fffM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4fffM$origingen))[unlist(soj4fffM$effs)!=0 & is.na(unlist(soj4fffM$soj))]/gens,0,0),
# 		ylim=c(-0.6,0.4),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# dev.off() 


# soj4bP <- sojT(test.4b$Plant)
# soj4bffP <- sojT(test.4bff$Plant)
# soj4aP <- sojT(test.4a$Plant)
# soj4affP <- sojT(test.4aff$Plant)
# soj4fP <- sojT(test.4f$Plant)
# soj4fffP <- sojT(test.4fff$Plant)
# 
# 
# 
# pdf(paste(Sys.getenv("SCRATCH"), "/sojsegT_plant_six.pdf",sep=""),height=5,width=7)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# plot( unlist(soj4bP$effs)[unlist(soj4bP$effs)!=0 & is.na(unlist(soj4bP$soj))] ~
# 	unlist(soj4bP$finalfrq) [ unlist(soj4bP$effs)!=0 & is.na(unlist(soj4bP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4bP$origingen))[unlist(soj4bP$effs)!=0 & is.na(unlist(soj4bP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# 	mtext("effect size",side=2,line=2,adj=-0.5)
# plot( unlist(soj4bffP$effs)[unlist(soj4bffP$effs)!=0 & is.na(unlist(soj4bffP$soj))] ~
# 	unlist(soj4bffP$finalfrq) [ unlist(soj4bffP$effs)!=0 & is.na(unlist(soj4bffP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4bffP$origingen))[unlist(soj4bffP$effs)!=0 & is.na(unlist(soj4bffP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# plot(unlist(soj4aP$effs)[unlist(soj4aP$effs)!=0 & is.na(unlist(soj4aP$soj))] ~
# 	 unlist(soj4aP$finalfrq) [ unlist(soj4aP$effs)!=0 & is.na(unlist(soj4aP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4aP$origingen))[unlist(soj4aP$effs)!=0 & is.na(unlist(soj4aP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("Same optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# plot( unlist(soj4affP$effs)[unlist(soj4affP$effs)!=0 & is.na(unlist(soj4affP$soj))] ~
# 	unlist(soj4affP$finalfrq) [ unlist(soj4affP$effs)!=0 & is.na(unlist(soj4affP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4affP$origingen))[unlist(soj4affP$effs)!=0 & is.na(unlist(soj4affP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
#  	mtext("Derived allele frequency",side=1,line=2)
# plot( unlist(soj4fP$effs)[unlist(soj4fP$effs)!=0 & is.na(unlist(soj4fP$soj))]~
# 	unlist(soj4fP$finalfrq) [ unlist(soj4fP$effs)!=0 & is.na(unlist(soj4fP$soj))],  
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4fP$origingen))[unlist(soj4fP$effs)!=0 & is.na(unlist(soj4fP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("Different optima",side=3,line=2)
# 	mtext("Equal links to fitness",side=3,line=0.5)
# plot(unlist(soj4fffP$effs)[unlist(soj4fffP$effs)!=0 & is.na(unlist(soj4fffP$soj))]~
# 	unlist(soj4fffP$finalfrq) [ unlist(soj4fffP$effs)!=0 & is.na(unlist(soj4fffP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4fffP$origingen))[unlist(soj4fffP$effs)!=0 & is.na(unlist(soj4fffP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# dev.off()  
# 
#EXTRA FIGURE
# soj4cM <- sojT(test.4c$Micr,type="micr")
# soj4cffM <- sojT(test.4cff$Micr,type="micr")
# soj4cP <- sojT(test.4c$Plant)
# soj4cffP <- sojT(test.4cff$Plant)
# pdf(paste(Sys.getenv("SCRATCH"), "/sojsegT_micr_extra.pdf",sep=""),height=5,width=3)
# layout(matrix(1:2,ncol=1,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# plot( unlist(soj4cM$effs)[unlist(soj4cM$effs)!=0 & is.na(unlist(soj4cM$soj))] ~
# 	unlist(soj4cM$finalfrq) [ unlist(soj4cM$effs)!=0 & is.na(unlist(soj4cM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4cM$origingen))[unlist(soj4cM$effs)!=0 & is.na(unlist(soj4cM$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("effect size",side=2,line=2,adj=-0.5)
# plot( unlist(soj4cffM$effs)[unlist(soj4cffM$effs)!=0 & is.na(unlist(soj4cffM$soj))] ~
# 	unlist(soj4cffM$finalfrq) [ unlist(soj4cffM$effs)!=0 & is.na(unlist(soj4cffM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4cffM$origingen))[unlist(soj4cffM$effs)!=0 & is.na(unlist(soj4cffM$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
#  	mtext("Derived allele frequency",side=1,line=2)
# dev.off()
#  pdf(paste(Sys.getenv("SCRATCH"), "/sojsegT_plant_extra.pdf",sep=""),height=5,width=3)
# layout(matrix(1:2,ncol=1,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# plot( unlist(soj4cP$effs)[unlist(soj4cP$effs)!=0 & is.na(unlist(soj4cP$soj))] ~
# 	unlist(soj4cP$finalfrq) [ unlist(soj4cP$effs)!=0 & is.na(unlist(soj4cP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4cP$origingen))[unlist(soj4cP$effs)!=0 & is.na(unlist(soj4cP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("effect size",side=2,line=2,adj=-0.5)
# plot( unlist(soj4cffP$effs)[unlist(soj4cffP$effs)!=0 & is.na(unlist(soj4cffP$soj))] ~
# 	unlist(soj4cffP$finalfrq) [ unlist(soj4cffP$effs)!=0 & is.na(unlist(soj4cffP$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4cffP$origingen))[unlist(soj4cffP$effs)!=0 & is.na(unlist(soj4cffP$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
#  	mtext("Derived allele frequency",side=1,line=2)
# dev.off()

# 
#   
# #fitness conflicts through time
# genvec <- c(25,75,300)
# traitfit.4b <- 	 lapply(genvec,function(gen)  getrelfitandtrait(test.4b  ,  gen+1,zoP=3,zoM=3,wP=0.75,wM=10,pfP=1,pfM=1))
# traitfit.4bff <- lapply(genvec,function(gen)  getrelfitandtrait(test.4bff  ,gen+1,zoP=3,zoM=3,wP=0.75,wM=10,pfP=0.6,pfM=0.6))
#   
# traitfit.4a <- 	 lapply(genvec,function(gen)  getrelfitandtrait(test.4a    ,gen+1,zoP=3,zoM=3,wP=0.75,wM=0.75,pfP=1,pfM=1))
# traitfit.4aff <- lapply(genvec,function(gen)  getrelfitandtrait(test.4aff  ,gen+1,zoP=3,zoM=3,wP=0.75,wM=0.75,pfP=0.6,pfM=0.6))
# 
# traitfit.4f <- 	 lapply(genvec,function(gen)  getrelfitandtrait(test.4f    ,gen+1,zoP=3,zoM=2,wP=0.75,wM=0.75,pfP=1,pfM=1))
# traitfit.4fff <- lapply(genvec,function(gen)  getrelfitandtrait(test.4fff  ,gen+1,zoP=3,zoM=2,wP=0.75,wM=0.75,pfP=0.6,pfM=0.6))
# 
# 
# pdf(paste(Sys.getenv("SCRATCH"),"/TraitFit_time_sims_sixdemos_n.pdf",sep=""),height=6,width=7)
# par(mfrow=c(3,3))
# par(mar=c(1,1,1,1))
# par(oma=c(4,5,2,1))
# for(i in 1:3){
#  xlims <- if(i==1){c(c(0,3.5))} else {c(1.4,3.6)}
#  plot( traitfit.4b[[i]]$rfitplnt~traitfit.4b[[i]]$traitexp,col=rgb(0,0.5,0,alpha=0.25), pch=16,xlim=xlims,ylim=c(0,0.0023),xaxt="n",cex=3)
# 	points( traitfit.4b[[i]]$rfitmicr~traitfit.4b[[i]]$traitexp,col=rgb(0.5,0,0.5,alpha=0.25),pch=16 ,cex=2.25)
# 	points( traitfit.4bff[[i]]$rfitplnt~traitfit.4bff[[i]]$traitexp,bg=rgb(0,0.5,0,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1.75)
# 	points( traitfit.4bff[[i]]$rfitmicr~traitfit.4bff[[i]]$traitexp,bg=rgb(0.5,0,0.5,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1)
# 	if(i==1){mtext("Microbe link none/indirect",side=2,line=3.5)}
# 	mtext(paste("Generation ",genvec[[i]],sep=""),side=3,line=1.5)
# 	abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=5)
# }	
# legend(1.3,0.0024,c("Host","Microbe"),fill=c(rgb(0,0.5,0),rgb(0.5,0,0.5)),bty="n")
# legend(1.3,0.0018,c("Direct only","With +FF"),pch=c(16,21),col=c(rgb(0,0,0,alpha=0.75),rgb(0,0,0)),pt.bg=rgb(0.5,0.5,0.5),pt.cex=c(2.5,1),bty="n")
# for(i in 1:3){
#  xlims <- if(i==1){c(c(0,3.5))} else {c(1.4,3.6)}
#  plot( traitfit.4a[[i]]$rfitplnt~traitfit.4a[[i]]$traitexp,col=rgb(0,0.5,0,alpha=0.25), pch=16,xlim=xlims,ylim=c(0,0.0023),xaxt="n",cex=3)
# 	points( traitfit.4a[[i]]$rfitmicr~traitfit.4a[[i]]$traitexp,col=rgb(0.5,0,0.5,alpha=0.25),pch=16 ,cex=2.25)
# 	points( traitfit.4aff[[i]]$rfitplnt~traitfit.4aff[[i]]$traitexp,bg=rgb(0,0.5,0,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1.75)
# 	points( traitfit.4aff[[i]]$rfitmicr~traitfit.4aff[[i]]$traitexp,bg=rgb(0.5,0,0.5,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1)
# 	if(i==1){mtext("relative fitness",side=2,line=2)}
# 	if(i==1){mtext("Same optima",side=2,line=3.5)}
# 	abline(v=3,col=rgb(0,0,0,alpha=0.5),lwd=5)
# }
# for(i in 1:3){
#  xlims <- if(i==1){c(c(0,3.5))} else {c(1.4,3.6)}
#  plot( traitfit.4f[[i]]$rfitplnt~traitfit.4f[[i]]$traitexp,col=rgb(0,0.5,0,alpha=0.25), pch=16,xlim=xlims,ylim=c(0,0.0023),cex=3)
# 	points( traitfit.4f[[i]]$rfitmicr~traitfit.4f[[i]]$traitexp,col=rgb(0.5,0,0.5,alpha=0.25),pch=16 ,cex=2.25)
# 	points( traitfit.4fff[[i]]$rfitplnt~traitfit.4fff[[i]]$traitexp,bg=rgb(0,0.5,0,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1.75)
# 	points( traitfit.4fff[[i]]$rfitmicr~traitfit.4fff[[i]]$traitexp,bg=rgb(0.5,0,0.5,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1)
# 	if(i==1){mtext("Different optima",side=2,line=3.5)}
# 	if(i==2){mtext("expressed trait value",side=1,line=2)}
# 	abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=5)
# 	abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=5)
# }
#  dev.off()
# 
# pdf(paste(Sys.getenv("SCRATCH"),"/Fit_corr_time_sims_sixdemos_n.pdf",sep=""),height=6,width=7)
# par(mfrow=c(3,3))
# par(mar=c(1,1,1,1))
# par(oma=c(4,5,2,1))
#  for(i in 1:3){
#  plot( traitfit.4b[[i]]$rfitplnt~traitfit.4b[[i]]$rfitmicr,col=rgb(0,0,0), pch=1,xlim=c(0,0.0023),ylim=c(0,0.0023),cex=2.5)
# 	points( traitfit.4bff[[i]]$rfitplnt~traitfit.4bff[[i]]$rfitmicr,col=rgb(0.5,0.5,0.5),pch=1 ,cex=1.5)
# 	if(i==1){mtext("Microbe link none/indirect",side=2,line=3.5)}
# # 	rhoD <- round(cor(traitfit.4b[[i]]$rfitplnt,traitfit.4b[[i]]$rfitmicr),digits=2)
# 	rhoI <- round(cor(traitfit.4bff[[i]]$rfitplnt,traitfit.4bff[[i]]$rfitmicr),digits=2)
# # 	text(0.0018,y=0.0005, bquote(rho==.(rhoD) ))
# 		#because there's no error here, this is actually perfectly correlated! the link is weaker than the error, but it can't be infinite, so...
# 	text(0.0018,y=0.0002, bquote(rho==.(rhoI) ),col=rgb(0.4,0.4,0.4))
# 	mtext(paste("Generation ",genvec[[i]],sep=""),side=3,line=1.5)
#  }
#  for(i in 1:3){
#  plot( traitfit.4a[[i]]$rfitplnt~traitfit.4a[[i]]$rfitmicr,col=rgb(0,0,0), pch=1,xlim=c(0,0.0023),ylim=c(0,0.0023),cex=2.5)
# 	points( traitfit.4aff[[i]]$rfitplnt~traitfit.4aff[[i]]$rfitmicr,col=rgb(0.5,0.5,0.5),pch=1 ,cex=1.5)
# 	if(i==1){mtext("relative fitness of host",side=2,line=2)}
# 	if(i==1){mtext("Same optima",side=2,line=3.5)}
# 	rhoD <- round(cor(traitfit.4a[[i]]$rfitplnt,traitfit.4a[[i]]$rfitmicr),digits=2)
# 	rhoI <- round(cor(traitfit.4aff[[i]]$rfitplnt,traitfit.4aff[[i]]$rfitmicr),digits=2)
# 	text(0.0018,y=0.0005, bquote(rho==.(rhoD) ))
# 	text(0.0018,y=0.0002, bquote(rho==.(rhoI) ),col=rgb(0.4,0.4,0.4))
#  }
#  for(i in 1:3){
#  plot( traitfit.4f[[i]]$rfitplnt~traitfit.4f[[i]]$rfitmicr,col=rgb(0,0,0), pch=1,xlim=c(0,0.0023),ylim=c(0,0.0023),cex=2.5)
# 	points( traitfit.4fff[[i]]$rfitplnt~traitfit.4fff[[i]]$rfitmicr,col=rgb(0.5,0.5,0.5),pch=1 ,cex=1.5)
# 	if(i==1){mtext("Different optima",side=2,line=3.5)}
# 	if(i==2){mtext("relative fitness of microbe",side=1,line=2)}
# 	rhoD <- round(cor(traitfit.4f[[i]]$rfitplnt,traitfit.4f[[i]]$rfitmicr),digits=2)
# 	rhoI <- round(cor(traitfit.4fff[[i]]$rfitplnt,traitfit.4fff[[i]]$rfitmicr),digits=2)
# 	text(0.0018,y=0.0005, bquote(rho==.(rhoD) ))
# 	text(0.0018,y=0.0002, bquote(rho==.(rhoI) ),col=rgb(0.4,0.4,0.4))
# # 	legend(0,0.0018,c("Direct only","With +FF"),pch=c(16,21),col=rgb(0,0,0),pt.bg=rgb(0.5,0.5,0.5),pt.cex=c(2.5,1),bty="n")
# 	if(i==3){legend(0,0.0024,c("Direct only","With +FF"),pch=1,pt.cex=c(2.5,1.5),col=c(rgb(0,0,0),rgb(0.5,0.5,0.5)) ,bty="n")}
#  }
# dev.off()
# 

#EXTRA FIGURE
# traitfit.4c <- 	 lapply(genvec,function(gen)  getrelfitandtrait(test.4c  ,  gen+1,zoP=3,zoM=3,wP=10,wM=0.75,pfP=1,pfM=1))
# traitfit.4cff <- lapply(genvec,function(gen)  getrelfitandtrait(test.4cff  ,gen+1,zoP=3,zoM=3,wP=10,wM=0.75,pfP=0.6,pfM=0.6))
# pdf(paste(Sys.getenv("SCRATCH"),"/TraitFit_time_sims_extra_n.pdf",sep=""),height=3,width=7)
# par(mfrow=c(1,3))
# par(mar=c(1,1,1,1))
# par(oma=c(4,5,2,1))
# for(i in 1:3){
#  xlims <- if(i==1){c(c(0,3.5))} else {c(1.4,3.75)}
#  plot( traitfit.4c[[i]]$rfitplnt~traitfit.4c[[i]]$traitexp,col=rgb(0,0.5,0,alpha=0.25), pch=16,xlim=xlims,ylim=c(0,0.0030),cex=3)
# 	points( traitfit.4c[[i]]$rfitmicr~traitfit.4c[[i]]$traitexp,col=rgb(0.5,0,0.5,alpha=0.25),pch=16 ,cex=2.25)
# 	points( traitfit.4cff[[i]]$rfitplnt~traitfit.4cff[[i]]$traitexp,bg=rgb(0,0.5,0,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1.75)
# 	points( traitfit.4cff[[i]]$rfitmicr~traitfit.4cff[[i]]$traitexp,bg=rgb(0.5,0,0.5,alpha=1),col=rgb(0,0,0),pch=21 ,cex=1)
# 	if(i==1){mtext("relative fitness",side=2,line=2)}
# 	if(i==1){mtext("Host link none/indirect",side=2,line=3.5)}
# 	if(i==2){mtext("expressed trait value",side=1,line=2)}
# 	mtext(paste("Generation ",genvec[[i]],sep=""),side=3,line=1.5)
# 	abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=5)
# }	
# dev.off()
# pdf(paste(Sys.getenv("SCRATCH"),"/Fit_corr_time_sims_extra_n.pdf",sep=""),height=3,width=7)
# par(mfrow=c(1,3))
# par(mar=c(1,1,1,1))
# par(oma=c(4,5,2,1))
#  for(i in 1:3){
#  plot( traitfit.4c[[i]]$rfitplnt~traitfit.4c[[i]]$rfitmicr,col=rgb(0,0,0), pch=1,xlim=c(0,0.0030),ylim=c(0,0.0030),cex=2.5)
# 	points( traitfit.4cff[[i]]$rfitplnt~traitfit.4cff[[i]]$rfitmicr,col=rgb(0.5,0.5,0.5),pch=1 ,cex=1.5)
# 	if(i==1){mtext("Host link none/indirect",side=2,line=3.5)}
# 	rhoI <- round(cor(traitfit.4cff[[i]]$rfitplnt,traitfit.4cff[[i]]$rfitmicr),digits=2)
# # 	text(0.0018,y=0.0005, bquote(rho==.(rhoD) ))
# 		#because there's no error here, this is actually perfectly correlated! the link is weaker than the error, but it can't be infinite, so...
# 	text(0.0018,y=0.0002, bquote(rho==.(rhoI) ),col=rgb(0.4,0.4,0.4))
# 	if(i==1){mtext("relative fitness of host",side=2,line=2)}
# 	if(i==2){mtext("relative fitness of microbe",side=1,line=2)}
# 	mtext(paste("Generation ",genvec[[i]],sep=""),side=3,line=1.5)
#  }
# dev.off()
#
# 
# 
# 
##EXTRA fig combo of main text figs
# pdf(paste(Sys.getenv("SCRATCH"), "/Simulation_Results_COMBOextrademos_n.pdf",sep=""),height=5,width=7)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,1,0))
# windowplot(1,gens+1, 3,test.4c,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[3])
# 	abline(h=3,col=rgb(0.5,0,0.5),lty=2,lwd=2)
# 	abline(h=0)
# 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("breeding values",side=2,line=2,adj=-1)
# 	legend(gens*-0.05,y=0,c("Microbe optimum"), lty = 2,lwd=2, col = rgb(0.5,0,0.5), bty="n")
# windowplot(1,gens+1, 3,test.4cff,ylim=c(-3,6),"",xlabs="",ylabs="")#"titles[4]")
# 	abline(h=3,col=rgb(0.5,0,0.5),lty=2,lwd=2)
# 	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	legend(gens*-0.05,y=0,c("Host","Microbe","Expressed"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
# 	mtext("generations",side=1,line=2)
# plot(VmVp4c$Vp ~ c(1:length(VmVp4c$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 	lines( 1:length(VmVp4c$Vp),VmVp4c$Vp,col=rgb(0,0.5,0)) 
# 	lines(1:length(VmVp4c$Vm),VmVp4c$Vm , col=rgb(0.5,0,0.5)) 
# 	lines(1:length(VmVp4c$Vb),VmVp4c$Vb , col=rgb(0.5,0.5,0.5)) 
# 	abline(v=which(win4c$dM==min(win4c$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# #  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("genetic variance",side=2,line=2,adj=-1.25)
# plot(VmVp4cff$Vp ~ c(1:length(VmVp4cff$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
# 	lines( 1:length(VmVp4cff$Vp),VmVp4cff$Vp,col=rgb(0,0.5,0)) 
# 	lines(1:length(VmVp4cff$Vm),VmVp4cff$Vm , col=rgb(0.5,0,0.5)) 
# 	lines(1:length(VmVp4cff$Vb),VmVp4cff$Vb , col=rgb(0.5,0.5,0.5)) 
# 	abline(v=which(win4cff$dM==min(win4cff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# # 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("generations",side=1,line=2)
# 	legend(gens*0.4,y=0.10,c(expression(V[H]),expression(V[M]),expression(V[A])),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
# plot( unlist(soj4cM$effs)[unlist(soj4cM$effs)!=0 & is.na(unlist(soj4cM$soj))] ~
# 	unlist(soj4cM$finalfrq) [ unlist(soj4cM$effs)!=0 & is.na(unlist(soj4cM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4cM$origingen))[unlist(soj4cM$effs)!=0 & is.na(unlist(soj4cM$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# #  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to host fitness",side=3,line=0.5)
# 	mtext("effect size",side=2,line=2,adj=-0.5)
# plot( unlist(soj4cffM$effs)[unlist(soj4cffM$effs)!=0 & is.na(unlist(soj4cffM$soj))] ~
# 	unlist(soj4cffM$finalfrq) [ unlist(soj4cffM$effs)!=0 & is.na(unlist(soj4cffM$soj))],
# 	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
# 	col = rgb( (unlist(soj4cffM$origingen))[unlist(soj4cffM$effs)!=0 & is.na(unlist(soj4cffM$soj))]/gens,0,0),
# 		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
# 	abline(h=0)
# # 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
#  	mtext("Derived allele frequency",side=1,line=2)
# dev.off()

###EXTRA COMBO FIG 2, reverse host and microbe in conflict scenario
VmVp4e <- extractVmVp(test.4e, 1,gens+1,1)
VmVp4eff <- extractVmVp(test.4eff, 1,gens+1,1)
win4e  <- extractwinning(test.4e,first=1,last=gens+1,1,zoP=3,zoM= 3)
win4eff  <- extractwinning(test.4eff,first=1,last=gens+1,1,zoP=3,zoM= 3)
soj4eM <- sojT(test.4e$Micr,type="micr")
soj4effM <- sojT(test.4eff$Micr,type="micr")
soj4eP <- sojT(test.4e$Plant)
soj4effP <- sojT(test.4eff$Plant)
pdf(paste(Sys.getenv("SCRATCH"), "/Simulation_Results_COMBOextrademosConfl_n.pdf",sep=""),height=5,width=7)
layout(matrix(1:6,ncol=3,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,1,0))
windowplot(1,gens+1, 3,test.4e,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[3])
	abline(h=2,col=rgb(0,0.5,0),lty=2,lwd=3); abline(h=3,col=rgb(0.5,0,0.5),lty=2,lwd=3); abline(h=2.5,col=rgb(0,0,0),lty=2,lwd=3) 
	abline(h=0)
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
	mtext("breeding values",side=2,line=2,adj=-1)
	legend(gens*-0.05,y=0,c("Host optimum","Microbe optimum","Midpoint"), lty = c(2,2,3),lwd=2, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0)), bty="n")
windowplot(1,gens+1, 3,test.4eff,ylim=c(-3,6),"",xlabs="",ylabs="")#"titles[4]")
	abline(h=2,col=rgb(0,0.5,0),lty=2,lwd=3); abline(h=3,col=rgb(0.5,0,0.5),lty=2,lwd=3); abline(h=2.5,col=rgb(0,0,0),lty=2,lwd=3) 
	abline(h=0)
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	legend(gens*-0.05,y=0,c("Host","Microbe","Expressed"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
	mtext("generations",side=1,line=2)
plot(VmVp4e$Vp ~ c(1:length(VmVp4e$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
	lines( 1:length(VmVp4e$Vp),VmVp4e$Vp,col=rgb(0,0.5,0)) 
	lines(1:length(VmVp4e$Vm),VmVp4e$Vm , col=rgb(0.5,0,0.5)) 
	lines(1:length(VmVp4e$Vb),VmVp4e$Vb , col=rgb(0.5,0.5,0.5)) 
	abline(v=which(win4e$dM==min(win4e$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
	mtext("genetic variance",side=2,line=2,adj=-1.25)
plot(VmVp4eff$Vp ~ c(1:length(VmVp4eff$Vp)),pch=NA,ylim=c(0,0.1),ylab="",xlab="")
	lines( 1:length(VmVp4eff$Vp),VmVp4eff$Vp,col=rgb(0,0.5,0)) 
	lines(1:length(VmVp4eff$Vm),VmVp4eff$Vm , col=rgb(0.5,0,0.5)) 
	lines(1:length(VmVp4eff$Vb),VmVp4eff$Vb , col=rgb(0.5,0.5,0.5)) 
	abline(v=which(win4eff$dM==min(win4eff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #first generation does microbe mean reach its optima or beyond
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("generations",side=1,line=2)
	legend(gens*0.4,y=0.10,c(expression(V[H]),expression(V[M]),expression(V[A])),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
plot( unlist(soj4eM$effs)[unlist(soj4eM$effs)!=0 & is.na(unlist(soj4eM$soj))] ~
	unlist(soj4eM$finalfrq) [ unlist(soj4eM$effs)!=0 & is.na(unlist(soj4eM$soj))],
	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
	col = rgb( (unlist(soj4eM$origingen))[unlist(soj4eM$effs)!=0 & is.na(unlist(soj4eM$soj))]/gens,0,0),
		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
	abline(h=0)
#  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 	mtext("No direct link to host fitness",side=3,line=0.5)
	mtext("effect size",side=2,line=2,adj=-0.5)
plot( unlist(soj4effM$effs)[unlist(soj4effM$effs)!=0 & is.na(unlist(soj4effM$soj))] ~
	unlist(soj4effM$finalfrq) [ unlist(soj4effM$effs)!=0 & is.na(unlist(soj4effM$soj))],
	  pch=1, #col=rgb(0,0.5,0,alpha=0.5),
	col = rgb( (unlist(soj4effM$origingen))[unlist(soj4effM$effs)!=0 & is.na(unlist(soj4effM$soj))]/gens,0,0),
		ylim=c(-0.3,0.7),xlim=c(0,1), ylab="",xlab="")
	abline(h=0)
# 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
 	mtext("Derived allele frequency",side=1,line=2)
dev.off()

# 
# 
# ##this figure not included in manuscript, but may be of interest.
# 
# pdf(paste(Sys.getenv("SCRATCH"), "/coefV_Results_sixdemos_n.pdf",sep=""),height=5,width=7)
# # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/coefV_Results_fourdemos_n.pdf",height=4.5,width=5.5)
# layout(matrix(1:6,ncol=3,byrow=F))
# par(mar=c(1.5,3,1,1))
# par(oma=c(2,2,3,0))
# 	plot(VmVp4b$cVp ~ c(1:length(VmVp4b$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
# 		lines( 1:length(VmVp4b$cVp),VmVp4b$cVp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4b$cVm),VmVp4b$cVm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4b$cVb),VmVp4b$cVb , col=rgb(0.5,0.5,0.5)) 
# 	 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# 		mtext("No direct link to microbe fitness",side=3,line=0.5)
# 		mtext("genetic variance scaled by expressed mean",side=2,line=2,adj=1.08)
# 		abline(v=which(win4b$dP==min(win4b$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# 	plot(VmVp4bff$cVp ~ c(1:length(VmVp4bff$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
# 		lines( 1:length(VmVp4bff$cVp),VmVp4bff$cVp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4bff$cVm),VmVp4bff$cVm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4bff$cVb),VmVp4bff$cVb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4bff$dP==min(win4bff$dP)),col=rgb(0,0.5,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# 		mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# 	plot(VmVp4a$cVp ~ c(1:length(VmVp4a$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
# 		lines( 1:length(VmVp4a$cVp),VmVp4a$cVp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4a$cVm),VmVp4a$cVm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4a$cVb),VmVp4a$cVb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4a$dP==min(win4a$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# 		mtext("Same optima",side=3,line=2)
# 		mtext("Equal links to fitness",side=3,line=0.5)
# 	plot(VmVp4aff$cVp ~ c(1:length(VmVp4aff$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
# 		lines( 1:length(VmVp4aff$cVp),VmVp4aff$cVp,col=rgb(0,0.5,0)) 
# 		lines(1:length(VmVp4aff$cVm),VmVp4aff$cVm , col=rgb(0.5,0,0.5)) 
# 		lines(1:length(VmVp4aff$cVb),VmVp4aff$cVb , col=rgb(0.5,0.5,0.5)) 
# 		abline(v=which(win4aff$dP==min(win4aff$dP)),col=rgb(0,0,0,alpha=0.5),lwd=5) #first generation does plant mean reach its optima or beyond
# 		mtext("generations",side=1,line=2)
# 	plot(VmVp4f$cVp ~ c(1:length(VmVp4f$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
# 		lines( 1:length(VmVp4f$cVp),VmVp4f$cVp,col=rgb(0,0.5,0))
# 		lines(1:length(VmVp4f$cVm),VmVp4f$cVm , col=rgb(0.5,0,0.5))
# 		lines(1:length(VmVp4f$cVb),VmVp4f$cVb , col=rgb(0.5,0.5,0.5))
# 		abline(v=which(win4f$dM==min(win4f$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# # 		text(x=85,y=0.12, labels =expression(Z[opt[M]]~reached))
# 		mtext("Different optima",side=3,line=2)
# 		mtext("Equal links to fitness",side=3,line=0.5)
# 	plot(VmVp4fff$cVp ~ c(1:length(VmVp4fff$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
# 		lines( 1:length(VmVp4fff$cVp),VmVp4fff$cVp,col=rgb(0,0.5,0))
# 		lines(1:length(VmVp4fff$cVm),VmVp4fff$cVm , col=rgb(0.5,0,0.5))
# 		lines(1:length(VmVp4fff$cVb),VmVp4fff$cVb , col=rgb(0.5,0.5,0.5))
# 		abline(v=which(win4fff$dM==min(win4fff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
# 		legend(gens*0.2,y=0.16,c("Host portion","Microbe portion","Interacting"),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
# dev.off()
# 
# 
# 
# # ## not needed for current figures
# # #characterize segregating variants
# # finalplant4f<- test.4f$Plant[,,,gens+1]
# # finalplant4fff<- test.4fff$Plant[,,,gens+1]
# # finalplant4a<- test.4a$Plant[,,,gens+1]
# # finalplant4aff<- test.4aff$Plant[,,,gens+1]
# # finalplant4b<- test.4b$Plant[,,,gens+1]
# # finalplant4bff<- test.4bff$Plant[,,,gens+1]
# # 
# # # finalmicr<- test.4fff$Microbe[,,gens+1]
# # finalmicr4f<- test.4f$Microbe[,,gens+1]
# # finalmicr4fff<- test.4fff$Microbe[,,gens+1]
# # finalmicr4a<- test.4a$Microbe[,,gens+1]
# # finalmicr4aff<- test.4aff$Microbe[,,gens+1]
# # finalmicr4b<- test.4b$Microbe[,,gens+1]
# # finalmicr4bff<- test.4bff$Microbe[,,gens+1]
# # #
# # unlistallP4fff <- getsegcounts(finalplant4fff,type="plant")
# # unlistallP4f <- getsegcounts(finalplant4f,type="plant")
# # unlistallP4aff <- getsegcounts(finalplant4aff,type="plant")
# # unlistallP4a <- getsegcounts(finalplant4a,type="plant")
# # unlistallP4bff <- getsegcounts(finalplant4bff,type="plant")
# # unlistallP4b <- getsegcounts(finalplant4b,type="plant")
# # unlistallM4fff <- getsegcounts(finalmicr4fff,type="microbe")
# # unlistallM4f <- getsegcounts(finalmicr4f,type="microbe")
# # unlistallM4aff <- getsegcounts(finalmicr4aff,type="microbe")
# # unlistallM4a <- getsegcounts(finalmicr4a,type="microbe")
# # unlistallM4bff <- getsegcounts(finalmicr4bff,type="microbe")
# # unlistallM4b <- getsegcounts(finalmicr4b,type="microbe")
# # 
# ####these figures not currently used, though may be of interest.
# # pdf(paste(Sys.getenv("SCRATCH"), "/PooledLociAndAlleles_effectsizedistr_six.pdf",sep=""),height=5,width=7)
# # # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PooledLociAndAlleles_effectsizedistr.pdf",height=4.5,width=5.5)
# # layout(matrix(1:6,ncol=3,byrow=F))
# # par(mar=c(1.5,3,1,1))
# # par(oma=c(2,2,3,0))
# # hist(abs(as.vector(unlistallP4b))[as.vector(unlistallP4b)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,25000))
# # 	hist(abs(as.vector(unlistallM4b))[as.vector(unlistallM4b)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# #  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# # 	mtext("Frequency",side=2,line=2,adj=-0.75)
# # hist(abs(as.vector(unlistallP4bff))[as.vector(unlistallP4bff)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,25000))
# # 	hist(abs(as.vector(unlistallM4bff))[as.vector(unlistallM4bff)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# # hist(abs(as.vector(unlistallP4a))[as.vector(unlistallP4a)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,25000))
# # 	hist(abs(as.vector(unlistallM4a))[as.vector(unlistallM4a)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 		mtext("Same optima",side=3,line=2)
# # 		mtext("Equal links to fitness",side=3,line=0.5)
# # hist(abs(as.vector(unlistallP4aff))[as.vector(unlistallP4aff)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,25000))
# # 	hist(abs(as.vector(unlistallM4aff))[as.vector(unlistallM4aff)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 	mtext("Effect Size",side=1,line=2)
# # hist(abs(as.vector(unlistallP4f))[as.vector(unlistallP4f)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,25000))
# # 	hist(abs(as.vector(unlistallM4f))[as.vector(unlistallM4f)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 		mtext("Different optima",side=3,line=2)
# # 		mtext("Equal links to fitness",side=3,line=0.5)
# # hist(abs(as.vector(unlistallP4fff))[as.vector(unlistallP4fff)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,25000))
# # 	hist(abs(as.vector(unlistallM4fff))[as.vector(unlistallM4fff)!=0],breaks=seq(from=0,to=1,by=0.05),
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # dev.off()
# #  #includes all frequencies -- e.g biallelic in 2x, non biallelic >2x, includes alleles that are still the 0 state.
# # 
# # pdf(paste(Sys.getenv("SCRATCH"), "/PooledLociAndAlleles_effectsizedistr_sixDens.pdf",sep=""),height=5,width=7)
# # # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PooledLociAndAlleles_effectsizedistr.pdf",height=4.5,width=5.5)
# # layout(matrix(1:6,ncol=3,byrow=F))
# # par(mar=c(1.5,3,1,1))
# # par(oma=c(2,2,3,0))
# # # hist(abs(as.vector(unlistallP4b))[as.vector(unlistallP4b)!=0],breaks=seq(from=0,to=1,by=0.05),
# # hist(abs(as.vector(unlistallP4b))[as.vector(unlistallP4b)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,20))
# # 	hist(abs(as.vector(unlistallM4b))[as.vector(unlistallM4b)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# #  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# # 	mtext("Probability Density",side=2,line=2,adj=-0.75)
# # hist(abs(as.vector(unlistallP4bff))[as.vector(unlistallP4bff)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,20))
# # 	hist(abs(as.vector(unlistallM4bff))[as.vector(unlistallM4bff)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# # hist(abs(as.vector(unlistallP4a))[as.vector(unlistallP4a)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,20))
# # 	hist(abs(as.vector(unlistallM4a))[as.vector(unlistallM4a)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 		mtext("Same optima",side=3,line=2)
# # 		mtext("Equal links to fitness",side=3,line=0.5)
# # hist(abs(as.vector(unlistallP4aff))[as.vector(unlistallP4aff)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,20))
# # 	hist(abs(as.vector(unlistallM4aff))[as.vector(unlistallM4aff)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 	mtext("Effect Size",side=1,line=2)
# # hist(abs(as.vector(unlistallP4f))[as.vector(unlistallP4f)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,20))
# # 	hist(abs(as.vector(unlistallM4f))[as.vector(unlistallM4f)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # 		mtext("Different optima",side=3,line=2)
# # 		mtext("Equal links to fitness",side=3,line=0.5)
# # hist(abs(as.vector(unlistallP4fff))[as.vector(unlistallP4fff)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 	ylab="",xlab="",main="",col=rgb(0,0.5,0,alpha=0.75),ylim=c(0,20))
# # 	hist(abs(as.vector(unlistallM4fff))[as.vector(unlistallM4fff)!=0],breaks=seq(from=0,to=1,by=0.05),freq=F,
# # 		add=T,col=rgb(0.5,0,0.5,alpha=0.75))
# # dev.off()
# # 
# # 
# # 
# # 
# 
# ####these figures not currently used, though may be of interest.
# # pdf(paste(Sys.getenv("SCRATCH"), "/sojT_plant_six.pdf",sep=""),height=5,width=7)
# # # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sojT_plant.pdf",width=5,height=5)
# # layout(matrix(1:6,ncol=3,byrow=F))
# # par(mar=c(1.5,3,1,1))
# # par(oma=c(2,2,1,0))
# # plot((unlist(soj4bP$soj))[unlist(soj4bP$effs)!=0 & !is.na(unlist(soj4bP$soj))]~
# # 	  unlist(soj4bP$effs)[unlist(soj4bP$effs)!=0 & !is.na(unlist(soj4bP$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01( unlist(soj4bP$origingen) [ unlist(soj4bP$effs)!=0 & !is.na(unlist(soj4bP$soj))] ),0,0),# alpha = c(0.1,1)[(unlist(soj4bP$fstate)[unlist(soj4bP$effs)!=0] == "fixed") +1 ]), 
# # 	cex = c(0.1,1)[(unlist(soj4bP$fstate)[unlist(soj4bP$effs)!=0 & !is.na(unlist(soj4bP$soj))] == "fixed") +1 ])
# # 	# 	col=rgb(0,0.5,0), cex = c(0.1,1)[(unlist(soj4bP$fstate)[unlist(soj4bP$effs)!=0] == "fixed") +1 ])
# # 	abline(v=0)
# #  	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# # 	mtext("Sojourn (gens)",side=2,line=2,adj=-1.08)
# # plot((unlist(soj4bffP$soj))[unlist(soj4bffP$effs)!=0 & !is.na(unlist(soj4bffP$soj))]~
# # 	  unlist(soj4bffP$effs)[unlist(soj4bffP$effs)!=0 & !is.na(unlist(soj4bffP$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # # 	col=rgb(0,0.5,0), cex = c(0.1,1)[(unlist(soj4bffP$fstate)[unlist(soj4bffP$effs)!=0] == "fixed") +1 ])
# # 	col=rgb(range01(unlist(soj4bffP$origingen)[unlist(soj4bffP$effs)!=0 & !is.na(unlist(soj4bffP$soj))]),0,0), 
# # 	cex = c(0.1,1)[(unlist(soj4bffP$fstate)[unlist(soj4bffP$effs)!=0 & !is.na(unlist(soj4bffP$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# # plot((unlist(soj4aP$soj))[unlist(soj4aP$effs)!=0 & !is.na(unlist(soj4aP$soj))]~
# # 	  unlist(soj4aP$effs)[unlist(soj4aP$effs)!=0 & !is.na(unlist(soj4aP$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01( unlist(soj4aP$origingen) [ unlist(soj4aP$effs)!=0 & !is.na(unlist(soj4aP$soj))] ),0,0),# alpha = c(0.1,1)[(unlist(soj4bP$fstate)[unlist(soj4bP$effs)!=0] == "fixed") +1 ]), 
# # 	cex = c(0.1,1)[(unlist(soj4aP$fstate)[unlist(soj4aP$effs)!=0 & !is.na(unlist(soj4aP$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("Different optima",side=3,line=3)
# # 	mtext("Equal links to fitness",side=3,line=0.5)
# # plot((unlist(soj4affP$soj))[unlist(soj4affP$effs)!=0 & !is.na(unlist(soj4affP$soj))]~
# # 	  unlist(soj4affP$effs)[unlist(soj4affP$effs)!=0 & !is.na(unlist(soj4affP$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4affP$origingen)[unlist(soj4affP$effs)!=0 & !is.na(unlist(soj4affP$soj))]),0,0), 
# # 	cex = c(0.1,1)[(unlist(soj4affP$fstate)[unlist(soj4affP$effs)!=0 & !is.na(unlist(soj4affP$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# #  	mtext("Effect size",side=1,line=2)
# # plot((unlist(soj4fP$soj))[unlist(soj4fP$effs)!=0 & !is.na(unlist(soj4fP$soj))]~
# # 	  unlist(soj4fP$effs)[unlist(soj4fP$effs)!=0 & !is.na(unlist(soj4fP$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4fP$origingen)[unlist(soj4fP$effs)!=0 & !is.na(unlist(soj4fP$soj))]),0,0), 
# # 	cex = c(0.1,1)[(unlist(soj4fP$fstate)[unlist(soj4fP$effs)!=0 & !is.na(unlist(soj4fP$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("Different optima",side=3,line=3)
# # 	mtext("Equal links to fitness",side=3,line=0.5)
# # plot((unlist(soj4fffP$soj))[unlist(soj4fffP$effs)!=0 & !is.na(unlist(soj4fffP$soj))]~
# # 	  unlist(soj4fffP$effs)[unlist(soj4fffP$effs)!=0 & !is.na(unlist(soj4fffP$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4fffP$origingen)[(unlist(soj4fffP$effs)!=0 & !is.na(unlist(soj4fffP$soj)))]),0,0), 
# # 	cex = c(0.1,1)[(unlist(soj4fffP$fstate)[unlist(soj4fffP$effs)!=0 & !is.na(unlist(soj4fffP$soj))] == "fixed") +1 ])
# # 	abline(v=0) 
# # dev.off()
# # 
# # pdf(paste(Sys.getenv("SCRATCH"), "/sojT_micr_six.pdf",sep=""),height=5,width=7)
# # # pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sojT_micr.pdf",width=5,height=5)
# # layout(matrix(1:6,ncol=3,byrow=F))
# # par(mar=c(1.5,3,1,1))
# # par(oma=c(2,2,1,0))
# # plot((unlist(soj4bM$soj))[unlist(soj4bM$effs)!=0 & !is.na(unlist(soj4bM$soj))]~
# # 	  unlist(soj4bM$effs)[unlist(soj4bM$effs)!=0 & !is.na(unlist(soj4bM$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4bM$origingen)[unlist(soj4bM$effs)!=0 & !is.na(unlist(soj4bM$soj))]),0,0),
# # 	cex = c(0.1,1)[(unlist(soj4bM$fstate)[unlist(soj4bM$effs)!=0 & !is.na(unlist(soj4bM$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
# # 	mtext("No direct link to microbe fitness",side=3,line=0.5)
# # 	mtext("Sojourn (gens)",side=2,line=2,adj=-1.08)
# # plot((unlist(soj4bffM$soj))[unlist(soj4bffM$effs)!=0 & !is.na(unlist(soj4bffM$soj))]~
# # 	  unlist(soj4bffM$effs)[unlist(soj4bffM$effs)!=0 & !is.na(unlist(soj4bffM$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4bffM$origingen)[unlist(soj4bffM$effs)!=0 & !is.na(unlist(soj4bffM$soj))]),0,0),
# # 	cex = c(0.1,1)[(unlist(soj4bffM$fstate)[unlist(soj4bffM$effs)!=0 & !is.na(unlist(soj4bffM$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
# # 
# # plot((unlist(soj4aM$soj))[unlist(soj4aM$effs)!=0 & !is.na(unlist(soj4aM$soj))]~
# # 	  unlist(soj4aM$effs)[unlist(soj4aM$effs)!=0 & !is.na(unlist(soj4aM$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4aM$origingen)[unlist(soj4aM$effs)!=0 & !is.na(unlist(soj4aM$soj))]),0,0),
# # 	cex = c(0.1,1)[(unlist(soj4aM$fstate)[unlist(soj4aM$effs)!=0 & !is.na(unlist(soj4aM$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("Same optima",side=3,line=3)
# # 	mtext("Equal links to fitness",side=3,line=0.5)
# # plot((unlist(soj4affM$soj))[unlist(soj4affM$effs)!=0 & !is.na(unlist(soj4affM$soj))]~
# # 	  unlist(soj4affM$effs)[unlist(soj4affM$effs)!=0 & !is.na(unlist(soj4affM$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4affM$origingen)[unlist(soj4affM$effs)!=0 & !is.na(unlist(soj4affM$soj))]),0,0),
# # 	cex = c(0.1,1)[(unlist(soj4affM$fstate)[unlist(soj4affM$effs)!=0 & !is.na(unlist(soj4affM$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# #  	mtext("Effect size",side=1,line=2)
# # plot((unlist(soj4fM$soj))[unlist(soj4fM$effs)!=0 & !is.na(unlist(soj4fM$soj))]~
# # 	  unlist(soj4fM$effs)[unlist(soj4fM$effs)!=0 & !is.na(unlist(soj4fM$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4fM$origingen)[unlist(soj4fM$effs)!=0 & !is.na(unlist(soj4fM$soj))]),0,0),
# # 	cex = c(0.1,1)[(unlist(soj4fM$fstate)[unlist(soj4fM$effs)!=0 & !is.na(unlist(soj4fM$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # 	mtext("Different optima",side=3,line=3)
# # 	mtext("Equal links to fitness",side=3,line=0.5)
# # plot((unlist(soj4fffM$soj))[unlist(soj4fffM$effs)!=0 & !is.na(unlist(soj4fffM$soj))]~
# # 	  unlist(soj4fffM$effs)[unlist(soj4fffM$effs)!=0 & !is.na(unlist(soj4fffM$soj))],
# # 	xlim=c(-0.7,0.7),ylim=c(0,gens), ylab="",xlab="",
# # 	col=rgb(range01(unlist(soj4fffM$origingen)[unlist(soj4fffM$effs)!=0 & !is.na(unlist(soj4fffM$soj))]),0,0),
# # 	cex = c(0.1,1)[(unlist(soj4fffM$fstate)[unlist(soj4fffM$effs)!=0 & !is.na(unlist(soj4fffM$soj))] == "fixed") +1 ])
# # 	abline(v=0)
# # dev.off()
# #  
#   
#   
# 
# 
# 
# ## possibly of interest, but not included:
# # #weaker variance correlations with fitness feedbacks. interesting. but the plots make it look like this could depend on timing of events, which change between sims/scenarios
# # cor(VmVp4fff$Vp,VmVp4fff$Vm) #[1] 0.2660599
# # cor(VmVp4fff$Vp,VmVp4fff$Vb,use="complete.obs") #[1] 0.8572949
# # cor(VmVp4fff$Vm,VmVp4fff$Vb,use="complete.obs") #[1] 0.6107672
# # cor(VmVp4fff$cVp,VmVp4fff$cVm,use="complete.obs") #[1] 0.9656252
# # cor(VmVp4fff$cVm,VmVp4fff$cVb,use="complete.obs") #[1] 0.9980781
# # cor(VmVp4fff$cVp,VmVp4fff$cVb,use="complete.obs")#[1] 0.9791138
# # cor(VmVp4f$Vp,VmVp4f$Vm) #[1] 0.3198037
# # cor(VmVp4f$Vp,VmVp4f$Vb) #[1] 0.9296124
# # cor(VmVp4f$Vm,VmVp4f$Vb) #[1] 0.5733526
# # cor(VmVp4f$cVp,VmVp4f$cVm,use="complete.obs")#[1] 0.998429
# # cor(VmVp4f$cVp,VmVp4f$cVb,use="complete.obs")#[1] 0.9997046
# # cor(VmVp4f$cVm,VmVp4f$cVb,use="complete.obs") #[1] 0.9994929
# # 
# # 
