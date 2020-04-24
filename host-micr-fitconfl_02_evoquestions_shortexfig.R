#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to produce illustrative graphics
	#run simulations across specific scenarios for illustrative purposes
##


trait.v1 <- sort(rnorm(1000,mean=1.5,sd=1))
teste1.micr   <- univar.fit( z=trait.v1, zopt=2, sd.fit=2)
teste1.plnt   <- univar.fit( z=trait.v1, zopt=3, sd.fit=0.75)
teste1.micrff   <- 0.6*teste1.plnt + (1-0.6)*teste1.micr
teste1.plntff   <- (1-0.6)*teste1.plnt + 0.6*teste1.micr
con.e1 <- c(trait.v1[which(teste1.micrff == max(teste1.micrff))],trait.v1[which(teste1.plntff == max(teste1.plntff))])
testf1.plnt   <- univar.fit( z=trait.v1, zopt=3, sd.fit=0.75)
testf1.micr   <- univar.fit( z=trait.v1, zopt=2, sd.fit=0.75)
testf1.plntff   <- 0.6*testf1.plnt + (1-0.6)*testf1.micr
testf1.micrff   <- (1-0.6)*testf1.plnt + 0.6*testf1.micr
con.f1 <- c(trait.v1[which(testf1.micrff == max(testf1.micrff))],trait.v1[which(testf1.plntff == max(testf1.plntff))])

trait.v2 <- sort(rnorm(1000,mean=2.5,sd=1))
teste2.micr   <- univar.fit( z=trait.v2, zopt=2, sd.fit=2)
teste2.plnt   <- univar.fit( z=trait.v2, zopt=3, sd.fit=0.75)
teste2.micrff   <- 0.6*teste2.plnt + (1-0.6)*teste2.micr
teste2.plntff   <- (1-0.6)*teste2.plnt + 0.6*teste2.micr
con.e2 <- c(trait.v2[which(teste2.micrff == max(teste2.micrff))],trait.v2[which(teste2.plntff == max(teste2.plntff))])
testf2.plnt   <- univar.fit( z=trait.v2, zopt=3, sd.fit=0.75)
testf2.micr   <- univar.fit( z=trait.v2, zopt=2, sd.fit=0.75)
testf2.plntff   <- 0.6*testf2.plnt + (1-0.6)*testf2.micr
testf2.micrff   <- (1-0.6)*testf2.plnt + 0.6*testf2.micr
con.f2 <- c(trait.v2[which(testf2.micrff == max(testf2.micrff))],trait.v2[which(testf2.plntff == max(testf2.plntff))])

which(testf2.plntff == max(testf2.plntff))


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/RelativeFitness_DirectandFeedbacks.pdf",height=5,width=5)
layout(matrix(1:6, ncol=2,byrow=F))
par(mar=c(3,3,0,0))
par(oma=c(0.5,2,3,1))
	plot(function(x) dnorm(x,mean=1,sd=1),min(trait.v1),max(trait.v1),xlim=c(-2,6),yaxt="n",ylab="",xlab="",lty=3)
		abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=2); abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=2)
		mtext("Frequency",side=2,line=1)
		mtext("Far from optima",side=3,line=1)
	plot(teste1.plnt~trait.v1,pch=NA,xlim=c(-2,6),ylim=c(0,0.0035),ylab="",xlab="")
# 		polygon(c(con.e1,rev(con.e1)), y = rep(c(-1,1),each=2),border=NA,col=rgb(0,0,0,alpha=0.5))
		abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=2); abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=2)
		lines(teste1.micr ~trait.v1,col=rgb(0.5,0,0.5)); lines(teste1.plnt ~trait.v1,col=rgb(0,0.5,0))
		lines(teste1.micrff ~trait.v1,col=rgb(0.5,0,0.5),lty=2); lines(teste1.plntff ~trait.v1,col=rgb(0,0.5,0),lty=2)
	#	mtext("Wm >> Wp",side=2,line=2)
		mtext("Trait affects B more",side=2,line=2,cex=0.75)
		mtext("Relative Fitness",side=2,line=3, adj=15)
	plot(testf1.plnt~trait.v1,pch=NA,xlim=c(-2,6),ylim=c(0,0.0035),ylab="",xlab="")
		polygon(c(con.f1,rev(con.f1)), y = rep(c(-1,1),each=2),border=NA,col=rgb(0,0,0,alpha=0.25))
		abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=2); abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=2)
	#	abline(v=3,col=rgb(0,0.5,0),lty=3); abline(v=2,col=rgb(0.5,0,0.5),lty=3)
		lines(testf1.micr ~trait.v1,col=rgb(0.5,0,0.5)); lines(testf1.plnt ~trait.v1,col=rgb(0,0.5,0))
		lines(testf1.micrff ~trait.v1,col=rgb(0.5,0,0.5),lty=2); lines(testf1.plntff ~trait.v1,col=rgb(0,0.5,0),lty=2)
# 		mtext("Wm = Wp",side=2,line=2)
		mtext("Equal link to fitness",side=2,line=2,cex=0.75)
	#####
	plot(function(x) dnorm(x,mean=2.5,sd=1),min(trait.v2),max(trait.v2),xlim=c(-2,6),yaxt="n",ylab="",xlab="",lty=3)
		abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=2); abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=2)
#		mtext("Frequency",side=2,line=1)
		mtext("Bewteen optima",side=3,line=1)
		text(x=c(0.75,4.25),y=c(0.3,0.3), labels =c(expression(Z[opt[A]]==2),expression(Z[opt[B]]==3)))
	plot(teste2.plnt~trait.v2,pch=NA,xlim=c(-2,6),ylim=c(0,0.0035),ylab="",xlab="")
# 		polygon(c(con.e2,rev(con.e2)), y = rep(c(-1,1),each=2),border=NA,col=rgb(0,0,0,alpha=0.5))
		abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=2); abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=2)
		lines(teste2.micr ~trait.v2,col=rgb(0.5,0,0.5)); lines(teste2.plnt ~trait.v2,col=rgb(0,0.5,0))
		lines(teste2.micrff ~trait.v2,col=rgb(0.5,0,0.5),lty=2); lines(teste2.plntff ~trait.v2,col=rgb(0,0.5,0),lty=2)
		legend(-2.5,0.0035,c("Species A","Species B","Conflict zone"),fill=c(rgb(0.5,0,0.5),rgb(0,0.5,0),rgb(0,0,0,alpha=0.25)),bty="n",border=NA)
	plot(testf2.plnt~trait.v2,pch=NA,xlim=c(-2,6),ylim=c(0,0.0035),ylab="",xlab="")
		polygon(c(con.f2,rev(con.f2)), y = rep(c(-1,1),each=2),border=NA,col=rgb(0,0,0,alpha=0.25))
		abline(v=3,col=rgb(0,0.5,0,alpha=0.5),lwd=2); abline(v=2,col=rgb(0.5,0,0.5,alpha=0.5),lwd=2)
		lines(testf2.micr ~trait.v2,col=rgb(0.5,0,0.5)); lines(testf2.plnt ~trait.v2,col=rgb(0,0.5,0))
		lines(testf2.micrff ~trait.v2,col=rgb(0.5,0,0.5),lty=2); lines(testf2.plntff ~trait.v2,col=rgb(0,0.5,0),lty=2)
		mtext("Trait value",side=1,line=2.5,adj=-0.5)
		legend(-2.5,0.0035,c("Direct only","+ feedback"),lty=c(1,2),bty="n")
dev.off()



# source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

gens <- 300
#matched optima, no ff
test.4a <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)#at about .1 for fiterr relative to sdfit of 1 is when error in fitness starts to obscure fitness differences, not calculated. just apprx guess.	
#mattched optima, ff
test.4aff <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)#at about .1 for fiterr relative to sdfit of 1 is when error in fitness starts to obscure fitness differences, not calculated. just apprx guess.	
#plant direct very weak
test.4b <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=10,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#plant direct very weak, indirect stronger
test.4bff <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=10,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#conflict, microbe link stronger,  direct links only
test.4d <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=2,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#conflict, microbe link stronger,  direct and indirect links
test.4dff <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=2,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#conflict, plant link stronger,  direct links only
test.4e <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=2,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#conflict, plant link stronger,  direct and indirect links
test.4eff <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=2,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#conflict, links =,  direct only
test.4f <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#conflict, links =,  direct and indirect links
test.4fff <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)


##
# #fourdemosims <- list(test.4a,test.4b,test.4d,test.4e)
# #save(fourdemosims,file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_PFF_expDFE.R",sep=""))
# fourdmosims <- load(paste(Sys.getenv("SCRATCH"),"/Simulation Results_PFF_expDFE.R",sep=""))
# test.4a <- fourdemosims[[1]]
# test.4b <- fourdemosims[[2]]
# test.4d <- fourdemosims[[3]]
# test.4e <- fourdemosims[[4]]

#pg <- colorRampPalette(c( rgb(0.5,0,0.5), rgb(1,1,1), rgb(0,0.5,0) ))


titles <- c("Zoptp = Zoptm, alphas = 1","Zoptp = Zoptm, alphas = 0.6",
			"~0 link to Fp, alphas = 1", "~0 link to Fp, alphas = 0.6",
			"Fm link > Fp link, alphas = 1", "Fm link > Fp link, alphas = 0.6",
			"Fp link > Fm link, alphas = 1", "Fp link > Fm link, alphas = 0.6",
			"Fp link = Fm link, alphas = 1", "Fp link = Fm link, alphas = 0.6")

# titles <- c("~0 link to microbe fitness","~0 link to plant fitness","microbe link > plant link","plant link > microbe link")

gens <- 200
# pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_tendemos.pdf",sep=""),height=6,width=6)
layout(matrix(1:10,ncol=5,byrow=F))
# layout(matrix(1:8,ncol=4,byrow=F))
par(mar=c(3,3,1,1))
par(oma=c(1,0,1,0))
windowplot(1,gens+1, 3,test.4b,ylim=c(-3,6),titles[3])
	abline(h=3,col=rgb(0.5,0,0.5),lty=2); abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
	abline(h=0)
windowplot(1,gens+1, 3,test.4bff,ylim=c(-3,6),titles[4])
	abline(h=3,col=rgb(0.5,0,0.5),lty=2); abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
	abline(h=0)
	legend(gens*-0.05,y=0,c("Plant","Microbe","Expressed"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
windowplot(1,gens+1,3,test.4d,ylim=c(-3,6),titles[5])
	abline(h=2,col=rgb(0,0.5,0),lty=2); abline(h=3,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot(1,gens+1,3,test.4dff,ylim=c(-3,6),titles[6])
	abline(h=2,col=rgb(0,0.5,0),lty=2); abline(h=3,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
	legend(gens*-0.05,y=0,c("Plant optima","Microbe optima","Midpoint"), lty = 2, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
windowplot(1,gens+1,3,test.4e,ylim=c(-3,6),titles[7])
	abline(h=3,col=rgb(0,0.5,0),lty=2); abline(h=2,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot(1,gens+1,3,test.4eff,ylim=c(-3,6),titles[8])
	abline(h=3,col=rgb(0,0.5,0),lty=2); abline(h=2,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot(1,gens+1,3,test.4f,ylim=c(-3,6),titles[9])
	abline(h=3,col=rgb(0,0.5,0),lty=2); abline(h=2,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot(1,gens+1,3,test.4fff,ylim=c(-3,6),titles[10])
	abline(h=3,col=rgb(0,0.5,0),lty=2); abline(h=2,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot(1,gens+1,3,test.4a,ylim=c(-3,6),titles[1])
	abline(h=3,col=rgb(0,0.5,0),lty=2);  abline(h=3.1,col=rgb(0.5,0,0.5),lty=2)#col=rgb(0,0,1),lty=2)
	abline(h=0)
windowplot(1,gens+1,3,test.4aff,ylim=c(-3,6),titles[2])
	abline(h=3,col=rgb(0,0.5,0),lty=2);  abline(h=3.1,col=rgb(0.5,0,0.5),lty=2)#col=rgb(0,0,1),lty=2)
	abline(h=0)
# dev.off()


# pdf(paste(Sys.getenv("SCRATCH"),"/Simulation Results_fourdemos.pdf",sep=""),height=6,width=6)
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/Simulation Results_fourdemos_n.pdf",height=5,width=5)
layout(matrix(1:4,ncol=2,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,0,0))
windowplot(1,gens+1, 3,test.4b,ylim=c(-3,6),"",xlabs="")#titles[3])
	abline(h=3,col=rgb(0.5,0,0.5),lty=2); abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
	abline(h=0)
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No link to plant fitness",side=3)
windowplot(1,gens+1, 3,test.4bff,ylim=c(-3,6),"")#"titles[4]")
	abline(h=3,col=rgb(0.5,0,0.5),lty=2); abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
	abline(h=0)
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	legend(gens*-0.05,y=0,c("Plant","Microbe","Expressed"), fill=c(rgb(0,0.5,0,alpha=0.5),rgb(0.5,0,0.5,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")#"Holobiont"	
windowplot(1,gens+1,3,test.4f,ylim=c(-3,6),"",xlabs="",ylabs="")#titles[9])
	abline(h=3,col=rgb(0,0.5,0),lty=2); abline(h=2,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
	mtext("Equal links to fitness",side=3)
windowplot(1,gens+1,3,test.4fff,ylim=c(-3,6),"",ylabs="")#titles[10])
	abline(h=3,col=rgb(0,0.5,0),lty=2); abline(h=2,col=rgb(0.5,0,0.5),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
	legend(gens*-0.05,y=0,c("Plant optima","Microbe optima","Midpoint"), lty = 2, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
dev.off()



#demo figure for total variance "G", variance G proportions of M and P, distance from optimum, local dynamics?
#take two simulations from above figure, and plot each below
FC4b <- getfitcon(10, gens+1, 1, test.4b,zoP=2,zoM=3, wP=5, wM=0.75,pfP=0.6,pfM=0.6)   
VmVp4b <- extractVmVp(test.4b, 1,gens+1,1)
win4b  <- extractwinning(test.4b,first=1,last=gens+1,1,zoP=2,zoM= 3)
dyn4b <- extractDyn(test.4b,first=1,last=gens+1,10)
FC4e <- getfitcon(10, gens+1, 1, test.4e,zoP=3,zoM=2, wP=0.75, wM=1,pfP=0.6,pfM=0.6)   
VmVp4e <- extractVmVp(test.4e, 1,gens+1,1)
win4e  <- extractwinning(test.4e,first=1,last=gens+1,1,zoP=3,zoM= 2)
dyn4e <- extractDyn(test.4e,first=1,last=gens+1,10)
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
		lines( 1:length(VmVp4b$PVp),VmVp4b$PVp,col=rgb(0,0.5,0)) 
		lines(1:length(VmVp4b$PVm),VmVp4b$PVm , col=rgb(0.5,0,0.5)) 
	plot(VmVp4e$Vp ~ c(1:length(VmVp4e$Vp)),pch=NA,ylim=c(0,1),ylab="",xlab="")
		lines( 1:length(VmVp4e$PVp),VmVp4e$PVp,col=rgb(0,0.5,0))
		lines(1:length(VmVp4e$PVm),VmVp4e$PVm , col=rgb(0.5,0,0.5))
	plot(win4b$dP ~ c(1:length(win4b$dP)),pch=NA,ylim=c(0,max(c(win4b$dP,win4b$dM))),ylab="Distance from optima",xlab="")
		lines(win4b$dP~c(1:(gens+1)),col=rgb(0,0.5,0))
		lines(win4b$dM~c(1:(gens+1)),col=rgb(0.5,0,0.5))
	plot(win4e$dP ~ c(1:length(win4e$dP)),pch=NA,ylim=c(0,max(c(win4e$dP,win4e$dM))),ylab="",xlab="")
		lines(win4e$dP~c(1:(gens+1)),col=rgb(0,0.5,0))
		lines(win4e$dM~c(1:(gens+1)),col=rgb(0.5,0,0.5))
	plot(dyn4b$tcoefP ~ c(1:length(dyn4b$tcoefP)),pch=NA,ylim=c(min(c(dyn4b$tcoefP,dyn4b$tcoefM)),max(c(dyn4b$tcoefP,dyn4b$tcoefM))),ylab="Local trait change",xlab="generations")
		abline(h=0)
		lines(dyn4b$tcoefP~c(1:length(dyn4b$tcoefP)),col=rgb(0,0.5,0))
		lines(dyn4b$tcoefM~c(1:length(dyn4b$tcoefM)),col=rgb(0.5,0,0.5))	
	plot(dyn4e$tcoefP ~ c(1:length(dyn4e$tcoefP)),pch=NA,ylim=c(min(c(dyn4e$tcoefP,dyn4e$tcoefM)),max(c(dyn4e$tcoefP,dyn4e$tcoefM))),ylab="Local trait change",xlab="generations")
 		abline(h=0)
		lines(dyn4e$tcoefP~c(1:length(dyn4e$tcoefP)),col=rgb(0,0.5,0))
		lines(dyn4e$tcoefM~c(1:length(dyn4e$tcoefM)),col=rgb(0.5,0,0.5))
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
dynVE <- extractDyn(testvarE,first=1,last=gens+1,10)

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
		lines( 1:length(VmVpVE$PVp),VmVpVE$PVp,col=rgb(0,0.5,0))
		lines(1:length(VmVpVE$PVm),VmVpVE$PVm , col=rgb(0.5,0,0.5))
        plot(winVE$dP ~ c(1:length(winVE$dP)),pch=NA,ylim=c(0,max(c(winVE$dP,winVE$dM))),ylab="Distance from optima",xlab="")
		lines(winVE$dP~c(1:(gens+1)),col=rgb(0,0.5,0))
		lines(winVE$dM~c(1:(gens+1)),col=rgb(0.5,0,0.5))
        plot(dynVE$tcoefP ~ c(1:length(dynVE$tcoefP)),pch=NA,ylim=c(min(c(dynVE$tcoefP,dynVE$tcoefM)),max(c(dynVE$tcoefP,dynVE$tcoefM))),ylab="Local trait change",xlab="generations")
		abline(h=0)
		lines(dynVE$tcoefP~c(1:length(dynVE$tcoefP)),col=rgb(0,0.5,0))
		lines(dynVE$tcoefM~c(1:length(dynVE$tcoefM)),col=rgb(0.5,0,0.5))
	plot(zvecs$PlantZ~c(1:length(zvecs$PlantZ)),ylim=c(-2,6),ylab="trait optima",xlab="generations",pch=NA)
		lines(zvecs$PlantZ~c(1:length(zvecs$PlantZ)),col=rgb(0,0.5,0))
        lines(zvecs$MicrZ~c(1:length(zvecs$MicrZ)),col=rgb(0.5,0,0.5))
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
                abline(h=3,col=rgb(0,0.5,0),lty=2)
                abline(h=2,col=rgb(0.5,0,0.5),lty=2)
                abline(h=0)
        windowplot(1,gens+1,5,testmanyMb,ylim=c(-2,6),"high cost of free-living, stronger link to plant fitness")
                abline(h=3,col=rgb(0,0.5,0),lty=2)
                abline(h=2,col=rgb(0.5,0,0.5),lty=2)
                abline(h=0)
	 plot(FCmMa$fitnesscorrelation~c(10:(gens+1)),main="",ylab="Fitness correlation",xlab="",ylim=c(-1,1),pch=NA)
                lines(FCmMa$fitnesscorrelation~c(10:(gens+1)),lty=2)
                abline(h=0)
        plot(FCmMb$fitnesscorrelation~c(10:(gens+1)),main="",ylab="",xlab="",ylim=c(-1,1),pch=NA)
                lines(FCmMb$fitnesscorrelation~c(10:(gens+1)),lty=2)
                abline(h=0)
       plot(VmVpmMa$Vp ~ c(1:length(VmVpmMa$Vp)),pch=NA,ylim=c(0,1),ylab="Proportion of genetic variance",xlab="generations")
                lines( 1:length(VmVpmMa$PVp),VmVpmMa$PVp,col=rgb(0,0.5,0))
                lines(1:length(VmVpmMa$PVm),VmVpmMa$PVm , col=rgb(0.5,0,0.5))
       plot(VmVpmMb$Vp ~ c(1:length(VmVpmMb$Vp)),pch=NA,ylim=c(0,1),ylab="",xlab="generations")
                lines( 1:length(VmVpmMb$PVp),VmVpmMb$PVp,col=rgb(0,0.5,0))
                lines(1:length(VmVpmMb$PVm),VmVpmMb$PVm , col=rgb(0.5,0,0.5))
dev.off()




 #for each scenario, show fitness-trait relationship for plants and microbes (not exactly relative fitness)
# trait.v <- seq(from=0,to=5,length.out=1000)
# trait.v <- c(seq(from=0,to=1,length.out=800),seq(from=1,to=5,length.out=200))
# trait.v <- sort(rnorm(1000,mean=1,sd=1))
# # zopt and w identical
# testa.all   <- univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# #zopt differ, plant almost 0, w/ and w/o ff
# testb.plnt   <- univar.fit( z=trait.v, zopt=2, sd.fit=10)
# testb.micr   <- univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# testb.plntff   <- 0.6*univar.fit( z=trait.v, zopt=2, sd.fit=10) + (1-0.6)*univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# testb.micrff   <- (1-0.6)*univar.fit( z=trait.v, zopt=2, sd.fit=10) + 0.6*univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# #zopt differ, micr link stronger, w/ and w/o ff
# testd.plnt   <- univar.fit( z=trait.v, zopt=2, sd.fit=2)
# testd.micr   <- univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# testd.plntff   <- 0.6*univar.fit( z=trait.v, zopt=2, sd.fit=2) + (1-0.6)*univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# testd.micrff   <- (1-0.6)*univar.fit( z=trait.v, zopt=2, sd.fit=2) + 0.6*univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# #zopt differ, micr link stronger, w/ and w/o ff
# teste.micr   <- univar.fit( z=trait.v, zopt=2, sd.fit=2)
# teste.plnt   <- univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# teste.micrff   <- 0.6*univar.fit( z=trait.v, zopt=2, sd.fit=2) + (1-0.6)*univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# teste.plntff   <- (1-0.6)*univar.fit( z=trait.v, zopt=2, sd.fit=2) + 0.6*univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# #zopt differ, w idential
# testf.plnt   <- univar.fit( z=trait.v, zopt=3, sd.fit=0.75)
# testf.micr   <- univar.fit( z=trait.v, zopt=2, sd.fit=0.75)
# testf.plntff   <- 0.6*univar.fit( z=trait.v, zopt=3, sd.fit=0.75) + (1-0.6)*univar.fit( z=trait.v, zopt=2, sd.fit=0.75)
# testf.micrff   <- (1-0.6)*univar.fit( z=trait.v, zopt=3, sd.fit=0.75) + 0.6*univar.fit( z=trait.v, zopt=2, sd.fit=0.75)
# 
# #plot(function(x) dnorm(x,mean=1,sd=1),min(trait.v),max(trait.v))
# # plot(trait.v,testa.all,pch=NA)
# # 	lines(testa.all~trait.v)
# # 	abline(v=3,lty=3)
# plot(testb.micr~trait.v,pch=NA)
# 	abline(v=2,col=rgb(0,0.5,0),lty=3); abline(v=3,col=rgb(0.5,0,0.5),lty=3)
# 	lines(testb.micr ~trait.v,col=rgb(0.5,0,0.5)); lines(testb.plnt ~trait.v,col=rgb(0,0.5,0))
# 	lines(testb.micrff ~trait.v,col=rgb(0.5,0,0.5),lty=2); lines(testb.plntff ~trait.v,col=rgb(0,0.5,0),lty=2)
# plot(testd.micr~trait.v,pch=NA)
# 	abline(v=2,col=rgb(0,0.5,0),lty=3); abline(v=3,col=rgb(0.5,0,0.5),lty=3)
# 	lines(testd.micr ~trait.v,col=rgb(0.5,0,0.5)); lines(testd.plnt ~trait.v,col=rgb(0,0.5,0))
# 	lines(testd.micrff ~trait.v,col=rgb(0.5,0,0.5),lty=2); lines(testd.plntff ~trait.v,col=rgb(0,0.5,0),lty=2)
# plot(teste.plnt~trait.v,pch=NA)
# 	abline(v=3,col=rgb(0,0.5,0),lty=3); abline(v=2,col=rgb(0.5,0,0.5),lty=3)
# 	lines(teste.micr ~trait.v,col=rgb(0.5,0,0.5)); lines(teste.plnt ~trait.v,col=rgb(0,0.5,0))
# 	lines(teste.micrff ~trait.v,col=rgb(0.5,0,0.5),lty=2); lines(teste.plntff ~trait.v,col=rgb(0,0.5,0),lty=2)
# plot(testf.plnt~trait.v,pch=NA)
# 	abline(v=3,col=rgb(0,0.5,0),lty=3); abline(v=2,col=rgb(0.5,0,0.5),lty=3)
# 	lines(testf.micr ~trait.v,col=rgb(0.5,0,0.5)); lines(testf.plnt ~trait.v,col=rgb(0,0.5,0))
# 	lines(testf.micrff ~trait.v,col=rgb(0.5,0,0.5),lty=2); lines(testf.plntff ~trait.v,col=rgb(0,0.5,0),lty=2)
# 

