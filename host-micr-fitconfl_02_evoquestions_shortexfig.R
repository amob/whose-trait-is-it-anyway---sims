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

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/Fitcor_DirectandFeedbacks.pdf",height=3.5,width=3.5)
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(3,3,1,1))
plot(teste1.plnt~teste1.micr,ylim=c(0,0.004),xlim=c(0,0.0025),pch=1,cex=0.5,col=rgb(0,0,0,alpha=1),ylab="",xlab="",xaxt="n")
	points(teste1.plntff~teste1.micrff,pch=8,cex=0.5,col=rgb(0.5,0.5,0.5,alpha=1))
	mtext("Plant relative fitness",side=2,line=2.5,adj=1.75)
plot(testf1.plnt~testf1.micr,ylim=c(0,0.004),xlim=c(0,0.0025),pch=1,cex=0.5,col=rgb(0,0,0,alpha=1),ylab="",xlab="",xaxt="n") #pch=16,cex=0.5,col=rgb(0,0,0,alpha=0.25),ylab="",xlab="",xaxt="n",yaxt="n")
	points(testf1.plntff~testf1.micrff,pch=16,cex=0.5,col=rgb(0.5,0.5,0.5,alpha=1))#,col=rgb(1,0,0,alpha=0.25))
plot(teste2.plnt~teste2.micr,ylim=c(0,0.004),xlim=c(0,0.0025),pch=1,cex=0.5,col=rgb(0,0,0,alpha=1),ylab="",xlab="",xaxt="n")
	points(teste2.plntff~teste2.micrff,pch=16,cex=0.5,col=rgb(0.5,0.5,0.5,alpha=1))
plot(testf2.plnt~testf2.micr,ylim=c(0,0.004),xlim=c(0,0.0025),pch=1,cex=0.5,col=rgb(0,0,0,alpha=1),ylab="",xlab="",xaxt="n")
	points(testf2.plntff~testf2.micrff,pch=16,cex=0.5,col=rgb(0.5,0.5,0.5,alpha=1))
	mtext("Microbe relative fitness",side=1,line=2.5,adj=1.75)
# 	abline(lm(testf2.plnt~testf2.micr),col=rgb(0,0,0),lty=2)
# 	abline(lm(testf2.plntff~testf2.micrff),col=rgb(1,0,0),lty=2)
legend(-0.00025,y=0.004,c("Direct only","+ feedbacks"),fill=c(rgb(0,0,0),rgb(0.5,0.5,0.5)),bty="n",border=NA )
dev.off()


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

gens <- 200
#matched optima, no ff
test.4a <-		sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)#
#mattched optima, ff
test.4aff <- 	sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=3,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#plant direct very weak
test.4b <- 		sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=10,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#plant direct very weak, indirect stronger
test.4bff <- 	sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=10,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#conflict, microbe link stronger,  direct links only
test.4d <- 		sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=2,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#conflict, microbe link stronger,  direct and indirect links
test.4dff <-	 sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=2,zoM=3,wP=2,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#conflict, plant link stronger,  direct links only
test.4e <- 		sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=2,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#conflict, plant link stronger,  direct and indirect links
test.4eff <- 	sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=2,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
#conflict, links =,  direct only
test.4f <- 		sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 1, pfM=1,FLFC=0.1)
#conflict, links =,  direct and indirect links
test.4fff <- 	sim.cotrait(NP=100,NM=100,nlP=100,nlM=200,nlnP=3,nlnM=3,zoP=3,zoM=2,wP=0.75,wM=0.75,timesteps=gens,Lambda=25,mutprb=0.0005,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)


##
tendemosims <- list(test.4a,test.4aff,test.4b,test.4bff,test.4d,test.4dff,test.4e,test.4eff,test.4f,test.4fff)
save(tendemosims,file=paste(Sys.getenv("SCRATCH"),"/Simulation Results_PFF_expDFE.R",sep=""))


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
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/Simulation_Results_fourdemos_n.pdf",height=5,width=5)
layout(matrix(1:4,ncol=2,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,0,0))
windowplot(1,gens+1, 3,test.4b,ylim=c(-3,6),"",xlabs="")#titles[3])
	abline(h=3,col=rgb(0.5,0,0.5),lty=2); # abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
	abline(h=0)
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No link to plant fitness",side=3)
windowplot(1,gens+1, 3,test.4bff,ylim=c(-3,6),"")#"titles[4]")
	abline(h=3,col=rgb(0.5,0,0.5),lty=2); #abline(h=2,col=rgb(0,0.5,0),lty=2)#col=rgb(1,0,0),lty=2)
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



getrelfitandtrait <- function(simdat,gen,zoP,zoM,wP,wM,pfP,pfM){
	traitend <- rowSums(colSums(simdat$Plant[,,,gen])) + colSums(simdat$Microbe[,,gen])
	micrfit   <- univar.fit( z=traitend, zopt=zoM, sd.fit=wM)
	hostfit   <- univar.fit( z=traitend, zopt=zoP, sd.fit=wP)
	rHfit <-     pfP*hostfit + (1-pfP)*micrfit
	rMfit <- (1-pfM)*hostfit +     pfM*micrfit
	return(data.frame(traitend = traitend,rfitplnt=rHfit,rfitmicr=rMfit))
}

trait.a 	<- getrelfitandtrait(test.4a  ,gens+1,zoP=3,zoM=3,wP=0.75,wM=0.75,pfP=1,pfM=1)
trait.aff 	<- getrelfitandtrait(test.4aff,gens+1,zoP=3,zoM=3,wP=0.75,wM=0.75,pfP=0.6,pfM=0.6)
trait.b 	<- getrelfitandtrait(test.4b  ,gens+1,zoP=2,zoM=3,wP=10  ,wM=0.75,pfP=1,pfM=1)
trait.bff 	<- getrelfitandtrait(test.4bff,gens+1,zoP=2,zoM=3,wP=10  ,wM=0.75,pfP=0.6,pfM=0.6)
trait.d 	<- getrelfitandtrait(test.4d  ,gens+1,zoP=2,zoM=3,wP=2   ,wM=0.75,pfP=1,pfM=1)
trait.dff 	<- getrelfitandtrait(test.4dff,gens+1,zoP=2,zoM=3,wP=2   ,wM=0.75,pfP=0.6,pfM=0.6)
trait.e 	<- getrelfitandtrait(test.4e  ,gens+1,zoP=3,zoM=2,wP=0.75,wM=2   ,pfP=1,pfM=1)
trait.eff 	<- getrelfitandtrait(test.4eff,gens+1,zoP=3,zoM=2,wP=0.75,wM=2   ,pfP=0.6,pfM=0.6)
trait.f 	<- getrelfitandtrait(test.4f  ,gens+1,zoP=3,zoM=2,wP=0.75,wM=0.75,pfP=1,pfM=1)
trait.fff 	<- getrelfitandtrait(test.4fff,gens+1,zoP=3,zoM=2,wP=0.75,wM=0.75,pfP=0.6,pfM=0.6)

conend.e <- c(trait.eff$traitend[which(trait.eff$rfitmicr == max(trait.eff$rfitmicr))],trait.eff$traitend[which(trait.eff$rfitplnt== max(trait.eff$rfitplnt))])
conend.f <- c(trait.fff$traitend[which(trait.fff$rfitmicr == max(trait.fff$rfitmicr))],trait.fff$traitend[which(trait.fff$rfitplnt== max(trait.fff$rfitplnt))])

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/RelativeandDirectandCorrelation_sims_fourdemos_n.pdf",height=5,width=5)
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
par(oma=c(3,5,1,1))
plot( trait.b$rfitplnt~trait.b$traitend,col=rgb(0,0.5,0), pch=1,xlim=c(1,4),ylim=c(0,0.015),xaxt="n",cex=0.5)
	points( trait.b$rfitmicr~trait.b$traitend,col=rgb(0.5,0,0.5),pch=1 ,cex=0.5)
	points( trait.bff$rfitplnt~trait.bff$traitend,col=rgb(0,0.5,0),pch=8 ,cex=0.5)
	points( trait.bff$rfitmicr~trait.bff$traitend,col=rgb(0.5,0,0.5),pch=8 ,cex=0.5)
	mtext("Plant relative fitness",side=2,adj=-3,line=2)
	mtext("No link to plant fitness",side=2,line=3.5)
plot( trait.b$rfitplnt~trait.b$rfitmicr,col=rgb(0,0,0),pch=1 ,xlim=c(0,0.015) ,ylim=c(0,0.015),xaxt="n",yaxt="n",cex=0.5)
#	points( trait.bff$rfitplnt~trait.bff$rfitmicr,col=rgb(0.5,0.5,0.5) ,cex=0.5)
	points( trait.bff$rfitplnt~trait.bff$rfitmicr,pch=8 ,cex=0.5)
plot( trait.f$rfitplnt~trait.f$traitend,col=rgb(0,0.5,0), pch=1,xlim=c(1,4),ylim=c(0,0.015),cex=0.5)
	points( trait.f$rfitmicr~trait.f$traitend,col=rgb(0.5,0,0.5),pch=1 ,cex=0.5)
	points( trait.fff$rfitplnt~trait.fff$traitend,col=rgb(0,0.5,0),pch=8 ,cex=0.5)
	points( trait.fff$rfitmicr~trait.fff$traitend,col=rgb(0.5,0,0.5),pch=8 ,cex=0.5)
	mtext("Equal links to fitness",side=2,line=3.5)
	mtext("Trait Value",side=1, line=2,cex=1)
	legend(-0.05,y=0.04,c("Plant","Microbe"),bty="n",border=NA, fill=c(rgb(0,0.5,0),rgb(0.5,0,0.5)) )
plot( trait.f$rfitplnt~trait.f$rfitmicr,col=rgb(0,0,0) ,pch=1,xlim=c(0,0.015) ,ylim=c(0,0.015),yaxt="n",cex=0.5)
	points( trait.fff$rfitplnt~trait.fff$rfitmicr,pch=8 ,cex=0.5)#,col=rgb(0.5,0.5,0.5)
	legend(0,y=0.004,c("Direct only","+ feedbacks"),pch=c(1,8),bty="n",border=NA) #fill=c(rgb(0,0,0),rgb(0.5,0.5,0.5)) )
	mtext("Microbe relative fitness",side=1, line=2,cex=1)
dev.off()



#demo figure for total variance "G", variance G proportions of M and P, distance from optimum, local dynamics?
#take two simulations from above figure, and plot each below

VmVp4b <- extractVmVp(test.4b, 1,gens+1,1)
VmVp4e <- extractVmVp(test.4e, 1,gens+1,1)
VmVp4bff <- extractVmVp(test.4bff, 1,gens+1,1)
VmVp4f <- extractVmVp(test.4f, 1,gens+1,1)
VmVp4fff <- extractVmVp(test.4fff, 1,gens+1,1)

#weaker variance correlations with fitness feedbacks. interesting. but the plots make it look like this could depend on timing of events, which change between sims/scenarios
cor(VmVp4fff$Vp,VmVp4fff$Vm) #[1] 0.2660599
cor(VmVp4fff$Vp,VmVp4fff$Vb,use="complete.obs") #[1] 0.8572949
cor(VmVp4fff$Vm,VmVp4fff$Vb,use="complete.obs") #[1] 0.6107672
cor(VmVp4fff$cVp,VmVp4fff$cVm,use="complete.obs") #[1] 0.9656252
cor(VmVp4fff$cVm,VmVp4fff$cVb,use="complete.obs") #[1] 0.9980781
cor(VmVp4fff$cVp,VmVp4fff$cVb,use="complete.obs")#[1] 0.9791138
cor(VmVp4f$Vp,VmVp4f$Vm) #[1] 0.3198037
cor(VmVp4f$Vp,VmVp4f$Vb) #[1] 0.9296124
cor(VmVp4f$Vm,VmVp4f$Vb) #[1] 0.5733526
cor(VmVp4f$cVp,VmVp4f$cVm,use="complete.obs")#[1] 0.998429
cor(VmVp4f$cVp,VmVp4f$cVb,use="complete.obs")#[1] 0.9997046
cor(VmVp4f$cVm,VmVp4f$cVb,use="complete.obs") #[1] 0.9994929

win4bff  <- extractwinning(test.4bff,first=1,last=gens+1,1,zoP=2,zoM= 3)
win4f  <- extractwinning(test.4f,first=1,last=gens+1,1,zoP=3,zoM= 2)
win4fff  <- extractwinning(test.4fff,first=1,last=gens+1,1,zoP=3,zoM= 2)
win4e  <- extractwinning(test.4e,first=1,last=gens+1,1,zoP=3,zoM= 2)
win4b  <- extractwinning(test.4b,first=1,last=gens+1,1,zoP=2,zoM= 3)

which(win4f$dM<0.01)

rangepheno.4b <- sapply(2:(gens+2), function(t) range( rowSums(colSums(test.4b$Plant[,,,t])) + colSums(test.4b$Microbe[,,t]) )  )
rangepheno.4bff <- sapply(2:(gens+2), function(t) range( rowSums(colSums(test.4bff$Plant[,,,t])) + colSums(test.4bff$Microbe[,,t]) )  )
meanpheno.4b <- sapply(2:(gens+2), function(t) mean( rowSums(colSums(test.4b$Plant[,,,t])) + colSums(test.4b$Microbe[,,t]) )  )
meanpheno.4bff <- sapply(2:(gens+2), function(t) mean( rowSums(colSums(test.4bff$Plant[,,,t])) + colSums(test.4bff$Microbe[,,t]) )  )
meanpheno.4f <- sapply(2:(gens+2), function(t) mean( rowSums(colSums(test.4f$Plant[,,,t])) + colSums(test.4f$Microbe[,,t]) )  )
meanpheno.4fff <- sapply(2:(gens+2), function(t) mean( rowSums(colSums(test.4fff$Plant[,,,t])) + colSums(test.4fff$Microbe[,,t]) )  )


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/Variance_Results_fourdemos_n.pdf",height=4.5,width=5.5)
layout(matrix(1:4,ncol=2,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,1,0))
	plot(VmVp4b$Vp ~ c(1:length(VmVp4b$Vp)),pch=NA,ylim=c(0,0.12),ylab="",xlab="")
		lines( 1:length(VmVp4b$Vp),VmVp4b$Vp,col=rgb(0,0.5,0)) 
		lines(1:length(VmVp4b$Vm),VmVp4b$Vm , col=rgb(0.5,0,0.5)) 
		lines(1:length(VmVp4b$Vb),VmVp4b$Vb , col=rgb(0.5,0.5,0.5)) 
#		abline(v=min(which(win4b$dM < 0.1)))# first generation where the expressed value is pretty close to the microbe optima
# 		abline(v=min(which(rangepheno.4b[2,]>=3)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=3)# first generation does any individual in the pop express the microbe optima or beyond i
	 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
		mtext("No link to plant fitness",side=3,line=0.5)
		mtext("genetic variance",side=2,line=2,adj=-1.75)
	plot(VmVp4bff$Vp ~ c(1:length(VmVp4bff$Vp)),pch=NA,ylim=c(0,0.12),ylab="",xlab="")
		lines( 1:length(VmVp4bff$Vp),VmVp4bff$Vp,col=rgb(0,0.5,0)) 
		lines(1:length(VmVp4bff$Vm),VmVp4bff$Vm , col=rgb(0.5,0,0.5)) 
		lines(1:length(VmVp4bff$Vb),VmVp4bff$Vb , col=rgb(0.5,0.5,0.5)) 
		abline(v=min(which(meanpheno.4bff>=3)), col=rgb(0.5,0,0.5,alpha=0.5),lwd=5)# first generation does microbe mean reach its optima or beyond
		mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	plot(VmVp4f$Vp ~ c(1:length(VmVp4f$Vp)),pch=NA,ylim=c(0,0.12),ylab="",xlab="")
		lines( 1:length(VmVp4f$Vp),VmVp4f$Vp,col=rgb(0,0.5,0))
		lines(1:length(VmVp4f$Vm),VmVp4f$Vm , col=rgb(0.5,0,0.5))
		lines(1:length(VmVp4f$Vb),VmVp4f$Vb , col=rgb(0.5,0.5,0.5))
		abline(v=which(win4f$dM==min(win4f$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
		text(x=85,y=0.11, labels =expression(Z[opt[M]]~reached))
		mtext("Equal links to fitness",side=3,line=0.5)
	plot(VmVp4fff$Vp ~ c(1:length(VmVp4fff$Vp)),pch=NA,ylim=c(0,0.12),ylab="",xlab="")
		lines( 1:length(VmVp4f$Vp),VmVp4fff$Vp,col=rgb(0,0.5,0))
		lines(1:length(VmVp4f$Vm),VmVp4fff$Vm , col=rgb(0.5,0,0.5))
		lines(1:length(VmVp4fff$Vb),VmVp4fff$Vb , col=rgb(0.5,0.5,0.5))
		abline(v=which(win4fff$dM==min(win4fff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
		mtext("generations",side=1,line=2,adj=-0.75)
		legend(gens*0.2,y=0.11,c("Plant portion","Microbe portion","Interacting"),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
dev.off()

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/coefV_Results_fourdemos_n.pdf",height=4.5,width=5.5)
layout(matrix(1:4,ncol=2,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,1,0))
	plot(VmVp4b$cVp ~ c(1:length(VmVp4b$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
		lines( 1:length(VmVp4b$cVp),VmVp4b$cVp,col=rgb(0,0.5,0)) 
		lines(1:length(VmVp4b$cVm),VmVp4b$cVm , col=rgb(0.5,0,0.5)) 
		lines(1:length(VmVp4b$cVb),VmVp4b$cVb , col=rgb(0.5,0.5,0.5)) 
#		abline(v=min(which(win4b$dM < 0.1)))# first generation where the expressed value is pretty close to the microbe optima
# 		abline(v=min(which(rangepheno.4b[2,]>=3)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=3)# first generation does any individual in the pop express the microbe optima or beyond i
	 	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
		mtext("No link to plant fitness",side=3,line=0.5)
		mtext("genetic variance scaled by expressed mean",side=2,line=2,adj=1.08)
	plot(VmVp4bff$cVp ~ c(1:length(VmVp4bff$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
		lines( 1:length(VmVp4bff$cVp),VmVp4bff$cVp,col=rgb(0,0.5,0)) 
		lines(1:length(VmVp4bff$cVm),VmVp4bff$cVm , col=rgb(0.5,0,0.5)) 
		lines(1:length(VmVp4bff$cVb),VmVp4bff$cVb , col=rgb(0.5,0.5,0.5)) 
		abline(v=min(which(meanpheno.4bff>=3)), col=rgb(0.5,0,0.5,alpha=0.5),lwd=5)# first generation does microbe mean reach its optima or beyond
		mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	plot(VmVp4f$cVp ~ c(1:length(VmVp4f$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
		lines( 1:length(VmVp4f$cVp),VmVp4f$cVp,col=rgb(0,0.5,0))
		lines(1:length(VmVp4f$cVm),VmVp4f$cVm , col=rgb(0.5,0,0.5))
		lines(1:length(VmVp4f$cVb),VmVp4f$cVb , col=rgb(0.5,0.5,0.5))
		abline(v=which(win4f$dM==min(win4f$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
		text(x=85,y=0.12, labels =expression(Z[opt[M]]~reached))
		mtext("Equal links to fitness",side=3,line=0.5)
	plot(VmVp4fff$cVp ~ c(1:length(VmVp4fff$cVp)),pch=NA,ylim=c(0,0.15),ylab="",xlab="")
		lines( 1:length(VmVp4f$cVp),VmVp4fff$cVp,col=rgb(0,0.5,0))
		lines(1:length(VmVp4f$cVm),VmVp4fff$cVm , col=rgb(0.5,0,0.5))
		lines(1:length(VmVp4fff$cVb),VmVp4fff$cVb , col=rgb(0.5,0.5,0.5))
		abline(v=which(win4fff$dM==min(win4fff$dM)),col=rgb(0.5,0,0.5,alpha=0.5),lwd=5) #at which generation does the lower optima partner cross its optimum? only makes sense for f and fff, and therefor M
		mtext("generations",side=1,line=2,adj=-0.75)
		legend(gens*0.2,y=0.15,c("Plant portion","Microbe portion","Interacting"),lty=1, col = c(rgb(0,0.5,0),rgb(0.5,0,0.5),rgb(0,0,0,alpha=0.5)), bty="n")
dev.off()








###############stopped here with updates###
#NO SHIFTING OPTIMA HAS BEEN RUN.

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

