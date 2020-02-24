#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

######what of the questions we discussed can the current simulator accomplish?
	###how much conflict does there have to be for us to measure/map it?
		#Not tested.
	##variation in link between how much of trait is controlled by one plant/micr or the other
		#YES (sort of, with number of loci,), but I think fewer loci can eventually accumulate strong effects per locus over time
	##if link of trait to fitness for one partner is stronger , does this affect which genome will accumulate heritability?
		#YES, and answer is yes
	##what if there are subtle conflicts at some loci, but overall fitness alignment at most? i.e. mutations at loci in either genome same sign eff on fitness of both
		#Not tested, I think this requires multiple traits? or a portion of fitness that is dependent on your partner having good fitness? -- i.e. the fitness of your partner impacts your own fitness
		#I this might be similar to just turning down the strength of the link between the trait and fitness. but not sure.

#flexibility allowed
	#pop size of micr and plant can differ? 
			##FIXED, but new questions about what to do with the extra microbes. 
	#diploid plant stage under selection, rather than haploid  FIXED.
	#there's currently no correlation between plant and microbe fitness, which should probably add, but that will compete with correlation between plant pheno and microbe fitness.
		#FIXED. have added pos fitness feedbacks., some fitness error remains
		#now some proportion of focal fitness comes from partner fitness, and some percentage from having the right phenotype
		#there is however, still a small portion of random fitness variance (i.e. unrelated to trait values )
	#exponenetial distribution of effect sizes of mutations? #ALLOWED.

##flexibility for further consideration:
	## Dominance?
	#Mapping issue: keep track of alleles? -- maybe I can just use unique effect size values within loci? although there might be similar ways to arrive at the same effect sizes....
		# could say that each mutational change in an allele's history is a new site (inf sites) but that all within a locus are non-recombining?
	#more than one microbe species?
	#way to save fewer timesteps while still keeping track of above? would save on memory...JUST RUN ON NIAGARA?
	#what to do with the extra microbes. microbes without plants still have some fitness? multiple partners per plant?
	#is exponential really the right DFE?
	#a lande-arnold multivariate traits perspective: different beta vectors for plant microbe fitness....what about a mutational inputs mv matrix too?
			#so this is two part: 1, multiple traits and a mv phenotypic optima. 
			# and 2, trait covariance -- main way this will occur is if mutational inputs have a correlated (?what about just also but independent) effect on other traits
	#microbes faster generation time?

####################
###priority goals
	##neutral loci
	##simulate several populations
	##draw individuals from those populations to genotype and run experiments on. -- calculate trait values
################	

	
library(abind)

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

###The current code.

#fitness function for trait, zopt
#based on normal distribution, but instead of pdf summing to one, fitness is 1 at the optimum, just remove the 1/sqrt(2*pi*sd) in the normal P.D.F., not sure who to cite for deterministic part of this formula, I feel like I have seen it a lot...
univar.fit <- function(z, zopt, sd.fit) {    rawfit <- exp( -1* ((z-zopt)^2) / (2 * (sd.fit^2)) )
											if(sum(rawfit>0)){
												adjfit <- rawfit /sum(rawfit)
											} else {
												adjfit <- rawfit
											}
										return(adjfit) }  # + rnorm(length(z), mean=0,sd=fiterr)  ) }
#see lecorre and kramer 2012, for example, but this is commonly not cited., looks like is from Haldane 1954

# multivar.fit <- function(z, zopt, sd.fit,fiterr) { 
# 
#  dmvnorm(z,means=zopt  ,sigma = diag(sd.fit))
# #   sapply(1:length(z), function(k) exp( -1* ((z[k]-zopt[k])^2) / (2 * (sd.fit[k]^2)) ) 
#    
#    + rnorm(length(z), mean=0,sd=fiterr) ) 
# 
# }
# #the same function, sort of, but vectorized, zopt has to be a matrix?
# #I found 
# # zdiff <- (z - zopt)^2
# #exp( -(1/2) * t(zdiff)*inv(V)*zdiff  )
# # Evolution. 2003 Aug;57(8):1761-75.
# #Multivariate stabilizing selection and pleiotropy in the maintenance of quantitative genetic variation.
# #Zhang XS1, Hill WG.
# 
# #from Zhang 2012: this one has a bit more info,,
# 
# #zeta^2 = 1/(2VS) and VS is the variance of the fitness profile on each trait
# #and from turelli 1984, which is cited. Vs = w2 + Ve; with the note : For convenience, it will be assumed hereafter that all measurements are scaled so that Ve = 1
# ##so I gues zeta is a vector?? but then from Zhang: Individuals are subject to real stabilizing selection, with independent and identical strength of selection on each trait, characterized by S = ζ^2 I
# # exp( -(1/2) * sum(1:ntraits,  (z[i]-zopt[i])^2 * (zeta^2)    ) )
# 
#
#so, sdfit controls the strength of the peak in the normal distribution, when higher, peaks are flatter, and genetic differences in fitness can be overwhelmed by random factors, including non-genetic sources of fitness differences: fiterr

## -- not sure what this comment refers to: ##higher sd.fit increases correlation between fitness and trait...BUT reduces variation in fitness.

###pos neg correct in mutations???
#added abs needs rerun
mutate <- function(nL,N,sd.eff,mean.eff, prbmut) { sapply( 1:N , function(z)
			abs(rnorm(nL,mean=mean.eff,sd=sd.eff)) * ifelse(rbinom(nL,size=1,prob=0.5), 1 , -1 ) * rbinom(nL, size=1, prob = prbmut) ) #size = 1 in the last argument implies haploid?
			  }# DistofAbsValueofTraitEff * ProbofPosvNegMutation * ProbMutOccurs(u*loci) 
#provides individuals in columns, loci in rows
#prbmut is the per locus mutation rate


mutate.exp <- function(nL,N,lambda, prbmut) { sapply( 1:N , function(z)
			abs(rexp(nL,rate=lambda)) * ifelse(rbinom(nL,size=1,prob=0.5), 1 , -1 ) * rbinom(nL, size=1, prob = prbmut) ) #size = 1 in the last argument implies haploid?
			  }# DistofAbsValueofTraitEff * ProbofPosvNegMutation * ProbMutOccurs(u*loci) 
#high values of lambda give LOWER effect sizes of mutations on average.


horizontal <- function(nL,N,prb,genomat) {  #genomat is row are loci, cols are individuals
			 transfer <- matrix( rbinom(nL*N, size=1, prob = prb) * sample( 1:N , nL*N, replace =T) , ncol=N, nrow=nL)# with probability prb, sample a single locus from another individual, so rate of locus transfer proportional to prevalence in pop 
			 NewGenomat <- matrix( sapply(1:N, function(ind) sapply(1:nL, function(loc) 
			  				ifelse(transfer[loc,ind] > 0, genomat[loc, transfer[loc,ind] ] , genomat[loc, ind] ) 
			  				) ), ncol = ,  , byrow=F )
			  return(NewGenomat)
			  }
#requires and provides ind. in colums, loci in rows


sim.cotrait <- function(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,timesteps,Lambda,mutprb,fiterrP,fiterrM ,prbHorz, pfP, pfM,FLFC){  #-- 
#popsize, number loci, optimal phenotype, shallowness of fitness decline away from trait opt, timesteps, average effect of mutation, sd of mutation effect distribution, prb of mutation,  fitness error (i.e. from sources random with respect to trait), prb of horizontal transfer,freeliving fitness cost
#note plant and 
	Pa.mat <- array(0, dim= c(nlP,NP,2,timesteps+1) ) #store mutations on genotypes, rows have genotypes, columns have individuals, timesteps in sep matrices
		#extra dimension for diploidy, i.e. second genome in two individuals
	Ma.mat <- array(0, dim= c(nlM,NM,timesteps+1) ) 

	Pneut.mat <-  array(0, dim= c(nlnP,NP,2,timesteps+1) )
	Mneut.mat <-  array(0, dim= c(nlnM,NM,timesteps+1) )
	
	for(i in 2:(timesteps+1) ) {
		if(NP<NM) { 
			symbiotic <- sample(1:NM,NP, replace=F)
			freeliving <- !((1:NM)%in%symbiotic)
		} else if(NP==NM){
		symbiotic <- 1:NM
			freeliving <- !((1:NM)%in%symbiotic)
		}else{print("pop of microbes, NM, must be > NP, pop of plants")}	
		
		#interact! since Pa.mat[,,i-1] and Ma.mat[,,i-1] are randomly arranged from reproduction the previous time, just pair up columns / symbiotic columns
		#joint trait value is evaluated with respect to distance from host optima and microbe optima -- then relative fitness is calculated based on the distance to both own trait optima
		 ### and how well the partner can do based on its optima, with some proportion -- i.e. assumption is that if this is a mutualism, when your partner is unfit it isn't as beneficial
		hostfit <-     pfP*univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,symbiotic,i-1]), zopt = zoP, sd.fit=wP) + 
				   (1-pfP)*univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,symbiotic,i-1]), zopt = zoM, sd.fit=wM)
			#ASSUMPTION by using rowSums across the colSums, I assume additivity between alleles at the same site (no dominance)
		micrfit <- c()
		micrfit[symbiotic] 	<- (1-pfM)*univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,symbiotic,i-1]), zopt = zoP, sd.fit=wP) +  
					               pfM*univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,symbiotic,i-1]), zopt = zoM, sd.fit=wM)
		micrfit[freeliving] <- min(micrfit[symbiotic])*FLFC
		 ##VERY ARTIFICIAL FREELIVING FITNESS COST  FLFC.  THE STRONGER IT IS THE STRONGER THE HOST-INDUCED single generation bottleneck.
		 #if none are freeliving, this won't change micrfit vector.
		#trait effects additive across loci and across plant and microbe loci, thus everything just added together for the phenotype
		wH.e <- range01(hostfit +rnorm(length(hostfit),mean=0,sd=fiterrP)) ##range01 is required due to the addition of random error. It does not change the relative fitness (used later) of individuals whether you use range01 or add the minimum value (to avoid negative fitness)
		wM.e <- range01(micrfit +rnorm(length(micrfit),mean=0,sd=fiterrM))
				#Plant reproductions
		#reproduce: random mating with respect to relative fitness, and unlinked loci, but finite popsize.
		seedP <- sample( 1:NP , NP ,replace=T , prob=  (wH.e)/ sum((wH.e)) ) # 
		pollP <- sample( 1:NP , NP ,replace=T , prob=  (wH.e)/ sum((wH.e)) ) #selfing is allowed.
 		#reproduce proportional to fitness and add next generation mutations in pollen, then seed donor genomes, paste together, 	
		pollPs <- sapply(1:NP, function(n) 
						sapply(1:nlP, function(l) 
							 Pa.mat[ l, pollP[n], sample(c(1,2),size=1), i-1] )) +  
							 mutate.exp(nL=nlP,N=NP,lambda = Lambda,prbmut=mutprb)  
							      #pull out the the pollen parents, the recombined genomes that those parents will pass on, across loci & #add mutations
		seedPs <- sapply(1:NP, function(n) 
							sapply(1:nlP, function(l) 
								 Pa.mat[ l, seedP[n], sample(c(1,2),size=1), i-1] )) + #same for seeds
				 mutate.exp(nL=nlP,N=NP,lambda = Lambda,prbmut=mutprb) 
		pollneut <- sapply(1:NP, function(n)  #same for neutral loci
						sapply(1:nlnP, function(l) 
							 Pneut.mat[ l, pollP[n], sample(c(1,2),size=1), i-1] )) +  
					mutate.exp(nL=nlnP,N=NP,lambda = Lambda,prbmut=mutprb) 
		seedneut <- sapply(1:NP, function(n)  #same for neutral loci
						sapply(1:nlnP, function(l) 
							 Pneut.mat[ l, seedP[n], sample(c(1,2),size=1), i-1] )) +  
					mutate.exp(nL=nlnP,N=NP,lambda = Lambda,prbmut=mutprb) 

		Pa.mat[,,,i] <- abind(pollPs,seedPs,along=3) 
		Pneut.mat[,,,i] <- abind(pollneut,seedneut,along=3) 
		#microbial reproduction is clonal, have horizontal gene exchange after selection at probability prbHorz (per nL*N).
		micrPs <- sample( 1:NM , NM ,replace=T , prob= (wM.e)/ sum((wM.e)) )
		Ma.mat[,,i]    <- horizontal(nL = nlM, N=NM,prb=prbHorz, genomat = Ma.mat[,,i-1][,micrPs]  + mutate.exp(nL=nlM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 	
		Mneut.mat[,,i] <- horizontal(nL = nlnM, N=NM,prb=prbHorz, genomat = Mneut.mat[,,i-1][,micrPs]  + mutate.exp(nL=nlnM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 	
			
	}
	result <- list(Pa.mat,Ma.mat,Pneut.mat,Mneut.mat)
	names(result) <- c("Plant","Microbe","P_neutral","M_neutral")
	return(result)
}


windowplot <- function(first, last, thinsize, simdat,ylim,main) {
	twindows <- seq(from=first, to = last,by=thinsize)
 	m.p <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) )  )
 	m.mp <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
 	m.m <- sapply(twindows, function(t) mean( colSums(simdat$Microbe[,,t]) )  )
	r.p <- sapply(twindows, function(t) range( rowSums(colSums(simdat$Plant[,,,t]))  )  )
	r.mp <- sapply(twindows, function(t) range( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	r.m <- sapply(twindows, function(t) range( colSums(simdat$Microbe[,,t]) )  )
	yrange <- range(c(as.vector(r.p),as.vector(r.m)))
	plot(r.p[2,]~twindows,ylim=ylim,pch=NA,,xlab="",ylab="",main=main, cex.main=1)
	lines(m.p~twindows); lines(m.mp~twindows); lines(m.m~twindows)
	mtext("breeding values",side = 2, line=2)
	mtext("generations",side = 1, line=2)
	polygon( c(twindows,rev(twindows)), c(r.mp[1,],rev(r.mp[2,])),col=rgb(0,0,0,alpha=0.25))
	polygon( c(twindows,rev(twindows)), c(r.p[1,],rev(r.p[2,])),col=rgb(0,0,1,alpha=0.5))
	polygon( c(twindows,rev(twindows)), c(r.m[1,],rev(r.m[2,])),col=rgb(1,0,0,alpha=0.5))
}


#one question: when plant and microbes are unevenly sized in pop, what does this do?
getfitcon <- function(first, last, thinsize, simdat,zoP,zoM, wP, wM,pfP,pfM) {
	twindows <- seq(from=first, to = last,by=thinsize)
		cor.p <-sapply(twindows, function(t) cor(
			pfP*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoP, sd.fit=wP) + 
		         (1-pfP)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoM, sd.fit=wM) ,
		(1-pfM)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoP, sd.fit=wP) +  
                   pfM*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoM, sd.fit=wM)
		) )
		slp.p <- sapply(twindows, function(t)  glm(  
			(pfP*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoP, sd.fit=wP) + 
		         (1-pfP)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoM, sd.fit=wM)) ~
		I((1-pfM)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoP, sd.fit=wP) +  
                   pfM*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoM, sd.fit=wM))
		)$coef[2] ) 
 	return(list(fitnesscorrelation=cor.p,fitnessslope=slp.p))
}

plotfitcon <- function(fitconobj,twindows,ylim,main) {
	plot(fitconobj[[2]]~twindows,main=main,ylim=ylim,type="l",ylab="",xlab="",lty=2)
	par(new=T)
	plot(fitconobj[[1]]~twindows,main=main,ylim=c(-1,1),type="l",axes=F, xlab="", ylab="",lty=3)
	axis(side=4)
}

plotjfitcor <- function(fitconobj,twindows,ylim,main) {
	plot(fitconobj[[1]]~twindows,main=main,ylim=c(-1,1),type="l", xlab="", ylab="",lty=3)
}


gens <- 1000
test.4a <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,zoP=3,zoM=2,wP=0.75,wM=5,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)#at about .1 for fiterr relative to sdfit of 1 is when error in fitness starts to obscure fitness differences, not calculated. just apprx guess.	
test.4b <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,zoP=2,zoM=3,wP=5,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
test.4d <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,zoP=2,zoM=3,wP=1,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
test.4e <- sim.cotraitd(NP=100,NM=100,nlP=100,nlM=100,zoP=3,zoM=2,wP=0.75,wM=1,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)

test.q <- sim.cotrait(NP=10,NM=10,nlP=5,nlM=5,nlnP=10,nlnM=10,zoP=2,zoM=3,wP=1,wM=0.75,timesteps=20,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)

titles <- c("~0 link to microbe fitness","~0 link to plant fitness","microbe link > plant link","plant link > microbe link")

pdf("~/Dropbox/Simulation Results_PFF_expDFE.pdf",height=6,width=6)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(oma=c(1,0,1,0))
windowplot(1,1000,5,test.4a,ylim=c(-2,6),titles[1])
	abline(h=3,col=rgb(0,0,1),lty=2)
	abline(h=0)
windowplot(1,gens, 5,test.4b,ylim=c(-2,6),titles[2])
	abline(h=3,col=rgb(1,0,0),lty=2)
	abline(h=0)
	legend(gens*-0.05,y=6,c("Plant","Microbe","Holobiont"), fill=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,0,0,alpha=0.5)),  bty="n")
	legend(gens*0.35,y=6,c("Plant optima","Microbe optima","Midpoint"), lty = 2, col = c(rgb(0,0,1),rgb(1,0,0),rgb(0,0,0,alpha=.5)), bty="n")
windowplot(1,gens,5,test.4d,ylim=c(-2,6),titles[3])
	abline(h=2,col=rgb(0,0,1),lty=2); abline(h=3,col=rgb(1,0,0),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
windowplot.diploid(1,gens,5,test.4e,ylim=c(-2,6),titles[4])
	abline(h=3,col=rgb(0,0,1),lty=2); abline(h=2,col=rgb(1,0,0),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
dev.off()

fitcon.4a <- getfitcon(2,gens,5,test.4a,zoP=3,zoM=2,wP=0.75,wM=5,pfM=0.6,pfP=0.6)
fitcon.4b <- getfitcon(2,gens,5,test.4b,zoP=2,zoM=3,wP=5,wM=0.75,pfM=0.6,pfP=0.6)
fitcon.4d <- getfitcon(2,gens,5,test.4d,zoP=2,zoM=3,wP=1,wM=0.75,pfM=0.6,pfP=0.6)
fitcon.4e <- getfitcon(2,gens,5,test.4e,zoP=3,zoM=2,wP=0.75,wM=1,pfM=0.6,pfP=0.6)

pdf("~/Dropbox/Simulation Results fitcor PFF expDFE.pdf",height=6,width=6)
par(mfrow=c(2,2))
par(mar=c(3,3,3,3))
par(oma=c(1,0,1,0))
plotjfitcor(fitcon.4a,seq(from=1,to=gens,by=5),ylim=0.1*c(-10,10),main=titles[1])
	mtext("fitness correlation",side = 2, line=2,adj=-2)
	abline(h=0,lty=1)
plotjfitcor(fitcon.4b,seq(from=1,to=gens,by=5),ylim=0.1*c(-10,10),main=titles[2])
	abline(h=0,lty=1)
plotjfitcor(fitcon.4d,seq(from=1,to=gens,by=5),ylim=0.1*c(-10,10),main=titles[3])
	abline(h=0,lty=1)
plotjfitcor(fitcon.4e,seq(from=1,to=gens,by=5),ylim=0.1*c(-10,10),main=titles[4])
	abline(h=0,lty=1)
	mtext("generations",side = 1, line=2, adj=-0.75)
dev.off()



###simulate multiple populations; run experiment on result; create genotypes. 

#list parameters; keep everything the same except optima?
plant.opt <- rep(c(1:5))

popset <- lapply(1:length(plant.opt), function(E)
	 sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,nlnP=100,nlnM=100,zoP=E,zoM=E,wP=1,wM=1,timesteps=200,Lambda=30,mutprb=0.0001,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
)




getalleles <- function(sel.finalT.P, sel.finalT.M){ #final timepoint genotypes as produced within function run.exp()

	nLP <- dim(sel.finalT.P[[1]])[1] 
	nLM <- dim(sel.finalT.M[[1]])[1]
	NP <- dim(sel.finalT.P[[1]])[2]
	NM <- dim(sel.finalT.M[[1]])[2]
	 ## the dimensions we're concerned with. NOTE: this is likely different from # loci and # individuals in sims, but abbrv the same
	uniquevalsP <- lapply(1:nLP, function(L)  sort(unique(unlist(lapply(1:length(sel.finalT.P), function(POP) as.vector(sel.finalT.P[[POP]][L,,]))  ) )  ))
	uniquevalsM <- lapply(1:nLM, function(L)  sort(unique(unlist(lapply(1:length(sel.finalT.M), function(POP) as.vector(sel.finalT.M[[POP]][L,]))  ) )  ))
	lociinblocksP <- unlist(lapply(uniquevalsP, function(z) length(z)-1 )) #number of new loci required
	lociinblocksM <- unlist(lapply(uniquevalsM, function(z) length(z)-1 )) #number of new loci required
	#if(FALSE) "yes" else "no"
	newlocinumP <- lapply(1:length(lociinblocksP), function(z) if(lociinblocksP[z]>0) c(1,1:lociinblocksP[z]) else 1  )
	newgenoP <- lapply(1:length(uniquevalsP), function(z)  c(0, rep(1,times=lociinblocksP[z]))   )
	newlocinumM <- lapply(1:length(lociinblocksM), function(z) if(lociinblocksM[z]>0) c(1,1:lociinblocksM[z]) else 1  )
	newgenoM <- lapply(1:length(uniquevalsM), function(z)  c(0, rep(1,times=lociinblocksM[z]))   )
	fixedincludedP <- lociinblocksP
	fixedincludedP[fixedincludedP==0] <- 1
	fixedincludedM <- lociinblocksM
	fixedincludedM[fixedincludedM==0] <- 1
	
	newgenomatP <- matrix(0,nrow=sum(fixedincludedP),ncol=NP*length(sel.finalT.P)*2) 	#I think it makes sense to make the default 0, each individual can have at most 2 non-zero loci in one linkage block
	colsallele1 <- seq(from=1, to = NP*length(sel.finalT.P)*2, by =2)
	colsallele2 <- seq(from=2, to = NP*length(sel.finalT.P)*2, by =2)
	
	for(l in 1:nLP){
		prevgenos1 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,,1]) ))# individuals are rows, alleles in columns
		prevgenos2 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,,2]) ))# individuals are rows, alleles in columns
#		prevgenos cbind(prevgenos1,prevgenos2)
		whichval1 <- sapply(1:length(prevgenos1), function(a) which(uniquevalsP[[l]]==prevgenos1[a]) )
		whichval2 <- sapply(1:length(prevgenos2), function(a) which(uniquevalsP[[l]]==prevgenos2[a]) )
		newloc1 <- sapply(1:length(prevgenos1), function(a)  newlocinumP[[l]][whichval1[a]]  )
		newloc2 <- sapply(1:length(prevgenos2), function(a)  newlocinumP[[l]][whichval2[a]]  )
		newallele1 <- sapply(1:length(prevgenos1), function(a)  newgenoP[[l]][whichval1[a]]  )
		newallele2 <- sapply(1:length(prevgenos2), function(a)  newgenoP[[l]][whichval2[a]]  )
		startrow <- if(l > 1) sum(fixedincludedP[1:(l-1)])+1 else 1
		for(i in 1:(NP*length(sel.finalT.P)) ){
			newgenomatP[ (startrow -1 + newloc1[i]) , colsallele1[i]  ] <- newallele1[i]
			newgenomatP[ (startrow -1 + newloc2[i]) , colsallele2[i]  ] <- newallele2[i]
		} 
	}

	newgenomatM <- matrix(0,nrow=sum(fixedincludedM),ncol=NM*length(sel.finalT.M)) 	#I think it makes sense to make the default 0, each individual can have at most 2 non-zero loci in one linkage block
	for(l in 1:nLM){
		prevgenos <- unlist(lapply(1:length(sel.finalT.M), function(POP) (sel.finalT.M[[POP]][l,]) ))# individuals are rows, alleles in columns
		whichval <- sapply(1:length(prevgenos), function(a) which(uniquevalsM[[l]]==prevgenos[a]) )
		newloc <- sapply(1:length(prevgenos), function(a)  newlocinumM[[l]][whichval[a]]  )
		newallele <- sapply(1:length(prevgenos), function(a)  newgenoM[[l]][whichval[a]]  )
		startrow <- if(l > 1) sum(fixedincludedM[1:(l-1)])+1 else 1
		for(i in 1:(NM*length(sel.finalT.M)) ){
			newgenomatM[ (startrow -1 + newloc[i]) , i  ] <- newallele[i]
		} 
	}

#add information about 
#add column to newgenomatP for linkage group; and also one for allele effect, and "all zero" effect -- e.g. because the first locus has 2 values.	
return(list(genoP = newgenomatP, genoM=newgenomatM) )
}


run.exp<- function(popsetobj,numperpop,exp.err){
		#SELECT INDIVIDUALS
		dimP <- dim(popsetobj[[1]]$Plant)
		dimM <- dim(popsetobj[[1]]$Microbe)	
		p.geno.s <- unlist( lapply(1:length(popsetobj), function(E)  sample(1:(dimP[2]),numperpop, repl=F)))
		p.pop.s <- rep(1:length(popsetobj), each = numperpop)
		m.geno.s <- unlist(lapply(1:length(popsetobj), function(E)  sample(1:(dimM[2]),numperpop, repl=F) ) )
		m.pop.s <- rep(1:length(popsetobj), each = numperpop)

		#PHENOTYPES			
		p.breed <- sapply(1:length(p.geno.s), function(G) sum(   popsetobj[[ p.pop.s[G] ]]$Plant[,p.geno.s[G],,dimP[4]]  ) )
		m.breed <- sapply(1:length(m.geno.s), function(G) sum(   popsetobj[[ m.pop.s[G] ]]$Microbe[,m.geno.s[G],dimP[4]] ) )
		#assume full factorial
		GbyGt <- unlist(sapply(1:length(p.breed), function(P) sapply(1:length(m.breed), function(M) 
					 p.breed[P] + m.breed[M] +	 rnorm(1,mean=0,sd=exp.err) ) )) #each Pbreed is a column in a matrix of rows = # mbreed, when unlist, first col, then second, then...
		GbyGpgrec <- rep(p.geno.s, each = length(m.geno.s)) 
		GbyGpprec <- rep(p.pop.s,  each = length(m.geno.s))
		GbyGmgrec <- rep(m.geno.s, times = length(p.geno.s))
		GbyGmprec <- rep(m.pop.s,  times = length(p.geno.s))
		dat <- data.frame(traitvalue = as.vector(sapply(1:length(p.breed), function(P) sapply(1:length(m.breed), function(M) 
					 p.breed[P] + m.breed[M] +	 rnorm(1,mean=0,sd=exp.err) ) )), #each Pbreed is a column in a matrix of rows = # mbreed, when unlist, first col, then second, then...
						  genoP = rep(p.geno.s, each = length(m.geno.s)) ,
						  popP = rep(p.pop.s,  each = length(m.geno.s)),
						  genoM = rep(m.geno.s, times = length(p.geno.s)),
						  popM = rep(m.pop.s,  times = length(p.geno.s)))

		sel.finalT.P <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$Plant[,p.geno.s[p.pop.s==POP],,dimP[4]]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.M  <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$Microbe[,m.geno.s[p.pop.s==POP],dimP[4]]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)
		sel.finalT.Pn <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$P_neutral[,p.geno.s[p.pop.s==POP],,dimP[4]]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.Mn  <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$M_neutral[,m.geno.s[p.pop.s==POP],dimP[4]]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)

		causalgenos <- getalleles(sel.finalT.P,sel.finalT.M)
		neutralgenos <- getalleles(sel.finalT.Pn,sel.finalT.Mn)
		ind.datP <- data.frame(
			ID = 1:length
		)
		ind.datPM = list(plant = ind.datP, microbe= ind.datM)
		return(list(expdat = dat, causalgenos = causalgenos, neutralgenos = neutralgenos, ind.dat = ind.datPM))
}



#export genotypes and phenotypes to GEMMA
###how to get genotypes for gemma? this is a weird case with (depending on loci, time and # pops) possibly many unique allele values per locus

###OPTION A
##for each trait locus in each organism, figure out how many unique values there are.  that apprx is the number of total mutuations; accurately is the number of unique alleles
## to convert to biallelic SNPs, needs to be the number of unique values , ceiling.
## then need to assign each a SNP value (literally just the quantity assigned by the locus) and location. perfectly linked within locus and unlinked across

getalleles <- function(sel.finalT.P, sel.finalT.M){ #final timepoint genotypes as produced within function run.exp()

	nLP <- dim(sel.finalT.P[[1]])[1] 
	nLM <- dim(sel.finalT.M[[1]])[1]
	NP <- dim(sel.finalT.P[[1]])[2]
	NM <- dim(sel.finalT.M[[1]])[2]
	 ## the dimensions we're concerned with. NOTE: this is likely different from # loci and # individuals in sims, but abbrv the same
	uniquevalsP <- lapply(1:nLP, function(L)  sort(unique(unlist(lapply(1:length(sel.finalT.P), function(POP) as.vector(sel.finalT.P[[POP]][L,,]))  ) )  ))
	uniquevalsM <- lapply(1:nLM, function(L)  sort(unique(unlist(lapply(1:length(sel.finalT.M), function(POP) as.vector(sel.finalT.M[[POP]][L,]))  ) )  ))
	lociinblocksP <- unlist(lapply(uniquevalsP, function(z) length(z)-1 )) #number of new loci required
	lociinblocksM <- unlist(lapply(uniquevalsM, function(z) length(z)-1 )) #number of new loci required
	#if(FALSE) "yes" else "no"
	newlocinumP <- lapply(1:length(lociinblocksP), function(z) if(lociinblocksP[z]>0) c(1,1:lociinblocksP[z]) else 1  )
	newgenoP <- lapply(1:length(uniquevalsP), function(z)  c(0, rep(1,times=lociinblocksP[z]))   )
	newlocinumM <- lapply(1:length(lociinblocksM), function(z) if(lociinblocksM[z]>0) c(1,1:lociinblocksM[z]) else 1  )
	newgenoM <- lapply(1:length(uniquevalsM), function(z)  c(0, rep(1,times=lociinblocksM[z]))   )
	fixedincludedP <- lociinblocksP
	fixedincludedP[fixedincludedP==0] <- 1
	fixedincludedM <- lociinblocksM
	fixedincludedM[fixedincludedM==0] <- 1
	
	newgenomatP <- matrix(0,nrow=sum(fixedincludedP),ncol=NP*length(sel.finalT.P)*2) 	#I think it makes sense to make the default 0, each individual can have at most 2 non-zero loci in one linkage block
	colsallele1 <- seq(from=1, to = NP*length(sel.finalT.P)*2, by =2)
	colsallele2 <- seq(from=2, to = NP*length(sel.finalT.P)*2, by =2)
	
	for(l in 1:nLP){
		prevgenos1 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,,1]) ))# individuals are rows, alleles in columns
		prevgenos2 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,,2]) ))# individuals are rows, alleles in columns
#		prevgenos cbind(prevgenos1,prevgenos2)
		whichval1 <- sapply(1:length(prevgenos1), function(a) which(uniquevalsP[[l]]==prevgenos1[a]) )
		whichval2 <- sapply(1:length(prevgenos2), function(a) which(uniquevalsP[[l]]==prevgenos2[a]) )
		newloc1 <- sapply(1:length(prevgenos1), function(a)  newlocinumP[[l]][whichval1[a]]  )
		newloc2 <- sapply(1:length(prevgenos2), function(a)  newlocinumP[[l]][whichval2[a]]  )
		newallele1 <- sapply(1:length(prevgenos1), function(a)  newgenoP[[l]][whichval1[a]]  )
		newallele2 <- sapply(1:length(prevgenos2), function(a)  newgenoP[[l]][whichval2[a]]  )
		startrow <- if(l > 1) sum(fixedincludedP[1:(l-1)])+1 else 1
		for(i in 1:(NP*length(sel.finalT.P)) ){
			newgenomatP[ (startrow -1 + newloc1[i]) , colsallele1[i]  ] <- newallele1[i]
			newgenomatP[ (startrow -1 + newloc2[i]) , colsallele2[i]  ] <- newallele2[i]
		} 
	}

	newgenomatM <- matrix(0,nrow=sum(fixedincludedM),ncol=NM*length(sel.finalT.M)) 	#I think it makes sense to make the default 0, each individual can have at most 2 non-zero loci in one linkage block
	for(l in 1:nLM){
		prevgenos <- unlist(lapply(1:length(sel.finalT.M), function(POP) (sel.finalT.M[[POP]][l,]) ))# individuals are rows, alleles in columns
		whichval <- sapply(1:length(prevgenos), function(a) which(uniquevalsM[[l]]==prevgenos[a]) )
		newloc <- sapply(1:length(prevgenos), function(a)  newlocinumM[[l]][whichval[a]]  )
		newallele <- sapply(1:length(prevgenos), function(a)  newgenoM[[l]][whichval[a]]  )
		startrow <- if(l > 1) sum(fixedincludedM[1:(l-1)])+1 else 1
		for(i in 1:(NM*length(sel.finalT.M)) ){
			newgenomatM[ (startrow -1 + newloc[i]) , i  ] <- newallele[i]
		} 
	}

	locusdatP <- data.frame(zerostate = unlist(    lapply(uniquevalsP, function(L) c(L[1], if(length(L)>2) rep(0, times = length(L)-2) else NULL  )  )  ), 
							reststate = unlist(    lapply(uniquevalsP, function(L)   if(length(L)>1) L[-1] else 0 ) )   ,
							linkage = unlist(    lapply(1:length(uniquevalsP), function(L)   rep(L, times = if(length(uniquevalsP[[L]])>1) length(uniquevalsP[[L]])-1 else 1 ) ) )   )
							locusdatP$location <-  unlist(sapply(1:length(unique(locusdatP$linkage)), function(z) seq( from = 0, to = 100, length.out = table(locusdatP$linkage)[z]) ))
							###LOCATION IS VERY FORCED AT THIS POINT. 
	locusdatM <- data.frame(zerostate = unlist(    lapply(uniquevalsM, function(L) c(L[1], if(length(L)>2) rep(0, times = length(L)-2) else NULL  )  )  ), 
							reststate = unlist(    lapply(uniquevalsM, function(L)   if(length(L)>1) L[-1] else 0 ) )   ,
							linkage = unlist(    lapply(1:length(uniquevalsM), function(L)   rep(L, times = if(length(uniquevalsM[[L]])>1) length(uniquevalsM[[L]])-1 else 1 ) ) )   )
							locusdatM$location <-  unlist(sapply(1:length(unique(locusdatM$linkage)), function(z) seq( from = 0, to = 100, length.out = table(locusdatM$linkage)[z]) ))
#add information about 
#add column to newgenomatP for linkage group; and also one for allele effect, and "all zero" effect -- e.g. because the first locus has 2 values.	
return(list(genoP = newgenomatP, genoM=newgenomatM, locusdatP = locusdatP, locusdatM = locusdatM) )
}

# 		for(nl in unique(newloc1) ){}
# 			newgenomatP[ c(startrow -1 + nl) , colsallele1  ] <- newallele1
# 		}
# 		for(nl in unique(newloc2) ){}
# 			newgenomatP[ c(startrow -1 + nl) , colsallele2 ] <- newallele2
# 		}
# 		


expset <- run.exp(popset, numperpop= 3,exp.err=0.05)




###OPTION B?
#alternatively. what if there are just linked and unlinked markers? -- i.e. a matrix where the same columns, rows (and allele if diploid) are drawn, but states are gatc
##but then this would have to have extra sites so there can be enough alleles. and how is it really different than assigning them after the fact?


##how to enact space in a way that makes sense.

##EXPORT TO GEMMA