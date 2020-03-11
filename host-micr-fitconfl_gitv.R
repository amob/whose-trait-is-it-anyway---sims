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
			##yes, but new questions about what to do with the extra microbes. 
	#diploid plant stage under selection, rather than haploid 
	#some proportion of focal fitness comes from partner fitness, and some percentage from having the right phenotype, with error on top
	#exponenetial distribution of effect sizes of mutations? #ALLOWED.

##flexibility for further consideration:
	## Dominance?
	#Mapping issue: keep track of alleles? -- maybe I can just use unique effect size values within loci? although there might be similar ways to arrive at the same effect sizes....
		# could say that each mutational change in an allele's history is a new site (inf sites) but that all within a locus are non-recombining?
	#more than one microbe species?
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

install.packages("abind")
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

mutate.expROUND2 <- function(nL,N,lambda, prbmut) { sapply( 1:N , function(z)
			round(rexp(nL,rate=lambda),digits=2) * ifelse(rbinom(nL,size=1,prob=0.5), 1 , -1 ) * rbinom(nL, size=1, prob = prbmut) ) #size = 1 in the last argument implies haploid?
			  }# DistofAbsValueofTraitEff * ProbofPosvNegMutation * ProbMutOccurs(u*loci) 
##this limits unique alleles a little bit.
# not the same as geneflow across population
##why use? might be a way to solve the problem of population structure = causal allele structure.
#we could apply it only to the neutral matrix?

horizontal <- function(nL,N,prb,genomat) {  #genomat is row are loci, cols are individuals
			 transfer <- matrix( rbinom(nL*N, size=1, prob = prb) * sample( 1:N , nL*N, replace =T) , ncol=N, nrow=nL)# with probability prb, sample a single locus from another individual, so rate of locus transfer proportional to prevalence in pop 
			 NewGenomat <- matrix( sapply(1:N, function(ind) sapply(1:nL, function(loc) 
			  				ifelse(transfer[loc,ind] > 0, genomat[loc, transfer[loc,ind] ] , genomat[loc, ind] ) 
			  				) ), ncol = ,  , byrow=F )
			  return(NewGenomat)
			  }
#requires and provides ind. in colums, loci in rows


sim.cotrait <- function(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,timesteps,Lambda,mutprb,fiterrP,fiterrM ,prbHorz, pfP, pfM,FLFC,startmats = "n"){  #-- 
#popsize, number loci, optimal phenotype, shallowness of fitness decline away from trait opt, timesteps, average effect of mutation, sd of mutation effect distribution, prb of mutation,  fitness error (i.e. from sources random with respect to trait), prb of horizontal transfer,freeliving fitness cost
#startmats must be list of matrices NAMED: Pa.mat, Ma.mat, Pneut.mat, Mneut.mat
#note plant and
 
	#initialize
	Pa.mat <- array(0, dim= c(nlP,NP,2,timesteps+1) ) #store mutations on genotypes, rows have genotypes, columns have individuals, timesteps in sep matrices
		#extra dimension for diploidy, i.e. second genome in two individuals
	Ma.mat <- array(0, dim= c(nlM,NM,timesteps+1) ) 

	Pneut.mat <-  array(0, dim= c(nlnP,NP,2,timesteps+1) )
	Mneut.mat <-  array(0, dim= c(nlnM,NM,timesteps+1) )
	
	#allow initialize with previous state
	if(startmats != "n"){
		Pa.mat[,,1,1] <- startmat$Pa.mat[,,1]
		Pa.mat[,,2,1] <- startmat$Pa.mat[,,2]
		Pneut.mat[,,1,1] <- startmat$Pneut.mat[,,1]
		Pneut.mat[,,2,1] <- startmat$Pneut.mat[,,2]
		Ma.mat[,,1] <- startmat$Ma.mat[,]
		Mneut.mat[,,1] <- startmat$Mneut.mat[,]
	} else{"assume empty starting matrices"} 
	#End initialize
	
	for(i in 2:(timesteps+1) ) {
		if(NP<NM) { 
			symbiotic <- 1:NP # should be no reason to further randomize 
			freeliving <- !((1:NM)%in%symbiotic) #this may not be the easiest solution.
		} else if(NP==NM){
		symbiotic <- 1:NM
			freeliving <- !((1:NM)%in%symbiotic)
		}else{print("pop of microbes, NM, must be >= NP, pop of plants")}	
		
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
							 mutate.expROUND2(nL=nlP,N=NP,lambda = Lambda,prbmut=mutprb)  
							      #pull out the the pollen parents, the recombined genomes that those parents will pass on, across loci & #add mutations
		seedPs <- sapply(1:NP, function(n) 
							sapply(1:nlP, function(l) 
								 Pa.mat[ l, seedP[n], sample(c(1,2),size=1), i-1] )) + #same for seeds
				 mutate.expROUND2(nL=nlP,N=NP,lambda = Lambda,prbmut=mutprb) 
		pollneut <- sapply(1:NP, function(n)  #same for neutral loci
						sapply(1:nlnP, function(l) 
							 Pneut.mat[ l, pollP[n], sample(c(1,2),size=1), i-1] )) +  
					mutate.expROUND2(nL=nlnP,N=NP,lambda = Lambda,prbmut=mutprb) 
		seedneut <- sapply(1:NP, function(n)  #same for neutral loci
						sapply(1:nlnP, function(l) 
							 Pneut.mat[ l, seedP[n], sample(c(1,2),size=1), i-1] )) +  
					mutate.expROUND2(nL=nlnP,N=NP,lambda = Lambda,prbmut=mutprb) 

		Pa.mat[,,,i] <- abind(pollPs,seedPs,along=3) 
		Pneut.mat[,,,i] <- abind(pollneut,seedneut,along=3) 
		#microbial reproduction is clonal, have horizontal gene exchange after selection at probability prbHorz (per nL*N).
		micrPs <- sample( 1:NM , NM ,replace=T , prob= (wM.e)/ sum((wM.e)) )
		Ma.mat[,,i]    <- horizontal(nL = nlM, N=NM,prb=prbHorz, genomat = Ma.mat[,,i-1][,micrPs]  + mutate.expROUND2(nL=nlM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 	
		Mneut.mat[,,i] <- horizontal(nL = nlnM, N=NM,prb=prbHorz, genomat = Mneut.mat[,,i-1][,micrPs]  + mutate.expROUND2(nL=nlnM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 	
			
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
		pfit <- sapply(twindows, function(t) 	pfP*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoP, sd.fit=wP) + 
		         (1-pfP)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoM, sd.fit=wM) )#twindows are columns
		mfit <- sapply(twindows, function(t) (1-pfM)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoP, sd.fit=wP) +  
                   pfM*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]), zopt = zoM, sd.fit=wM) )
		cor.p <-sapply(1:ncol(pfit), function(t) cor(pfit[,t]/sum(pfit[,t]),mfit[,t]/sum(mfit[,t])))
		slp.p <- sapply(1:ncol(pfit), function(t)  glm(pfit[,t] ~ mfit[,t])$coef[2] ) 
 	return(list(fitnesscorrelation=cor.p,fitnessslope=slp.p))
} #this function ASSUMES NO ERROR IN FITNESS, since it cannot be applied exactly as was simulated. alternative option is to record observed fitnesses.
#partners are the same as in the simulation, however, because interact function just pairs the columns of genotypes. 

plotfitcon <- function(fitconobj,twindows,main,ylim="n") {
	if(ylim == "n") {ylim= c(min(c(tqFC[[1]],tqFC[[2]])),max(c(tqFC[[1]],tqFC[[2]])))} else{ylim= ylim}
	plot(fitconobj[[2]]~twindows,main=main,ylim=ylim,type="l",ylab="",xlab="",lty=2)
	par(new=T)
	plot(fitconobj[[1]]~twindows,main=main,ylim=c(-1,1),type="l",axes=F, xlab="", ylab="",lty=3)
	axis(side=4)

} 

#how much conflict is there?
plotjfitcor <- function(fitconobj,twindows,ylim,main) {
	plot(fitconobj[[1]]~twindows,main=main,ylim=c(-1,1),type="l", xlab="", ylab="",lty=3)
}
#who is winning?
extractwinning <- function(simdat,first,last,eachNth,zoM,zoP){
	twindows <- seq(from=first, to = last,by=eachNth)
	z <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	dP <- abs(zoM-z)
	dM <- abs(zoP-z)
	return(list(dP=dP,dM=dM))
}

#which genome has more var?
extractVmVp <- function(simdat,first,last,eachNth){
	twindows <- seq(from=first, to = last,by=eachNth)
	Vp <- sapply(twindows, function(t) var( rowSums(colSums(simdat$Plant[,,,t])) )  )
	Vm <- sapply(twindows, function(t) var( colSums(simdat$Microbe[,,t]) )  )
	Vb <- sapply(twindows, function(t) var( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	return(list(Vp=Vp, Vm = Vm,Vb=Vb,PVp=Vp / (Vp+Vm), PVm = Vm/(Vp+Vm)))
}

#how do the following parameters change ans to above:
	#strength of link between trait and fitness; wP and wM, lower is stronger link -- with lower links, PFF becomes more important
		wP.v <- c(0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
		wM.v <- c(0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75?
	#drift and selection within each species; # uhhh. wP and wM vs fiterrP and fiterrM
		#N and s
		fiterrP.v <- c(0.0001, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		fiterrM.v <- c(0.0001, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		popsz.v <- c(10, 25, 50, 75, 100, 150, 200, 500) #note turning up M without N is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
	# mutational inputs;  Lambda
		Lambda.v <- seq(from = 40, to = 20, by =-2) #base 30
	#symmetries in conflict between partners (are these poss? I think just change pff for one only)
		pfP.v <- seq(from = 0.05, to =0.95, by =0.1) #base 0.6	


#short-term questions
#equilibrium?
extractDyn <- function(simdat,first,last,windowsize){
	twindows <- seq(from=first, to = last,by=windowsize)
	tcoefP <- sapply(twindows, function(tw) 
			glm(sapply( (tw):(tw+windowsize-1), function(tstp)  mean(rowSums(colSums(simdat$Plant[,,,tstp])) )  )~c(1:windowsize))$coef[2] )
	 
	tcoefM <- sapply(twindows, function(tw) 
			glm(sapply( (tw):(tw+windowsize-1), function(tstp)  mean( colSums(simdat$Microbe[,,tstp])) )  ~c(1:windowsize))$coef[2] )
	return(list(tcoefP = tcoefP, tcoefM=tcoefM ))
}
	
		
#e.g. super conserved important traits  will  have  very  little  variation.   Also  those  under  strong  (recent) directional selection



gens <- 500
test.4a <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,zoP=3,zoM=2,wP=0.75,wM=5,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)#at about .1 for fiterr relative to sdfit of 1 is when error in fitness starts to obscure fitness differences, not calculated. just apprx guess.	
test.4b <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,zoP=2,zoM=3,wP=5,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
test.4d <- sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,zoP=2,zoM=3,wP=1,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
test.4e <- sim.cotraitd(NP=100,NM=100,nlP=100,nlM=100,zoP=3,zoM=2,wP=0.75,wM=1,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)


# tq.pmu <- sim.cotrait(NP=100,NM=100,nlP=50,nlM=100,nlnP=3,nlnM=3,zoP=-2,zoM=-3,wP=1,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
# tqnFC<- getfitcon(10, gens, 1, tq.pmu,zoP=-2,zoM=-3, wP=1, wM=0.75,pfP=0.6,pfM=0.6) 
# tqnVmVp <- extractVmVp(tq.pmu, 1,gens,1)
# tqnwin <- extractwinning(tq.pmu,first=1,last=gens,1,zoP=-2,zoM= -3)
# tqndyn <- extractDyn(tq.pmu,first=1,last=gens,5)
# par(mfrow=c(2,3))
# windowplot(1,gens,1,tq.pmu,ylim=c(-5,0),main="tqpmu")
# plotjfitcor(tqnFC,10:gens,main="tqn",ylim=c(-1,1)); abline(h=0)
# plot(tqnVmVp$Vp ~ c(1:length(tqnVmVp$Vp)),pch=NA,ylim=c(0,max(c(tqnVmVp$Vp,tqnVmVp$Vm)))); lines( 1:length(tqnVmVp$Vp),tqnVmVp$Vp,col=rgb(0,0,1)) ; lines(1:length(tqnVmVp$Vm),tqnVmVp$Vm , col=rgb(1,0,0)) 
# plot(tqnVmVp$PVp ~ c(1:length(tqnVmVp$Vp)),pch=NA,ylim=c(0,1)); lines( 1:length(tqnVmVp$Vp),tqnVmVp$PVp,col=rgb(0,0,1)) ; lines(1:length(tqnVmVp$Vm),tqnVmVp$PVm , col=rgb(1,0,0)) 
# plot(tqnwin$dP ~ c(1:length(tqnwin$dP)),pch=NA,ylim=c(0,max(c(tqnwin$dP,tqnwin$dP)))); lines(tqnwin$dP~c(1:gens),col=rgb(0,0,1)); lines(tqnwin$dM~c(1:gens),col=rgb(1,0,0))
# plot(tqndyn$tcoefP ~ c(1:length(tqndyn$tcoefP)),pch=NA,ylim=c(min(c(tqndyn$tcoefP,tqndyn$tcoefM)),max(c(tqndyn$tcoefP,tqndyn$tcoefM)))); lines(tqndyn$tcoefP~c(1:length(tqndyn$tcoefP)),col=rgb(0,0,1)); lines(tqndyn$tcoefM~c(1:length(tqndyn$tcoefM)),col=rgb(1,0,0))
# 
# test.q <- sim.cotrait(NP=100,NM=100,nlP=50,nlM=50,nlnP=10,nlnM=10,zoP=2,zoM=3,wP=1,wM=0.75,timesteps=gens,Lambda=30,mutprb=0.0005,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
# 
# tqFC<- getfitcon(1, gens, 1, test.q,zoP=2,zoM=3, wP=1, wM=0.75,pfP=0.6,pfM=0.6) 
# tqVmVp <- extractVmVp(test.q, 1,gens,1)
# tqwin <- extractwinning(test.q,first=1,last=gens,1,zoP=2,zoM= 3)
# par(mfrow=c(2,3))
# windowplot(1,gens,1,test.q,ylim=c(0,3),main="tq")
# plotjfitcor(tqFC,1:gens,main="tq",ylim=c(-1,1)); abline(h=0)
# plot(tqVmVp$Vp ~ c(1:length(tqVmVp$Vp)),pch=NA,ylim=c(0,max(c(tqVmVp$Vp,tqVmVp$Vm)))); lines( 1:length(tqVmVp$Vp),tqVmVp$Vp,col=rgb(0,0,1)) ; lines(1:length(tqVmVp$Vm),tqVmVp$Vm , col=rgb(1,0,0)) 
# plot(tqVmVp$PVp ~ c(1:length(tqVmVp$Vp)),pch=NA,ylim=c(0,1)); lines( 1:length(tqVmVp$Vp),tqVmVp$PVp,col=rgb(0,0,1)) ; lines(1:length(tqVmVp$Vm),tqVmVp$PVm , col=rgb(1,0,0)) 
# plot(tqwin$dP ~ c(1:length(tqwin$dP)),pch=NA,ylim=c(0,max(c(tqwin$dP,tqwin$dP)))); lines(tqwin$dP~c(1:gens),col=rgb(0,0,1)); lines(tqwin$dM~c(1:gens),col=rgb(1,0,0))

titles <- c("~0 link to microbe fitness","~0 link to plant fitness","microbe link > plant link","plant link > microbe link")

pdf("~/Simulation Results_PFF_expDFE.pdf",height=6,width=6)
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
windowplot(1,gens,5,test.4e,ylim=c(-2,6),titles[4])
	abline(h=3,col=rgb(0,0,1),lty=2); abline(h=2,col=rgb(1,0,0),lty=2); abline(h=2.5,col=rgb(0,0,0,alpha=.5),lty=2) 
	abline(h=0)
dev.off()

fitcon.4a <- getfitcon(2,gens,5,test.4a,zoP=3,zoM=2,wP=0.75,wM=5,pfM=0.6,pfP=0.6)
fitcon.4b <- getfitcon(2,gens,5,test.4b,zoP=2,zoM=3,wP=5,wM=0.75,pfM=0.6,pfP=0.6)
fitcon.4d <- getfitcon(2,gens,5,test.4d,zoP=2,zoM=3,wP=1,wM=0.75,pfM=0.6,pfP=0.6)
fitcon.4e <- getfitcon(2,gens,5,test.4e,zoP=3,zoM=2,wP=0.75,wM=1,pfM=0.6,pfP=0.6)

pdf("~/Simulation Results fitcor PFF expDFE.pdf",height=6,width=6)
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



##EXPERIMENT FUNCTION


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
	# 	GbyGpgrec <- rep(p.geno.s, each = length(m.geno.s)) 
# 		GbyGpprec <- rep(p.pop.s,  each = length(m.geno.s))
# 		GbyGmgrec <- rep(m.geno.s, times = length(p.geno.s))
# 		GbyGmprec <- rep(m.pop.s,  times = length(p.geno.s))
		dat <- data.frame(traitvalue = as.vector(sapply(1:length(p.breed), function(P) sapply(1:length(m.breed), function(M) 
					 p.breed[P] + m.breed[M] +	 rnorm(1,mean=0,sd=exp.err) ) )), #each Pbreed is a column in a matrix of rows = # mbreed, when unlist, first col, then second, then...
						  genoP = rep(p.geno.s, each = length(m.geno.s)) ,
						  popP = rep(p.pop.s,  each = length(m.geno.s)),
						  genoM = rep(m.geno.s, times = length(p.geno.s)),
						  popM = rep(m.pop.s,  times = length(p.geno.s)),
						  IDPG = rep(1:length(p.geno.s), each = length(m.geno.s)),
						  IDMG = rep(1:length(m.geno.s), times = length(p.geno.s)) )

		sel.finalT.P <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$Plant[,p.geno.s[p.pop.s==POP],,dimP[4]]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.M  <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$Microbe[,m.geno.s[p.pop.s==POP],dimP[4]]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)
		sel.finalT.Pn <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$P_neutral[,p.geno.s[p.pop.s==POP],,dimP[4]]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.Mn  <- lapply(1:length(popsetobj), function(POP)  popsetobj[[POP]]$M_neutral[,m.geno.s[p.pop.s==POP],dimP[4]]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)

		causalgenos <- getalleles(sel.finalT.P,sel.finalT.M,numperpop)
		neutralgenos <- getalleles(sel.finalT.Pn,sel.finalT.Mn,numperpop)
		ind.datP <- data.frame(ID = 1:length(p.geno.s), POP = p.pop.s,	withinPOPgeno = p.geno.s)
		ind.datM <- data.frame(ID = 1:length(m.geno.s), POP = m.pop.s,	withinPOPgeno = m.geno.s)
		ind.datPM = list(plant = ind.datP, microbe= ind.datM)
		return(list(expdat = dat, causalgenos = causalgenos, neutralgenos = neutralgenos, ind.dat = ind.datPM))
}



##A FUNCTION TO GET GENOTYPES TO USE IN GEMMA (called in experiment function above)

###how to get genotypes for gemma? this is a weird case with (depending on loci, time and # pops) possibly many unique allele values per locus
###OPTION A
##for each trait locus in each organism, figure out how many unique values there are.  that apprx is the number of total mutuations; accurately is the number of unique alleles
## to convert to biallelic SNPs, needs to be the number of unique values , ceiling.
## then need to assign each a SNP value (literally just the quantity assigned by the locus) and location. perfectly linked within locus and unlinked across
###OPTION B?
#alternatively. what if there are just linked and unlinked markers? -- i.e. a matrix where the same columns, rows (and allele if diploid) are drawn, but states are gatc
##but then this would have to have extra sites so there can be enough alleles. and how is it really different than assigning them after the fact?
####RUNNING WITH OPTION A

##how to throw markers into genome space in a way that makes sense - just some random splitting. shouldn't affect GWAS, so this is trivial.
#what to do abou pop structure??  -- semi solution with rounding above


###THERE IS A PROBLEM WHEN THE LOWEST ALLELIC EFFECT IS NEGATIVE AS OPPOSED TO 0, I think. results in some pops appearing to share identical mutations...
#seems not not make bonus category when effects are negative.
#because I relied on 0 being first in code. 
getalleles <- function(sel.finalT.P, sel.finalT.M,numperpop){ #final timepoint genotypes as produced within function run.exp()
	if(numperpop>1){
		nLP <- dim(sel.finalT.P[[1]])[1] 
		nLM <- dim(sel.finalT.M[[1]])[1]
		NP <- dim(sel.finalT.P[[1]])[2]
		NM <- dim(sel.finalT.M[[1]])[2]
		 ## the dimensions we're concerned with. NOTE: this is likely different from # loci and # individuals in sims, but abbrv the same
		uniquevalsP <- lapply(1:nLP, function(L)  unique(unlist(lapply(1:length(sel.finalT.P), function(POP) as.vector(sel.finalT.P[[POP]][L,,]))  ) )  )
		uniquevalsM <- lapply(1:nLM, function(L)  unique(unlist(lapply(1:length(sel.finalT.M), function(POP) as.vector(sel.finalT.M[[POP]][L,]))  ) )  )
	}
	else {
		nLP <- dim(sel.finalT.P[[1]])[1] 
		nLM <- length(sel.finalT.M[[1]])[1] 
		NP <- 1
		NM <- 1
		uniquevalsP <- lapply(1:nLP, function(L)  unique(unlist(lapply(1:length(sel.finalT.P), function(POP) as.vector(sel.finalT.P[[POP]][L,]))  ) )  )
		uniquevalsM <- lapply(1:nLM, function(L)  unique(unlist(lapply(1:length(sel.finalT.M), function(POP) as.vector(sel.finalT.M[[POP]][L]))  ) )  )
	}
	uniquevalsP.s <- lapply(uniquevalsP, function(z)  z[order(abs(z))])
	uniquevalsM.s <- lapply(uniquevalsM, function(z)  z[order(abs(z))])
	nozeroP <- unlist(lapply(uniquevalsP, function(z) length(which(z == 0))!=1 ))
	nozeroM <- unlist(lapply(uniquevalsM, function(z) length(which(z == 0))!=1 ))

	lociinblocksP <- unlist(lapply(uniquevalsP.s, function(z) length(z) -1  )) + sapply(nozeroP, function(z) if(z) 1 else 0) #number of loci required  
	# THE PROBLEM IS: if the ancestral "0" isn't first, than plants have to be either of two derived alleles according to this plan, which is not good.
	# If this is 0 is still present, resorting the loci helps. if it isn't there, then you need new locus for every allele.
	lociinblocksM <- unlist(lapply(uniquevalsM.s, function(z) length(z) -1  )) + sapply(nozeroM, function(z) if(z) 1 else 0)  #number of loci required

	#if(FALSE) "yes" else "no"
	newlocinumP <- lapply(1:length(lociinblocksP), function(z) if(lociinblocksP[z]==0) 1 else if(nozeroP[z]) c(1:lociinblocksP[z]) else c(1,1:lociinblocksP[z])  )
	newgenoP <- lapply(1:length(uniquevalsP.s), function(z) if(nozeroP[z]) rep(1,times=lociinblocksP[z]) else c(0, rep(1,times=lociinblocksP[z]))    )
#	newlocinumM <- lapply(1:length(lociinblocksM), function(z) if(lociinblocksM[z]>0) c(1,1:lociinblocksM[z]) else 1  )
	newlocinumM <- lapply(1:length(lociinblocksM), function(z) if(lociinblocksM[z]==0) 1 else if(nozeroM[z]) c(1:lociinblocksM[z]) else c(1,1:lociinblocksM[z])  )
# 	newgenoM <- lapply(1:length(uniquevalsM), function(z)  c(0, rep(1,times=lociinblocksM[z]))   )
	newgenoM <- lapply(1:length(uniquevalsM.s), function(z) if(nozeroM[z]) rep(1,times=lociinblocksM[z]) else c(0, rep(1,times=lociinblocksM[z]))    )
	fixedincludedP <- lociinblocksP
	fixedincludedP[fixedincludedP==0] <- 1
	fixedincludedM <- lociinblocksM
	fixedincludedM[fixedincludedM==0] <- 1
	
	newgenomatP <- matrix(0,nrow=sum(fixedincludedP),ncol=NP*length(sel.finalT.P)*2) 	#I think it makes sense to make the default 0, each individual can have at most 2 non-zero loci in one linkage block
	colsallele1 <- seq(from=1, to = NP*length(sel.finalT.P)*2, by =2)
	colsallele2 <- seq(from=2, to = NP*length(sel.finalT.P)*2, by =2)
	
	for(l in 1:nLP){
		if(numperpop>1){
		prevgenos1 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,,1]) ))# individuals are rows, alleles in columns
		prevgenos2 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,,2]) ))# individuals are rows, alleles in columns
		} else {
		prevgenos1 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,1]) ))# individuals are rows, alleles in columns
		prevgenos2 <- unlist(lapply(1:length(sel.finalT.P), function(POP) (sel.finalT.P[[POP]][l,2]) ))# individuals are rows, alleles in columns
		}
#		prevgenos cbind(prevgenos1,prevgenos2)
		whichval1 <- sapply(1:length(prevgenos1), function(a) which(uniquevalsP.s[[l]]==prevgenos1[a]) )
		whichval2 <- sapply(1:length(prevgenos2), function(a) which(uniquevalsP.s[[l]]==prevgenos2[a]) )
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
		if(numperpop>1){
		prevgenos <- unlist(lapply(1:length(sel.finalT.M), function(POP) (sel.finalT.M[[POP]][l,]) ))# individuals are rows, alleles in columns
		} else {
		prevgenos <- unlist(lapply(1:length(sel.finalT.M), function(POP) (sel.finalT.M[[POP]][l]) ))# individuals are rows, alleles in columns
		}
		whichval <- sapply(1:length(prevgenos), function(a) which(uniquevalsM.s[[l]]==prevgenos[a]) )
		newloc <- sapply(1:length(prevgenos), function(a)  newlocinumM[[l]][whichval[a]]  )
		newallele <- sapply(1:length(prevgenos), function(a)  newgenoM[[l]][whichval[a]]  )
		startrow <- if(l > 1) sum(fixedincludedM[1:(l-1)])+1 else 1
		for(i in 1:(NM*length(sel.finalT.M)) ){
			newgenomatM[ (startrow -1 + newloc[i]) , i  ] <- newallele[i]
		} 
	}


	locusdatP <- data.frame(zerostate = unlist(    lapply(1:length(uniquevalsP.s), function(L)
											 if(nozeroP[L]) rep(0, times=length(uniquevalsP.s[[L]])) else (  c( uniquevalsP.s[[L]][1], if(length(uniquevalsP.s[[L]])>2) rep(0, times = length(uniquevalsP.s[[L]])-2) else NULL )  )  
											 )  ), 
							reststate = unlist(    lapply(1:length(uniquevalsP.s), function(L)  if(nozeroP[L]) uniquevalsP.s[[L]] else (if(length(uniquevalsP.s[[L]])>1) uniquevalsP.s[[L]][-1] else 0) ) )   ,
							linkage = unlist(    lapply(1:length(uniquevalsP.s), function(L)    rep(L, times = if(nozeroP[L]) length(uniquevalsP.s[[L]]) else if(length(uniquevalsP.s[[L]])>1) length(uniquevalsP.s[[L]])-1 else 1 ) ) )   )
							locusdatP$location <-  round(unlist(sapply(1:length(unique(locusdatP$linkage)), function(z) seq( from = 0, to = 100, length.out = table(locusdatP$linkage)[z]) )))
							###LOCATION IS VERY FORCED AT THIS POINT. 
	locusdatM <- data.frame(zerostate = unlist(    lapply(1:length(uniquevalsM.s), function(L)
											 if(nozeroM[L]) rep(0, times=length(uniquevalsM.s[[L]])) else (  c( uniquevalsM.s[[L]][1], if(length(uniquevalsM.s[[L]])>2) rep(0, times = length(uniquevalsM.s[[L]])-2) else NULL )  )  
												 )  ), 
							reststate = unlist(    lapply(1:length(uniquevalsM.s), function(L)  if(nozeroM[L]) uniquevalsM.s[[L]] else (if(length(uniquevalsM.s[[L]])>1) uniquevalsM.s[[L]][-1] else 0) ) )   ,
							linkage = unlist(    lapply(1:length(uniquevalsM.s), function(L)    rep(L, times = if(nozeroM[L]) length(uniquevalsM.s[[L]]) else if(length(uniquevalsM.s[[L]])>1) length(uniquevalsM.s[[L]])-1 else 1 ) ) )   )
							locusdatM$location <-  round(unlist(sapply(1:length(unique(locusdatM$linkage)), function(z) seq( from = 0, to = 100, length.out = table(locusdatM$linkage)[z]) )))
#add information about 
#add column to newgenomatP for linkage group; and also one for allele effect, and "all zero" effect -- e.g. because the first locus has 2 values.	
return(list(genoP = newgenomatP, genoM=newgenomatM, locusdatP = locusdatP, locusdatM = locusdatM) )
}


####GENERATE OUTPUT FILES FOR SIMULATIONS WITH GWAS
###simulate multiple populations; run experiment on result; create genotypes. 

#list parameters; keep everything the same except optima?
plant.opt <- rep(c(5:15))

popset <- lapply(1:length(plant.opt), function(E)
	 sim.cotrait(NP=100,NM=100,nlP=100,nlM=100,nlnP=100,nlnM=100,zoP=E,zoM=E,wP=1,wM=1,timesteps=200,Lambda=30,mutprb=0.0001,fiterrP=0.001,fiterrM=0.001,prbHorz = 0.2,pfP = 0.6, pfM=0.6,FLFC=0.1)
)

expset <- run.exp(popset, numperpop= 1,exp.err=0.05)

#how should the GWAS be run exaxtly?
#genotypes and phenotypes together, with only one phenotype column so genotypes must be repeated
#presumably that means 1 file each for phenotypes with plant genotypes and phenotypes with microbe genotypes
# -- OR we could paste together genos of HOSTS and MICROBES A-LA experimental design?
#trying this second option for now.

##MAKE OUTPUT FILES FOR EXPORT TO GEMMA
#GEMMA needs plink format for genotype/phenotype information
##MAP PED FOR EACH.

#PED
#tab or space delimited
#individuals are rows and loci are columns
#first columns are
#    Family ID
#      Individual ID
#      Paternal ID
#      Maternal ID
#      Sex (1=male; 2=female; other=unknown)
#      Phenotype
#then columns 7 onward MUST be paired alleles for each locus
#alleles must be in 2 column format even for haploids. alleles are 1 and 2 not 0 and 1
#read-in options allow you to specify that there is missing info for family, mother, father, and sex; but do need Individual ID and Phenotype
# --no-fid  --no-parents  --no-sex --no-pheno 
#but this is the same as having them missing (set to 0) in the first place
#NO header row
# it needs to be reordered so that each individual is a row, but each locus still has two columns
RGP <- t( sapply(seq(from=1, to =ncol(expset$causalgenos$genoP),by=2), function(Ind) 
				c(t(  rbind(expset$causalgenos$genoP,expset$neutralgenos$genoP)[,c(Ind,Ind+1)] )) ) )  #concatenating plant genotype matrix, merging individual columns so that loci rather than individuals are repeated, repeating across loci and transposing to rowxcol = individuals x loci (2 cols per locus)

HOLOgeno <- cbind( paste("p",expset$expdat$IDPG,"m",expset$expdat$IDMG,sep=""),  #Family ID when missing, plink wants this= to individual ID
		paste("p",expset$expdat$IDPG,"m",expset$expdat$IDMG,sep=""), # Individual ID
		rep(0,times=nrow(expset$expdat)) , rep(0,times=nrow(expset$expdat)) ,rep(0,times=nrow(expset$expdat)) , # 
		expset$expdat$traitvalue,
	    (RGP+1) [expset$expdat$IDPG,] , # adding 1 so bt 1 and 2, and then grabbing rows as they are used in the experiment.
 			#then
 	   t(rbind( expset$causalgenos$genoM[rep(1:nrow(expset$causalgenos$genoM)  ,each=2), ], 
 	   			## duplicating rows for pretend diploid microbe genotype matrix (lociin2rowsxindividuals), 
 	   		   expset$neutralgenos$genoM[rep(1:nrow(expset$neutralgenos$genoM) ,each=2), ])+1) [expset$expdat$IDMG,]  
 	   		   #concatenating neutral loci treated the same way, transposing to ind x loci(2 cols per), adding 1 so bt 1 and 2, 
 	   		   #and then grabbing individuals in rows as they are used in the experiment.
 	   )
PLANTgeno <- cbind( paste("p",expset$expdat$IDPG,"m",expset$expdat$IDMG,sep=""), 
		paste("p",expset$expdat$IDPG,"m",expset$expdat$IDMG,sep=""),
		rep(0,times=nrow(expset$expdat)) , rep(0,times=nrow(expset$expdat)) ,rep(0,times=nrow(expset$expdat)) ,
		expset$expdat$traitvalue,
	    RGP [expset$expdat$IDPG,]  # adding 1 so bt 1 and 2, and then grabbing rows as they are used in the experiment.
  	   )
MICRgeno <- cbind( paste("p",expset$expdat$IDPG,"m",expset$expdat$IDMG,sep=""), 
		paste("p",expset$expdat$IDPG,"m",expset$expdat$IDMG,sep=""),
		rep(0,times=nrow(expset$expdat)) , rep(0,times=nrow(expset$expdat)) ,rep(0,times=nrow(expset$expdat)) ,
		expset$expdat$traitvalue,
	 	t(rbind( expset$causalgenos$genoM[rep(1:nrow(expset$causalgenos$genoM)  ,each=2), ], 
 	   			## duplicating rows for pretend diploid microbe genotype matrix (lociin2rowsxindividuals), 
 	   		   expset$neutralgenos$genoM[rep(1:nrow(expset$neutralgenos$genoM) ,each=2), ])+1) [expset$expdat$IDMG,]  
 	   		   #concatenating neutral loci treated the same way, transposing to ind x loci(2 cols per), adding 1 so bt 1 and 2, 
 	   		   #and then grabbing individuals in rows as they are used in the experiment.
 	   )
write.table(HOLOgeno,"~/HOLOevosims.ped",quote=F,sep="\t",row.names=F,col.names=F)
write.table(PLANTgeno,"~/PLANTevosims.ped",quote=F,sep="\t",row.names=F,col.names=F)
write.table(MICRgeno,"~/MICRevosims.ped",quote=F,sep="\t",row.names=F,col.names=F)

#MAP
#The markers in the PED file do not need to be in genomic order: 
#******BUT the order MAP file should align with the order of the PED file markers# ******
# By default, each line of the MAP file describes a single marker and must contain exactly 4 columns:
#NO COLUMN NAMES, I think???

HOLOmap <- rbind( cbind( paste("scaffold_",c(expset$causalgenos$locusdatP$linkage,expset$neutralgenos$locusdatP$linkage),sep=""),  #      chromosome (1-22, X, Y or 0 if unplaced)
#linkage group must be alphanumeric not numeric if > 22 chromosomes. so paste "scaffold" first "--allow-extra-chr" flag to make plink accept them.
			c(paste("cP",1:nrow(expset$causalgenos$locusdatP),sep=""),paste("nP",1:nrow(expset$neutralgenos$locusdatP),sep="") ),#      rs# or snp identifier
	 		rep(0,times=nrow(expset$causalgenos$locusdatP)+ nrow(expset$neutralgenos$locusdatP)), #      Genetic distance (morgans) # For basic association testing, the genetic distance column can be set at 0. 
	 		c(expset$causalgenos$locusdatP$location,expset$neutralgenos$locusdatP$location) ), #      Base-pair position (bp units)
	 		#
	   		cbind( paste("scaffold_",c(expset$causalgenos$locusdatM$linkage,expset$neutralgenos$locusdatM$linkage),sep=""),  #      chromosome (1-22, X, Y or 0 if unplaced)
#linkage group must be alphanumeric not numeric if > 22 chromosomes. so paste "scaffold" first "--allow-extra-chr" flag to make plink accept them.
		 	c(paste("cM",1:nrow(expset$causalgenos$locusdatM),sep=""),paste("nM",1:nrow(expset$neutralgenos$locusdatM),sep="") ),#      rs# or snp identifier
		 	rep(0,times=nrow(expset$causalgenos$locusdatM)+ nrow(expset$neutralgenos$locusdatM)), #      Genetic distance (morgans) # For basic association testing, the genetic distance column can be set at 0. 
		 	c(expset$causalgenos$locusdatM$location,expset$neutralgenos$locusdatM$location) ) #      Base-pair position (bp units)
			)
PLANTmap <-  cbind( paste("scaffold_",c(expset$causalgenos$locusdatP$linkage,expset$neutralgenos$locusdatP$linkage),sep=""),  #      chromosome (1-22, X, Y or 0 if unplaced)
#linkage group must be alphanumeric not numeric if > 22 chromosomes. so paste "scaffold" first "--allow-extra-chr" flag to make plink accept them.
			c(paste("cP",1:nrow(expset$causalgenos$locusdatP),sep=""),paste("nP",1:nrow(expset$neutralgenos$locusdatP),sep="") ),#      rs# or snp identifier
	 		rep(0,times=nrow(expset$causalgenos$locusdatP)+ nrow(expset$neutralgenos$locusdatP)), #      Genetic distance (morgans) # For basic association testing, the genetic distance column can be set at 0. 
	 		c(expset$causalgenos$locusdatP$location,expset$neutralgenos$locusdatP$location) ) #      Base-pair position (bp units)
MICRmap <- cbind( paste("scaffold_",c(expset$causalgenos$locusdatM$linkage,expset$neutralgenos$locusdatM$linkage),sep=""),  #      chromosome (1-22, X, Y or 0 if unplaced)
#linkage group must be alphanumeric not numeric if > 22 chromosomes. so paste "scaffold" first "--allow-extra-chr" flag to make plink accept them.
		 	c(paste("cM",1:nrow(expset$causalgenos$locusdatM),sep=""),paste("nM",1:nrow(expset$neutralgenos$locusdatM),sep="") ),#      rs# or snp identifier
		 	rep(0,times=nrow(expset$causalgenos$locusdatM)+ nrow(expset$neutralgenos$locusdatM)), #      Genetic distance (morgans) # For basic association testing, the genetic distance column can be set at 0. 
		 	c(expset$causalgenos$locusdatM$location,expset$neutralgenos$locusdatM$location) ) #      Base-pair position (bp units)
write.table(HOLOmap,"~/HOLOevosims.map",quote=F,sep="\t",row.names=F,col.names=F)
write.table(PLANTmap,"~/PLANTevosims.map",quote=F,sep="\t",row.names=F,col.names=F)
write.table(MICRmap,"~/MICRevosims.map",quote=F,sep="\t",row.names=F,col.names=F)
  
#GEMMA can include SNP covariate information, e.g. ; the file could include whether SNP is host or microbial; whether it is causal or not
#this is stored in SNP names for now
#GEMMA can include individual/phenotype covariate information.
#might use if go with different way of dealing with 2 genome issue.


##PASS .ped and .map to PLINK to get binary formats .bim .fam .bed
##PASS to GEMMA to run models.