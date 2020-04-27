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


#must install on scratch node first. can't be done from within script, or write over existing install?
#install.packages("abind")
library(abind)

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

###The current code.

#fitness function for trait, zopt
#based on normal distribution, but instead of pdf summing to one, fitness is 1 at the optimum, just remove the 1/sqrt(2*pi*sd) in the normal P.D.F., not sure who to cite for deterministic part of this formula, I feel like I have seen it a lot...
univar.fit <- function(z, zopt, sd.fit) {   
										rawfit <-   exp( -1* ((z-zopt)^2) / (2 * (sd.fit^2)) )  
										adjfit <- rawfit /sum(rawfit)
										return(adjfit)
										} 
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


mutate.exp <- function(nL,N,lambda, prbmut) { sapply( 1:N , function(z)
			abs(rexp(nL,rate=lambda)) * ifelse(rbinom(nL,size=1,prob=0.5), 1 , -1 ) * rbinom(nL, size=1, prob = prbmut) ) #size = 1 in the last argument implies haploid?
			  }# DistofAbsValueofTraitEff * ProbofPosvNegMutation * ProbMutOccurs(u*loci) 
#high values of lambda give LOWER effect sizes of mutations on average.


#This version of horizontal transmission one is like picking up random fragments from the environment
horizontal <- function(nL,N,prb,genomat) {  #genomat is row are loci, cols are individuals
			 transfer <- matrix( rbinom(nL*N, size=1, prob = prb) * sample( 1:N , nL*N, replace =T) , ncol=N, nrow=nL)# with probability prb, sample a single locus from another individual, so rate of locus transfer proportional to prevalence in pop 
			 NewGenomat <- matrix( sapply(1:N, function(ind) sapply(1:nL, function(loc) 
			  				ifelse(transfer[loc,ind] > 0, genomat[loc, transfer[loc,ind] ] , genomat[loc, ind] ) 
			  				) ), ncol = N,  , byrow=F )
			  return(NewGenomat)
			  }
#requires and provides ind. in colums, loci in rows

#this one is more like conjugation. where some are randomly on plasmids			
# horizontal <- function(nL,N,prb,genomat) {  #genomat is row are loci, cols are individuals
# 			 transfer <- matrix( rbinom(nL*N, size=1, prob = prb) * rep(sample( 1:N , N, replace =T),each=nL ), ncol=N, nrow=nL,byrow=F)# with probability prb, sample a single locus from another individual, so rate of locus transfer proportional to prevalence in pop 
# 			 NewGenomat <- matrix( sapply(1:N, function(ind) sapply(1:nL, function(loc) 
# 			  				ifelse(transfer[loc,ind] > 0, genomat[loc, transfer[loc,ind] ] , genomat[loc, ind] ) 
# 			  				) ), ncol = N,  , byrow=F )
# 			  return(NewGenomat)
# 			  }



sim.cotrait <- function(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,timesteps,Lambda,mutprb ,prbHorz, pfP, pfM,FLFC,startmats = "n",zoptvects = "n"){  #-- 
#popsize, number loci, optimal phenotype, shallowness of fitness decline away from trait opt, timesteps, average effect of mutation, 
	#exponential rate of mutation effect distribution, prb of mutation, prb of horizontal transfer,freeliving fitness cost
	#startmats, if not "n", must be a list of matrices NAMED: Pa.mat, Ma.mat, Pneut.mat, Mneut.mat
	#zoptvects, if not "n", must be a list of vectors NAMED: PlantZ and MicrZ
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
	} else{print("assuming empty starting matrices")} 
	#End initialize
	
	for(i in 2:(timesteps+1) ) {
		if(length(zoptvects)==1){ #shortcut. there are too many warnings with better logical statement:  & zoptvects=="n"){
				zoP = zoP
				zoM = zoM} else{
				zoP = zoptvects$PlantZ[i-1]
				zoM = zoptvects$MicrZ[i-1]
				}
# 		if(NP<NM) #{ print("unequal pop sizes currently not allowed")
# 			symbiotic <- 1:NP # should be no reason to further randomize 
# 			freeliving <- !((1:NM)%in%symbiotic) #this may not be the easiest solution.
# 		} else if(NP==NM){
# 		symbiotic <- 1:NM
# 			freeliving <- !((1:NM)%in%symbiotic)
# 		} #else{print("pop of microbes, NM, must be >= NP, pop of plants")}	
		
		#interact! since Pa.mat[,,i-1] and Ma.mat[,,i-1] are randomly arranged from reproduction the previous time, just pair up columns / symbiotic columns
		#joint trait value is evaluated with respect to distance from host optima and microbe optima -- then relative fitness is calculated based on the distance to both own trait optima
		 ### and how well the partner can do based on its optima, with some proportion -- i.e. assumption is that if this is a mutualism, when your partner is unfit it isn't as beneficial
		hostfit <-     univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,,i-1]), zopt = zoP, sd.fit=wP) #, fiterr=fiterrP) 
		micrfit <-     univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,,i-1]), zopt = zoM, sd.fit=wM)# fiterr=fiterrM)#, freeliving = freeliving,FLCL=FLCL) 
		wH.e <-     pfP*hostfit + (1-pfP)*micrfit
		wM.e <- (1-pfM)*hostfit +     pfM*micrfit
	##NOTE: if wanted to do antagonistic interactions, do pfP*hostfit + (1-pfP)*((1-micrfit)/sum(1-micrfit))  #e.g. direct impact of trait on fitness, and impact passed through the inverse (e.g. perfectly negatively correlated) of impact on antagonists fitness
		#Plant reproductions
		#reproduce: random mating with respect to relative fitness, and unlinked loci, but finite popsize.
		seedP <- sample( 1:NP , NP ,replace=T , prob=  wH.e) # (wH.e)/ sum((wH.e)) ) # 
		pollP <- sample( 1:NP , NP ,replace=T , prob=  wH.e)#(wH.e)/ sum((wH.e)) ) #selfing is allowed.
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
		micrPs <- sample( 1:NM , NM ,replace=T , prob= wM.e)#(wM.e)/ sum((wM.e)) )
		Ma.mat[,,i]    <- horizontal(nL = nlM, N=NM,prb=prbHorz, genomat = Ma.mat[,,i-1][,micrPs]  + mutate.exp(nL=nlM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 	
		Mneut.mat[,,i] <- horizontal(nL = nlnM, N=NM,prb=prbHorz, genomat = Mneut.mat[,,i-1][,micrPs]  + mutate.exp(nL=nlnM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 			
	}
	result <- list(Pa.mat,Ma.mat,Pneut.mat,Mneut.mat)
	names(result) <- c("Plant","Microbe","P_neutral","M_neutral")
	return(result)
}


windowplot <- function(first, last, thinsize, simdat,ylim,main,ylabs="breeding values",xlabs="generations") {
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
	mtext(ylabs,side = 2, line=2)
	mtext(xlabs,side = 1, line=2)
	polygon( c(twindows,rev(twindows)), c(r.mp[1,],rev(r.mp[2,])),col=rgb(0,0,0,alpha=0.25),border=NA)
	polygon( c(twindows,rev(twindows)), c(r.p[1,],rev(r.p[2,])),col=rgb(0,0.5,0,alpha=0.5),border=NA)
	polygon( c(twindows,rev(twindows)), c(r.m[1,],rev(r.m[2,])),col=rgb(0.5,0,0.5,alpha=0.5),border=NA)
}



#one question: when plant and microbes are unevenly sized in pop, what does this do?
#now it discards unpartnered microbes. before, unclear what it did
getfitcon <- function(first, last, thinsize, simdat,zoP,zoM, wP, wM,pfP,pfM) {
	twindows <- seq(from=first, to = last,by=thinsize)
	if(dim(simdat$Plant)[2]<(dim(simdat$Microbe)[2])){
		mcols <- 1:(dim(simdat$Plant)[2])
	} else {mcols <- 1:(dim(simdat$Microbe)[2])}
		pfit <- sapply(twindows, function(t) 	pfP*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,mcols,t]), zopt = zoP, sd.fit=wP) + 
		         (1-pfP)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,mcols,t]), zopt = zoM, sd.fit=wM) )#twindows are columns
		mfit <- sapply(twindows, function(t) (1-pfM)*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,mcols,t]), zopt = zoP, sd.fit=wP) +  
                   pfM*univar.fit(rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,mcols,t]), zopt = zoM, sd.fit=wM) )
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
extractwinning <- function(simdat,first,last,eachNth,zoM,zoP,zoptvects="n"){ #same requirements for zoptvects as in sim.cotrait()
	twindows <- seq(from=first, to = last,by=eachNth)
	z <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	if(length(zoptvects)>1){
		dP <-abs(c(0,zoptvects$PlantZ)-z)
		dM <-abs(c(0,zoptvects$MicrZ)-z)
	} else{
		dP <- abs(zoP-z)
		dM <- abs(zoM-z)}
	return(list(dP=dP,dM=dM))
}

#which genome has more var?
extractVmVp <- function(simdat,first,last,eachNth){
	twindows <- seq(from=first, to = last,by=eachNth)
	Vp <- sapply(twindows, function(t) var( rowSums(colSums(simdat$Plant[,,,t])) )  )
	Vm <- sapply(twindows, function(t) var( colSums(simdat$Microbe[,,t]) )  )
	Vb <- sapply(twindows, function(t) var( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	return(list(Vp=Vp, Vm = Vm,Vb=Vb,PVp=Vp / (Vp+Vm), PVm = Vm/(Vp+Vm)))
}#currently pVx is a ratio of each to the sum, but not to the breeding value variance.

#short-term questions
#equilibrium?
extractDyn <- function(simdat,first,last,windowsize){
	twindows <- seq(from=first+1, to = last,by=windowsize)
	tcoefP <- sapply(twindows, function(tw) 
			glm(sapply( (tw):(tw+windowsize-1), function(tstp)  mean(rowSums(colSums(simdat$Plant[,,,tstp])) )  )~c(1:windowsize))$coef[2] )
	 
	tcoefM <- sapply(twindows, function(tw) 
			glm(sapply( (tw):(tw+windowsize-1), function(tstp)  mean( colSums(simdat$Microbe[,,tstp])) )  ~c(1:windowsize))$coef[2] )
	return(list(tcoefP = tcoefP, tcoefM=tcoefM ))
}

#unanswered q?
#e.g. super conserved important traits  will  have  very  little  variation.   Also  those  under  strong  (recent) directional selection

