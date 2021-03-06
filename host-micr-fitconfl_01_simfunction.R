#######
##This script is the simulation function, and simulation processing functions, these function are called and used in many other scripts
##Contains code for procedures described in "Simulation Details" of manuscript
#######
library(abind)

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}


#fitness function for trait, zopt
#based on normal distribution, but instead of pdf summing to one, fitness is 1 at the optimum, just remove the 1/sqrt(2*pi*sd) in the normal P.D.F., 
univar.fit <- function(z, zopt, sd.fit) {   #note that sd.fit is "omega" when this function is used in the simulation
										rawfit <-   exp( -1* ((z-zopt)^2) / (2 * (sd.fit^2)) )  
										adjfit <- rawfit /sum(rawfit) #relativise depending on other trait values in population
										return(adjfit)
										} 
#see lecorre and kramer 2012, for example, but from Haldane 1954 (full citations in manuscript)


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


sim.cotrait <- function(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,timesteps,Lambda,mutprb ,prbHorz, pfP, pfM,FLFC,startmats = "n",zoptvects = "n"){  #-- 
#popsize, number loci, optimal phenotype, shallowness of fitness decline away from trait opt (this is omega, note that w in code is represents either omega or fitness),
    # timesteps, average effect of mutation, 
	#exponential rate of mutation effect distribution, prb of mutation, prb of horizontal transfer,freeliving fitness cost
	#startmats, if not "n", must be a list of matrices NAMED: Pa.mat, Ma.mat, Pneut.mat, Mneut.mat
	#zoptvects, if not "n", must be a list of vectors NAMED: PlantZ and MicrZ
 
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
		if(length(zoptvects)==1){ 
				zoP = zoP
				zoM = zoM} else{
				zoP = zoptvects$PlantZ[i-1]
				zoM = zoptvects$MicrZ[i-1]
				}
		
		#interact! since Pa.mat[,,i-1] and Ma.mat[,,i-1] are randomly arranged from reproduction the previous time, just pair up columns / symbiotic columns
		#joint trait value is evaluated with respect to distance from host optima and microbe optima -- then relative fitness is calculated based on the distance to both own trait optima
		 ### and how well the partner can do based on its optima, with some proportion -- i.e. assumption is that if this is a mutualism, when your partner is unfit it isn't as beneficial
		hostfit <-     univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,,i-1]), zopt = zoP, sd.fit=wP) #individual trait-based fitness components
		micrfit <-     univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,,i-1]), zopt = zoM, sd.fit=wM) #individual trait-based fitness components
		wH.e <-     pfP*hostfit + (1-pfP)*micrfit #relative fitness
		wM.e <- (1-pfM)*hostfit +     pfM*micrfit #relative fitness
		#Plant reproductions
		#reproduce: random mating with respect to relative fitness, and unlinked loci, but finite popsize.
		seedP <- sample( 1:NP , NP ,replace=T , prob=  wH.e) # 
		pollP <- sample( 1:NP , NP ,replace=T , prob=  wH.e)# #selfing is allowed.
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
		micrPs <- sample( 1:NM , NM ,replace=T , prob= wM.e)
		Ma.mat[,,i]    <- horizontal(nL = nlM, N=NM,prb=prbHorz, genomat = Ma.mat[,,i-1][,micrPs]  + mutate.exp(nL=nlM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 	
		Mneut.mat[,,i] <- horizontal(nL = nlnM, N=NM,prb=prbHorz, genomat = Mneut.mat[,,i-1][,micrPs]  + mutate.exp(nL=nlnM,N=NM,lambda=Lambda,prbmut=mutprb) ) # 			
	}
	result <- list(Pa.mat,Ma.mat,Pneut.mat,Mneut.mat)
	names(result) <- c("Plant","Microbe","P_neutral","M_neutral")
	return(result)
}


####Functions for describing the results of a run of the simulation function above

windowplot <- function(first, last, thinsize, simdat,ylim,main,ylabs="breeding values",xlabs="generations") {
	twindows <- seq(from=first, to = last,by=thinsize)
 	m.p <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) )  )
 	m.mp <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
 	m.m <- sapply(twindows, function(t) mean( colSums(simdat$Microbe[,,t]) )  )
	r.p <- sapply(twindows, function(t) range( rowSums(colSums(simdat$Plant[,,,t]))  )  )
	r.mp <- sapply(twindows, function(t) range( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	r.m <- sapply(twindows, function(t) range( colSums(simdat$Microbe[,,t]) )  )
	sd.p <- sapply(twindows, function(t) sd( rowSums(colSums(simdat$Plant[,,,t]))  )  )
	sd.mp <- sapply(twindows, function(t) sd( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	sd.m <- sapply(twindows, function(t) sd( colSums(simdat$Microbe[,,t]) )  )
	sdw.p <- sapply(1:length(twindows), function(t) c( m.p[t] + sd.p[t], m.p[t] - sd.p[t] )  )
	sdw.mp <- sapply(1:length(twindows), function(t) c( m.mp[t] + sd.mp[t], m.mp[t] - sd.mp[t] )  )
	sdw.m <- sapply(1:length(twindows), function(t) c( m.m[t] + sd.m[t], m.m[t] - sd.m[t] )  )
	yrange <- range(c(as.vector(r.p),as.vector(r.m)))
	plot(r.p[2,]~twindows,ylim=ylim,pch=NA,,xlab="",ylab="",main=main, cex.main=1)
	lines(m.p~twindows); lines(m.mp~twindows); lines(m.m~twindows)
	mtext(ylabs,side = 2, line=2)
	mtext(xlabs,side = 1, line=2)
	polygon( c(twindows,rev(twindows)), c(r.mp[1,],rev(r.mp[2,])),col=rgb(0,0,0,alpha=0.25),border=NA)
	polygon( c(twindows,rev(twindows)), c(r.p[1,],rev(r.p[2,])),col=rgb(0,0.5,0,alpha=0.5),border=NA)
	polygon( c(twindows,rev(twindows)), c(r.m[1,],rev(r.m[2,])),col=rgb(0.5,0,0.5,alpha=0.5),border=NA)
}



# when plant and microbes are unevenly sized in pop, this discards unpartnered microbes. 
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
} #this function ASSUMES NO ERROR IN FITNESS, since it cannot be applied exactly as was simulated. 
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
	mup <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) )  )
	Vm <- sapply(twindows, function(t) var( colSums(simdat$Microbe[,,t]) )  )
	mum <- sapply(twindows, function(t) mean( colSums(simdat$Microbe[,,t]) )  )
	Vb <- sapply(twindows, function(t) var( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	mub <- sapply(twindows, function(t) mean( rowSums(colSums(simdat$Plant[,,,t])) + colSums(simdat$Microbe[,,t]) )  )
	return(list(Vp=Vp, Vm = Vm,Vb=Vb,mup = mup, mum=mum, mub=mub,
				cVp = abs(sqrt(Vp)/mub), cVm = abs(sqrt(Vm)/mub), cVb = abs(sqrt(Vb)/mub), 
				PVp=Vp / (Vp+Vm), PVm = Vm/(Vp+Vm)))
}#

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


#get counts of alleles segregating in the population, takes simulation output matrices subset to a particular point in time.
getsegcounts <- function(finaltime,type="plant"){
	nloci <- dim(finaltime)[[1]]
	if(type=="plant"){
		unlistall <- t(sapply(1:nloci, function(l) as.vector(finaltime[l,,])  ))
		whichsegallele <- which(sapply(1:nloci, function(l) length(unique(unlistall[l,]))) > 1)
		unlistallseg <- unlistall[whichsegallele,]
	}
	if(type=="microbe"){
		unlistall <- t(sapply(1:nloci, function(l) as.vector(finaltime[l,])  ))
		whichsegallele <- which(sapply(1:nloci, function(l) length(unique(unlistall[l,]))) > 1)
		unlistallseg <- unlistall[whichsegallele,]
	}
	return(unlistallseg)
} #type can be microbe



#plot trajectories of alleles through time; takes input subset to a particular organism -- here either host (plant) or microbe
plottraj <- function(genotimemat,type="plant",maxpos,maxneg){
	nloci <- dim(genotimemat)[1]
	nind <- dim(genotimemat)[2]
	if(type=="plant"){
		ntime <- dim(genotimemat)[4]
		fullunique <- sapply(1:nloci, function(l) sort(unique(as.vector(genotimemat[l,,,]))) )
		locusbytime <- lapply(1:nloci, function(l) sapply(2:(ntime), function(t)  sapply(1:length(fullunique[[l]]), function(a)  sum(as.vector(genotimemat[l,,,t])==fullunique[[l]][a]) ) ) )
		plot(1~c(1),pch=NA,ylim=c(0,nind*2),xlim=c(0,ntime-1),ylab="",xlab="")
	}
	if(type=="micr"){
		ntime <- dim(genotimemat)[3]
		fullunique <- sapply(1:nloci, function(l) sort(unique(as.vector(genotimemat[l,,]))) )
		locusbytime <- lapply(1:nloci, function(l) sapply(2:(ntime), function(t)  sapply(1:length(fullunique[[l]]), function(a)  sum(as.vector(genotimemat[l,,t])==fullunique[[l]][a]) ) ) )
		plot(1~c(1),pch=NA,ylim=c(0,nind),xlim=c(0,ntime-1),ylab="",xlab="")
	}
	maxall <- max(abs(c(maxpos,maxneg)))
	for(l in 1:length(locusbytime)){
		freqs <- locusbytime[[l]]
		sapply(1:nrow(freqs), function(a) lines(1:(ntime-1),freqs[a,], 
			col= ifelse(fullunique[[l]][a]>=0,
			rgb(0,0,1,alpha=abs(fullunique[[l]][a])/maxall),
			rgb(1,0,0,alpha=abs(fullunique[[l]][a])/maxall) )  ) )
	}
}



#get sojurn time for fixed alleles
sojT <- function(genotimemat,type="plant"){
	nloci <- dim(genotimemat)[1]
	nind <- dim(genotimemat)[2]
	copies <- ifelse(type=="plant",nind*2,nind)
	if(type=="plant"){
		ntime <- dim(genotimemat)[4]
		fullunique <- sapply(1:nloci, function(l) sort(unique(as.vector(genotimemat[l,,,]))) )
		locusbytime <- lapply(1:nloci, function(l) sapply(2:(ntime), function(t)  sapply(1:length(fullunique[[l]]), function(a)  sum(as.vector(genotimemat[l,,,t])==fullunique[[l]][a]) ) ) )
		finalfrq <- lapply(1:nloci, function(l) sapply(fullunique[[l]], function(a) length(which(genotimemat[l,,,ntime]==a)) )/copies )
	}
	if(type=="micr"){
		ntime <- dim(genotimemat)[3]
		fullunique <- sapply(1:nloci, function(l) sort(unique(as.vector(genotimemat[l,,]))) )
		locusbytime <- lapply(1:nloci, function(l) sapply(2:(ntime), function(t)  sapply(1:length(fullunique[[l]]), function(a)  sum(as.vector(genotimemat[l,,t])==fullunique[[l]][a]) ) ) )
		finalfrq <- lapply(1:nloci, function(l) sapply(fullunique[[l]], function(a) length(which(genotimemat[l,,ntime]==a)) )/copies )
	}
	soj <- list()
	fstate <- list()
	origingen <- list()
	for(l in 1:length(locusbytime)){
		freqs <- locusbytime[[l]]
		soj[[l]] <- c(rep(NA, times=nrow(freqs)))
		fstate[[l]] <- c(rep(NA, times=nrow(freqs)))
		origingen[[l]] <- c(rep(NA, times=nrow(freqs)))
		for(a in 1:nrow(freqs)){
			first <- which(freqs[a,]>0)[1]  
			iffix1 <- which(freqs[a,first:ncol(freqs)]==copies)[1] 
			iffix0 <- which(freqs[a,first:ncol(freqs)]==0)[1]  
			last <- ifelse(any(c(iffix1,iffix0)>0), min(c(first+iffix1,first+iffix0),na.rm=T),NA) #individually can be NA if segregating or 
			origingen[[l]][[a]]<- first # NOTE gives first generation EVEN for alleles segregating at the final timestep
			soj[[l]][a] <- last-first #returns NA when segregating at final time
			seg <- is.na(soj[[l]][a] )
			if(seg==T){
				fstate[[l]][a] <- "seg"
			} else if(!is.na(iffix1) & last == (iffix1+first)) {
				fstate[[l]][a] <- "fixed" #at least at one point
			} else{fstate[[l]][a] <- "lost"} 
		}
		#get first > 0, get first == 1 or next ==0  and then diff bt 2
	}
	return(list(soj=soj, effs = fullunique, fstate = fstate, origingen = origingen,finalfrq =finalfrq))
}


getrelfitandtrait <- function(simdat,gen,zoP,zoM,wP,wM,pfP,pfM){
	traitexp <- rowSums(colSums(simdat$Plant[,,,gen])) + colSums(simdat$Microbe[,,gen])
	micrfit   <- univar.fit( z=traitexp, zopt=zoM, sd.fit=wM)
	hostfit   <- univar.fit( z=traitexp, zopt=zoP, sd.fit=wP)
	rHfit <-     pfP*hostfit + (1-pfP)*micrfit
	rMfit <- (1-pfM)*hostfit +     pfM*micrfit
	return(data.frame(traitexp = traitexp,rfitplnt=rHfit,rfitmicr=rMfit))
}

