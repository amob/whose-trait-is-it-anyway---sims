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



#Migration: generate the rate at which, the new generation will draw individuals from another population. 
	##individuals moving this way will be incorporated using origin pop parents (t-1) chosen with that pop's fitness function, in sim.cotraitV
#Method is similar to simulations from Nuismer et al 1999, 2000, O'Brien et al 2018, but is stochastic
# migrate <- function(m,npops,sizepopV){
# 	geneflow <- matrix(NA,ncol=npops,nrow=npops)
# 	#this is a migration vector, that would be multiplied by the allele freqs across the j populations to give the frequencies in pop i (repeated each generation, below)
# 	for(i in 1:npops){
# 	individs <- sample(1:npops,size=sizepopV[i], replace=T,
# 		prob = sapply(1:npops, function(j) exp( -.5*( ( (j-i)/m)^2 ) ) ) / sum( sapply(1:npops, function(j) exp( -.5*( ( (j-i)/m)^2 ) ) )) 
# 		) #function of j which matches notation in nuismer et al 2000
# 	geneflow[,i]<- sapply(1:npops, function(pop) length(which(individs==pop)) )
# 	}
# 	return(geneflow)
# }

migrate <- function(m,npops,sizepopV,gfmat=NULL){
	geneflow <- list()
	#this is a migration vector, that would be multiplied by the allele freqs across the j populations to give the frequencies in pop i (repeated each generation, below)
	for(i in 1:npops){
		if(is.null(gfmat)) {
			prbs <- sapply(1:npops, function(j) exp( -.5*( ( (j-i)/m)^2 ) ) ) / sum( sapply(1:npops, function(j) exp( -.5*( ((j-i)/m)^2 ) ) ) )   #function of j which matches notation in nuismer et al 2000
		} else{ prbs <- gfmat[,i] }
		geneflow[[i]] <- sample(1:npops,size=sizepopV[i], replace=T, prob =prbs) 	
	}
	return(geneflow)
} #returns a list of individual assignments to origin pops for the next generation, allows unequally sized pops

###pos neg correct in mutations???
mutate.exp <- function(nL,N,lambda, prbmut) { sapply( 1:N , function(z)
			abs(rexp(nL,rate=lambda)) * ifelse(rbinom(nL,size=1,prob=0.5), 1 , -1 ) * rbinom(nL, size=1, prob = prbmut) ) #size = 1 in the last argument implies haploid?
			  }# DistofAbsValueofTraitEff * ProbofPosvNegMutation * ProbMutOccurs(u*loci) 
#high values of lambda give LOWER effect sizes of mutations on average.

horizontal <- function(nL,N,prb,genomat) {  #genomat is row are loci, cols are individuals
			 transfer <- matrix( rbinom(nL*N, size=1, prob = prb) * sample( 1:N , nL*N, replace =T) , ncol=N, nrow=nL)# with probability prb, sample a single locus from another individual, so rate of locus transfer proportional to prevalence in pop 
			 NewGenomat <- matrix( sapply(1:N, function(ind) sapply(1:nL, function(loc) 
			  				ifelse(transfer[loc,ind] > 0, genomat[loc, transfer[loc,ind] ] , genomat[loc, ind] ) 
			  				) ), ncol = N,  , byrow=F )
			  return(NewGenomat)
			  }			  
#requires and provides ind. in colums, loci in rows


sim.cotraitV <- function(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,timesteps,Lambda,mutprb ,prbHorz, pfP, pfM,ratemigr,npops,GFmat=NULL) {#,FLFC  ,startmats = "n",zoptvects = "n"){  #-- 
#SOME are VECTORS length the number of pops indicated with (V): 
	#popsize (V, note NP must equal NM currently, see below), number loci: causal then neutral, optimal phenotype (V), shallowness of fitness decline away from trait opt (V), 
			#timesteps, Lambda describes distribution of effect size of mutations, prb of mutation,  
			#fitness error (i.e. from sources random with respect to trait), prb of horizontal transfer, partner fitness feedback, r
			#rate of migration, number of pops, optional matrix of geneflow rates between pops that does not have to be symmetrical -- otherwise pops are assumed to be in straight line with IBD
#CURRENTLY NOT INCL	#startmats, if not "n", must be a list of matrices NAMED: Pa.mat, Ma.mat, Pneut.mat, Mneut.mat
#CURRENTLY NOT INCL	#zoptvects, if not "n", must be a list of vectors NAMED: PlantZ and MicrZ
#CURRENTLY NOT INCL : FLFC free living fitness cost

	#initialize
	Pa.list <- list()
	Ma.list <- list()
	Pneut.list <- list()
	Mneut.list <- list()
	for(i in 1:npops ) {
		Pa.list[[i]] <- array(0, dim= c(nlP,NP[i],2,timesteps+1) ) #store mutations on genotypes, rows have genotypes, columns have individuals, timesteps in sep matrices
			#extra dimension for diploidy, i.e. second genome in two individuals
		Ma.list[[i]] <- array(0, dim= c(nlM,NM[i],timesteps+1) ) 

		Pneut.list[[i]] <-  array(0, dim= c(nlnP,NP[i],2,timesteps+1) )
		Mneut.list[[i]] <-  array(0, dim= c(nlnM,NM[i],timesteps+1) )
	}

# 	
	for(i in 2:(timesteps+1) ) {
		
		#interact! since Pa.mat[,,i-1] and Ma.mat[,,i-1] are randomly arranged from reproduction the previous time, just pair up columns / symbiotic columns
		#joint trait value is evaluated with respect to distance from host optima and microbe optima -- then relative fitness is calculated based on the distance to both own trait optima
		 ### and how well the partner can do based on its optima, with some proportion -- i.e. assumption is that if this is a mutualism, when your partner is unfit it isn't as beneficial
		hostfit <-     univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,,i-1]), zopt = zoP, sd.fit=wP)#, fiterr=fiterrP) 
		micrfit <-     univar.fit(rowSums(colSums(Pa.mat[,,,i-1])) + colSums(Ma.mat[,,i-1]), zopt = zoM, sd.fit=wM)#, fiterr=fiterrM)#, freeliving = freeliving,FLCL=FLCL) 
		wH.e <-     pfP*hostfit + (1-pfP)*micrfit
		wM.e <- (1-pfM)*hostfit +     pfM*micrfit
	##NOTE: if wanted to do antagonistic interactions, do pfP*hostfit + (1-pfP)*((1-micrfit)/sum(1-micrfit))  #e.g. direct impact of trait on fitness, and impact passed through the inverse (e.g. perfectly negatively correlated) of impact on antagonists fitness
		#migrate
		if(is.null(GFmat)){
			migrP <- migrate(ratemigr, npops,NP)	#each population is a column and shows individuals sampled from each population (home population is the diagonal)	
			migrM <- migrate(ratemigr, npops,NM)
		} else{
			migrP <- migrate(ratemigr, npops,NP,gfmat=GFmat)	#
			migrM <- migrate(ratemigr, npops,NM,gfmat=GFmat)		
		}
		#Plant reproductions		
		#reproduce: random mating with respect to relative fitness, and unlinked loci, but finite popsize.
		#I would like each individual to choose the population of its parents as described by migrtabs, then parents within pops based on wH & wM
		seedInd <- lapply(1:npops, function(dpop) sapply(1:(NP[dpop]), function(ind) 
					sample(1:( NP[ migrP[[dpop]][ind] ] ),size=1, prob=  wH.e[[ migrP[[dpop]][ind] ]] )
					 )) #separately draw individuals based on fitness in their source population
		pollInd <- lapply(1:npops, function(dpop) sapply(1:(NP[dpop]), function(ind) 
					sample(1:( NP[ migrP[[dpop]][ind] ] ),size=1, prob=  wH.e[[ migrP[[dpop]][ind] ]] )
					 )) #separately draw individuals based on fitness

		seedPs <- lapply( 1:npops, function(pop)  #each element of seedPs should be an array of allele values [locus,individual]
					sapply(1:NP[pop], function(ind) sapply(1:nlP, function(l)  #repeat over destination pop, individuals of dest pop, and loci
						Pa.list[[ (migrP[[pop]][ind]) ]] [ l, #pick the source population for the seed parent out of the migrate lists, lth locus
														  seedInd[[pop]][ind] , #get the source parent previously drawn based on fitness for the indth dest pop individual
														sample(1:2,size=1), #randomly get one allele
														i-1 ] )) +
														 mutate.exp(nL=nlP,N=NP[pop],lambda = Lambda,prbmut=mutprb)   )
		pollPs <- lapply( 1:npops, function(pop) #repeat for prev. selected pollen parents
					sapply(1:NP[pop], function(ind) sapply(1:nlP, function(l)  
						Pa.list[[ (migrP[[pop]][ind]) ]] [ l, 
														  pollInd[[pop]][ind] , 
														sample(1:2,size=1), 
														i-1 ] )) +
														 mutate.exp(nL=nlP,N=NP[pop],lambda = Lambda,prbmut=mutprb)   )
		seedneut <- lapply( 1:npops, function(pop)  #repeat seed and pollen for neutral alleles
					sapply(1:NP[pop], function(ind) sapply(1:nlnP, function(l)  
						Pneut.list[[ (migrP[[pop]][ind]) ]] [ l, 
														  seedInd[[pop]][ind] , 
														sample(1:2,size=1), 
														i-1 ] )) +
														 mutate.exp(nL=nlnP,N=NP[pop],lambda = Lambda,prbmut=mutprb)   )
		pollneut <- lapply( 1:npops, function(pop) 
					sapply(1:NP[pop], function(ind) sapply(1:nlnP, function(l)  
						Pneut.list[[ (migrP[[pop]][ind]) ]] [ l, 
														  pollInd[[pop]][ind] , 
														sample(1:2,size=1), 
														i-1 ] )) +
														 mutate.exp(nL=nlnP,N=NP[pop],lambda = Lambda,prbmut=mutprb)   )
																
		#make a diploid offspring
		for(pop in 1:npops){    Pa.list[[pop]][,,,i] <- abind(pollPs[[pop]],    seedPs[[pop]],along=3) }
		for(pop in 1:npops){ Pneut.list[[pop]][,,,i] <- abind(pollneut[[pop]],seedneut[[pop]],along=3) }
		
		#microbial reproduction is clonal, have horizontal gene exchange after selection at probability prbHorz (per nL*N).

		micrInd <- lapply(1:npops, function(dpop) sapply(1:(NM[dpop]), function(ind) 
					sample(1:( NM[ migrM[[dpop]][ind] ] ),size=1, prob=  wM.e[[ migrM[[dpop]][ind] ]] )
					 )) #separately draw individuals based on fitness, so that can use same individuals for QTL and neutral loci
		micrgeno <- lapply(1:npops, function(dpop) sapply(1:(NM[dpop]), function(ind)		
			Ma.list[[ (migrM[[dpop]][ind]) ]] [ , #get the source pop for the indth individual in the destination pop , get all the loci
												micrInd[[pop]][ind]  , #sample a clonal parent from the source pop for the indth individual, with probability based on fitness
														i-1 ] ) +
														 mutate.exp(nL=nlM,N=NM[pop],lambda = Lambda,prbmut=mutprb)  )  
		micrgenoneut <- lapply(1:npops, function(dpop) sapply(1:(NM[dpop]), function(ind)		
			Mneut.list[[ (migrM[[dpop]][ind]) ]] [ ,
												micrInd[[pop]][ind]  , #
														i-1 ] ) +
														 mutate.exp(nL=nlnM,N=NM[pop],lambda = Lambda,prbmut=mutprb)  )  

		for(pop in 1:npops){ Ma.list[[pop]][,,i]    <- horizontal(nL = nlM,  N=NM[pop],prb=prbHorz, genomat =  micrgeno[[pop]] ) } # 	
		for(pop in 1:npops){ Mneut.list[[pop]][,,i] <- horizontal(nL = nlnM, N=NM[pop],prb=prbHorz, genomat =  micrgenoneut[[pop]] ) } # 	
	}
	result <- list(Pa.list,Ma.list,Pneut.list,Mneut.list)
	names(result) <- c("Plant","Microbe","P_neutral","M_neutral")
	return(result)
}



#######These are not adjusted for the multipop scenario, and so results must be limited to single pop outcomes 



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




# thetamat <- matrix(c(0.98 , 0.00 , 0    , 0.01 , 0.01 ,
# 					 0    , 0.98 , 0.01 , 0    , 0.01 , 
# 					 0    , 0.01 , 0.98 , 0.01 , 0    ,
# 					 0.01 , 0    , 0.01 , 0.98 , 0    ,
# 					 0.01 , 0.01 , 0    , 0    , 0.98 ), ncol=5,byrow=T)
# 
# tmpsim <-sim.cotraitV(NP=rep(50,times=5),NM=rep(50,times=5),nlP=10,nlM=20,nlnP=60,nlnM=60,
# 					zoP=c(1:5),zoM=c(2:6),wP=rep(1,times=5),wM=rep(1,times=5),timesteps=500,
# 					Lambda=10,mutprb=0.001,prbHorz=0.1,
# 					pfP=0.7,pfM=0.7,ratemigr= 0.5,npops=5,GFmat=thetamat) #note ratemigr doesn't matter if thetamat specified
# #
# # 50*5*20*2*0.001 = 10 per generation
# #unique alleles
# # sapply(1:20,function(l) table(unlist(sapply(1:5, function(pop) tmpsim$Plant[[pop]][l,,,301] ))))
# # sapply(1:20,function(l) table(unlist(sapply(1:5, function(pop) tmpsim$Microbe[[pop]][l,,301]))))
# # tmpsim$Microbe[[5]][,,301]
# # tmpsim$Plant[[5]][,,1,301]
# 
# length(unlist(sapply(1:10,function(l) table(unlist(sapply(1:5, function(pop) tmpsim$Plant[[pop]][l,,,501] )))))) #80
# length(unlist(sapply(1:20,function(l) table(unlist(sapply(1:5, function(pop) tmpsim$Microbe[[pop]][l,,501] )))))) #57
# length(unlist(sapply(1:10,function(l) table(unlist(sapply(1:5, function(pop) tmpsim$P_neutral[[pop]][l,,,501] )))))) #63
# length(unlist(sapply(1:20,function(l) table(unlist(sapply(1:5, function(pop) tmpsim$M_neutral[[pop]][l,,501] )))))) #68
# 
# 
# allelesM <- lapply(1:20,function(l) unique(as.vector(sapply(1:5, function(pop) tmpsim$Microbe[[pop]][l,,501] ))))
# allelesP <- lapply(1:10,function(l) unique(as.vector(sapply(1:5, function(pop) tmpsim$Plant[[pop]][l,,,501] ))))
# sapply(1:length(allelesP), function(l) sapply(1:5,function(p) sapply(1:length(allelesP[[l]]), function(a) length(which( unlist(tmpsim$Plant[[p]][l,,,501])==(allelesP[[l]][a]) ))/100  )) )
# sapply(1:length(allelesM), function(l) sapply(1:5,function(p) sapply(1:length(allelesM[[l]]), function(a) length(which( unlist(tmpsim$Microbe[[p]][l,,501])==(allelesM[[l]][a]) ))/50  )) )
# 
# allelesMn <- lapply(1:60,function(l) unique(as.vector(sapply(1:5, function(pop) tmpsim$M_neutral[[pop]][l,,501] ))))
# allelesPn <- lapply(1:60,function(l) unique(as.vector(sapply(1:5, function(pop) tmpsim$P_neutral[[pop]][l,,,501] ))))
# sapply(1:length(allelesPn), function(l) sapply(1:5,function(p) sapply(1:length(allelesPn[[l]]), function(a) length(which( unlist(tmpsim$P_neutral[[p]][l,,,301])==(allelesPn[[l]][a]) ))/100  )) )
# #let this serve as a reminder that if pop structure and environment are colinear there will be trouble!
# 
# lapply(1:5,function(p) mean( rowSums(colSums(tmpsim$Plant[[p]][,,,501])) + colSums(tmpsim$Microbe[[p]][,,501]) ) )
# lapply(1:5,function(p) mean( rowSums(colSums(tmpsim$Plant[[p]][,,,501]))  ) )
# lapply(1:5,function(p) mean( colSums(tmpsim$Microbe[[p]][,,501]) ) )
# 
# reordertmpsim <- lapply( 1:length(tmpsim[[1]]), function(pop) list(Plant=tmpsim$Plant[[pop]],Microbe =tmpsim$Microbe[[pop]],P_neutral=tmpsim$P_neutral[[pop]],M_neutral=tmpsim$M_neutral[[pop]]) )
# quartz()
# par(mfrow=c(2,3))
# lapply(1:length(reordertmpsim), function(pop) windowplot(2,501,10,reordertmpsim[[pop]] , ylim=c(-5,10),main=pop))

#unanswered q?
#e.g. super conserved important traits  will  have  very  little  variation.   Also  those  under  strong  (recent) directional selection

