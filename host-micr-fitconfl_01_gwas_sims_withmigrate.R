#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goals of this script are
	# to simulate an experiment from simulated pops 
	# to convert simulated loci to biallelic SNPs in a reasonable manner
	# to repeat across parameters
	# to export data in a format readable by GEMMA.
#might want to make parameter ranges shell/system variables in some way
##


source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_sim_multipopwithmigrate.R',sep="")) 
#source('~/Dropbox/host microbe trait evo and gwas/host-micr-fitconfl_01_sim_multipopwithmigrate.R') 

##EXPERIMENT FUNCTION
run.exp.allbyall<- function(popsetobj,numperpop,exp.err,nreps=1){
		#SELECT INDIVIDUALS
		numsites <- length(popsetobj$Plant)
		dimP <- lapply(1:numsites, function(pop) dim(popsetobj$Plant[[pop]]) )
		dimM <- lapply(1:numsites, function(pop) dim(popsetobj$Microbe[[pop]]) )
		finaltime <- dimP[[1]][4]
		p.geno.s <- unlist( lapply(1:numsites, function(pop)  sample(1:(dimP[[pop]][2]),numperpop, repl=F))) #which individual is sampled
		p.pop.s <- rep(1:numsites, each = numperpop) #which pop is sampled
		m.geno.s <- unlist(lapply(1:numsites, function(pop)  sample(1:(dimM[[pop]][2]),numperpop, repl=F) ) )
		m.pop.s <- rep(1:numsites, each = numperpop)

		#PHENOTYPES			
		p.breed <- sapply(1:length(p.geno.s), function(G) sum(   popsetobj$Plant[[ p.pop.s[G] ]] [,p.geno.s[G],, finaltime ] ) )
		m.breed <- sapply(1:length(m.geno.s), function(G) sum(   popsetobj$Microbe[[ m.pop.s[G] ]] [,m.geno.s[G], finaltime ] ) )
		#assume full factorial
		dat <- data.frame(traitvalue = as.vector(sapply(1:length(p.breed), function(P) sapply(1:length(m.breed), function(M) 
#					 p.breed[P] + m.breed[M] +	 rnorm(1,mean=0,sd=exp.err) ) )), #each Pbreed is a column in a matrix of rows = # mbreed, when unlist, first col, then second, then...
					 mean(p.breed[P] + m.breed[M] +	 rnorm(nreps,mean=0,sd=exp.err)) ) )), #each Pbreed is a column in a matrix of rows = # mbreed, when unlist, first col, then second, then...
						#add replicates
						  genoP = rep(p.geno.s, each = length(m.geno.s)) ,
						  popP = rep(p.pop.s,  each = length(m.geno.s)),
						  genoM = rep(m.geno.s, times = length(p.geno.s)),
						  popM = rep(m.pop.s,  times = length(p.geno.s)),
						  IDPG = rep(1:length(p.geno.s), each = length(m.geno.s)),
						  IDMG = rep(1:length(m.geno.s), times = length(p.geno.s)) )

		sel.finalT.P   <- lapply(1:numsites, function(POP)  popsetobj$Plant[[POP]] [,p.geno.s[p.pop.s==POP],,finaltime]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.M   <- lapply(1:numsites, function(POP)  popsetobj$Microbe[[POP]] [,m.geno.s[p.pop.s==POP],finaltime]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)
		sel.finalT.Pn  <- lapply(1:numsites, function(POP)  popsetobj$P_neutral[[POP]] [,p.geno.s[p.pop.s==POP],,finaltime]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.Mn  <- lapply(1:numsites, function(POP)  popsetobj$M_neutral[[POP]] [,m.geno.s[p.pop.s==POP],finaltime]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)

		causalgenos <- getalleles(sel.finalT.P,sel.finalT.M,numperpop)
		neutralgenos <- getalleles(sel.finalT.Pn,sel.finalT.Mn,numperpop)
		ind.datP <- data.frame(ID = 1:length(p.geno.s), POP = p.pop.s,	withinPOPgeno = p.geno.s)
		ind.datM <- data.frame(ID = 1:length(m.geno.s), POP = m.pop.s,	withinPOPgeno = m.geno.s)
		ind.datPM = list(plant = ind.datP, microbe= ind.datM)
		return(list(expdat = dat, 
				causalgenos = causalgenos, neutralgenos = neutralgenos,
				origcausalgenos = list(genoP=sel.finalT.P, genoM=sel.finalT.M), origneutralgenos = list(genoP=sel.finalT.Pn,genoM=sel.finalT.Mn),
				 ind.dat = ind.datPM))
}

run.exp.allbyone<- function(popsetobj,numperpop,exp.err,nreps=1) {
		#SELECT INDIVIDUALS
		numsites <- length(popsetobj$Plant)
		dimP <- lapply(1:numsites, function(pop) dim(popsetobj$Plant[[pop]]) )
		dimM <- lapply(1:numsites, function(pop) dim(popsetobj$Microbe[[pop]]) )
		finaltime <- dimP[[1]][4]
		p.geno.s <- unlist( lapply(1:numsites, function(pop)  sample(1:(dimP[[pop]][2]),numperpop, repl=F))) #which individual is sampled
		p.pop.s <- rep(1:numsites, each = numperpop) #which pop is sampled
		m.geno.s <- unlist(lapply(1:numsites, function(pop)  sample(1:(dimM[[pop]][2]),numperpop, repl=F) ) )
		m.pop.s <- rep(1:numsites, each = numperpop)
		
		whichP <- sample(1:length(p.geno.s),1) #pick the geno to be the "one" in allbyone
		whichM <- sample(1:length(m.geno.s),1) #pick the geno to be the "one" in allbyone
				
		#PHENOTYPES			
		p.breed <- sapply(1:length(p.geno.s), function(G) sum(   popsetobj$Plant[[ p.pop.s[G] ]] [,p.geno.s[G],, finaltime ] ) )
		m.breed <- sapply(1:length(m.geno.s), function(G) sum(   popsetobj$Microbe[[ m.pop.s[G] ]] [,m.geno.s[G], finaltime ] ) )
		#assume full factorial
		dat <- data.frame(traitvalue = c(sapply(1:length(p.breed), function(P) mean(p.breed[P] + m.breed[whichM] + rnorm(nreps,mean=0,sd=exp.err))  ) ,
										 sapply(1:length(m.breed), function(M) mean(p.breed[whichP] + m.breed[M] + rnorm(nreps,mean=0,sd=exp.err)) ) 
										) , #the all by one and one by all design
						  genoP = c(p.geno.s, rep(p.geno.s[whichP], times = length(m.geno.s))) ,
						  popP = c(p.pop.s, rep(p.pop.s[whichP], times = length(m.pop.s))) ,
						  genoM = c( rep(m.geno.s[whichM], times = length(p.geno.s)) , m.geno.s) ,
						  popM = c(rep(m.pop.s[whichM], times = length(p.pop.s)), m.pop.s) ,
						  IDPG = c(  1:length(p.geno.s), rep(whichP, times = length(m.geno.s))),
						  IDMG = c( rep(whichM, times = length(p.geno.s)), 1:length(m.geno.s)) )

		sel.finalT.P   <- lapply(1:numsites, function(POP)  popsetobj$Plant[[POP]] [,p.geno.s[p.pop.s==POP],,finaltime]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.M   <- lapply(1:numsites, function(POP)  popsetobj$Microbe[[POP]] [,m.geno.s[p.pop.s==POP],finaltime]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)
		sel.finalT.Pn  <- lapply(1:numsites, function(POP)  popsetobj$P_neutral[[POP]] [,p.geno.s[p.pop.s==POP],,finaltime]  ) # produces a list, each item is each pop. within these are 2 numloci X numind matrices, 1 for first allele and 1 for next
		sel.finalT.Mn  <- lapply(1:numsites, function(POP)  popsetobj$M_neutral[[POP]] [,m.geno.s[p.pop.s==POP],finaltime]  ) # produces a list, each item is each pop. within these are 1 numloci X numind matrices, (haploid!)

		causalgenos <- getalleles(sel.finalT.P,sel.finalT.M,numperpop)
		neutralgenos <- getalleles(sel.finalT.Pn,sel.finalT.Mn,numperpop)
		ind.datP <- data.frame(ID = 1:length(p.geno.s), POP = p.pop.s,	withinPOPgeno = p.geno.s)
		ind.datM <- data.frame(ID = 1:length(m.geno.s), POP = m.pop.s,	withinPOPgeno = m.geno.s)
		ind.datPM = list(plant = ind.datP, microbe= ind.datM)
		return(list(expdat = dat, expdatPexp = dat[1:(numsites*numperpop),],expdatMexp = dat[(numsites*numperpop+1):(numperpop*numsites*2),],  #this line also assumes a particular all by one and one by all design
				causalgenos = causalgenos, 	neutralgenos = neutralgenos,
				origcausalgenos = list(genoP=sel.finalT.P, genoM=sel.finalT.M), origneutralgenos = list(genoP=sel.finalT.Pn,genoM=sel.finalT.Mn),
				 ind.dat = ind.datPM))

}



##A FUNCTION TO GET GENOTYPES TO USE IN GEMMA (called in experiment function above)
###NOTES ON PROCESS
###how to get genotypes for gemma? this is a weird case with (depending on loci, time and # pops) possibly many unique allele values per locus
	##for each trait locus in each organism, figure out how many unique values there are.  that apprx is the number of total mutuations; accurately is the number of unique alleles
	## to convert to biallelic SNPs, needs to be the number of unique values , ceiling.
	## then need to assign each a SNP value (literally just the quantity assigned by the locus) and location. perfectly linked within locus and unlinked across
##just to throw markers onto genome space. with some random splitting. shouldn't affect GWAS, so this is trivial.
#what to do abou pop structure? 
	#well, we now have neutral alleles that tell us about geneflow. If we have many more neutral than causal in sims, we could just use standard methods
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
	
#	#	
	newgenomatP <- matrix(0,nrow=sum(fixedincludedP),ncol=NP*length(sel.finalT.P)*2) 	#I think it makes sense to make the default 0, each individual can have at most 2 non-zero loci in one linkage block
	colsallele1 <- seq(from=1, to = NP*length(sel.finalT.P)*2, by =2)
	colsallele2 <- seq(from=2, to = NP*length(sel.finalT.P)*2, by =2) #vectors of indices to set sequential columns as alleles for the same locus
#	#
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
							###LOCATION IS FORCED to be within 0 to 100 on a linkage group 
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


##needs to do something about the plant v microbe kin issue but is otherwise ok.

##MAKE OUTPUT FILES FOR EXPORT TO GEMMA
#NOTES
#how should the GWAS be run exaxtly?
#genotypes and phenotypes together, with only one phenotype column so genotypes must be repeated
	#presumably that means 1 file each for phenotypes with plant genotypes and phenotypes with microbe genotypes
# -- OR we could paste together genos of HOSTS and MICROBES A-LA experimental design?
	#trying this second option for now.
#GEMMA can include SNP covariate information, e.g. ; the file could include whether SNP is host or microbial; whether it is causal or not
	#this is stored in SNP names for now
#GEMMA can include individual/phenotype covariate information.
	#might use if go with different way of dealing with 2 genome issue.
#GEMMA needs plink format for genotype/phenotype information
	##MAP PED FOR EACH.
	#
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
	# 
	#MAP
		#The markers in the PED file do not need to be in genomic order: 
		#******BUT the order MAP file should align with the order of the PED file markers# ******
		# By default, each line of the MAP file describes a single marker and must contain exactly 4 columns:
		#


makegwasfiles <- function(expset,expdat,type="HOLO"){
	if(type=="HOLO"){
	 RGP <- t( sapply(seq(from=1, to =ncol(expset$causalgenos$genoP),by=2), function(Ind) 
					c(t(  rbind(expset$causalgenos$genoP,expset$neutralgenos$genoP)[,c(Ind,Ind+1)] )) ) )  #concatenating plant genotype matrix, merging individual columns so that loci rather than individuals are repeated, repeating across loci and transposing to rowxcol = individuals x loci (2 cols per locus)
	 geno <- cbind( paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),  #Family ID when missing, plink wants this= to individual ID
			paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""), # Individual ID
			rep(0,times=nrow(expdat)) , rep(0,times=nrow(expdat)) ,rep(0,times=nrow(expdat)) , # 
			expdat$traitvalue,
		    (RGP+1) [expdat$IDPG,] , # adding 1 so bt 1 and 2, and then grabbing rows as they are used in the experiment.
 			#then
	 	   t(rbind( expset$causalgenos$genoM[rep(1:nrow(expset$causalgenos$genoM)  ,each=2), ], 
	 	   			## duplicating rows for pretend diploid microbe genotype matrix (lociin2rowsxindividuals), 
	 	   		   expset$neutralgenos$genoM[rep(1:nrow(expset$neutralgenos$genoM) ,each=2), ])+1) [expdat$IDMG,]  
	 	   		   #concatenating neutral loci treated the same way, transposing to ind x loci(2 cols per), adding 1 so bt 1 and 2, 
	 	   		   #and then grabbing individuals in rows as they are used in the experiment. 
	 	   )
	 map <- rbind( cbind( paste("scaffold_",c(expset$causalgenos$locusdatP$linkage,expset$neutralgenos$locusdatP$linkage),sep=""),  #      chromosome (1-22, X, Y or 0 if unplaced)
	#linkage group must be alphanumeric not numeric if > 22 chromosomes. so paste "scaffold" first "--allow-extra-chr" flag to make plink accept them.
			c(paste("cP",1:nrow(expset$causalgenos$locusdatP),sep=""),paste("nP",1:nrow(expset$neutralgenos$locusdatP),sep="") ),#      rs# or snp identifier
	 		rep(0,times=nrow(expset$causalgenos$locusdatP)+ nrow(expset$neutralgenos$locusdatP)), #      Genetic distance (morgans) # For basic association testing, the genetic distance column can be set at 0. 
	 		c(expset$causalgenos$locusdatP$location,expset$neutralgenos$locusdatP$location) ), #      Base-pair position (bp units)
	 		#
	   		cbind( paste("scaffold_",c(expset$causalgenos$locusdatM$linkage,expset$neutralgenos$locusdatM$linkage),sep=""),  #same for microbe
		 	c(paste("cM",1:nrow(expset$causalgenos$locusdatM),sep=""),paste("nM",1:nrow(expset$neutralgenos$locusdatM),sep="") ),
		 	rep(0,times=nrow(expset$causalgenos$locusdatM)+ nrow(expset$neutralgenos$locusdatM)), 
		 	c(expset$causalgenos$locusdatM$location,expset$neutralgenos$locusdatM$location) ) 
			)
	 return(list(map=map,geno=geno))
	} else if(type=="PLANT"){
	 RGP <- t( sapply(seq(from=1, to =ncol(expset$causalgenos$genoP),by=2), function(Ind) 
					c(t(  rbind(expset$causalgenos$genoP,expset$neutralgenos$genoP)[,c(Ind,Ind+1)] )) ) )  #concatenating plant genotype matrix, merging individual columns so that loci rather than individuals are repeated, repeating across loci and transposing to rowxcol = individuals x loci (2 cols per locus)
	 geno <- cbind( paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""), #
			paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),
			rep(0,times=nrow(expdat)) , rep(0,times=nrow(expdat)) ,rep(0,times=nrow(expdat)) ,
			expdat$traitvalue,
		    (RGP+1) [expdat$IDPG,]  # adding 1 so bt 1 and 2, and then grabbing rows as they are used in the experiment.
	  	   )
	 map <-  cbind( paste("scaffold_",c(expset$causalgenos$locusdatP$linkage,expset$neutralgenos$locusdatP$linkage),sep=""),  #  now separately for plants
			c(paste("cP",1:nrow(expset$causalgenos$locusdatP),sep=""),paste("nP",1:nrow(expset$neutralgenos$locusdatP),sep="") ),
	 		rep(0,times=nrow(expset$causalgenos$locusdatP)+ nrow(expset$neutralgenos$locusdatP)),
	 		c(expset$causalgenos$locusdatP$location,expset$neutralgenos$locusdatP$location) ) 
	 return(list(map=map,geno=geno))
	} else if(type=="MICR"){
	 geno <- cbind( paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),  #microbes again, but alone
			paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),
			rep(0,times=nrow(expdat)) , rep(0,times=nrow(expdat)) ,rep(0,times=nrow(expdat)) ,
			expdat$traitvalue,
		 	t(rbind( expset$causalgenos$genoM[rep(1:nrow(expset$causalgenos$genoM)  ,each=2), ], 
	 	   		   expset$neutralgenos$genoM[rep(1:nrow(expset$neutralgenos$genoM) ,each=2), ])+1) [expdat$IDMG,]  
	 	   )
	 map <- cbind( paste("scaffold_",c(expset$causalgenos$locusdatM$linkage,expset$neutralgenos$locusdatM$linkage),sep=""), #and separately for microbes
		 	c(paste("cM",1:nrow(expset$causalgenos$locusdatM),sep=""),paste("nM",1:nrow(expset$neutralgenos$locusdatM),sep="") ),
		 	rep(0,times=nrow(expset$causalgenos$locusdatM)+ nrow(expset$neutralgenos$locusdatM)), 
		 	c(expset$causalgenos$locusdatM$location,expset$neutralgenos$locusdatM$location) ) 
	 return(list(map=map,geno=geno))
	} else if(type=="ALL"){
	 RGP <- t( sapply(seq(from=1, to =ncol(expset$causalgenos$genoP),by=2), function(Ind) 
					c(t(  rbind(expset$causalgenos$genoP,expset$neutralgenos$genoP)[,c(Ind,Ind+1)] )) ) )  #concatenating plant genotype matrix, merging individual columns so that loci rather than individuals are repeated, repeating across loci and transposing to rowxcol = individuals x loci (2 cols per locus)
	 HOLOgeno <- cbind( paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),  #Family ID when missing, plink wants this= to individual ID
			paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""), # Individual ID
			rep(0,times=nrow(expdat)) , rep(0,times=nrow(expdat)) ,rep(0,times=nrow(expdat)) , # 
			expdat$traitvalue,
		    (RGP+1) [expdat$IDPG,] , # adding 1 so bt 1 and 2, and then grabbing rows as they are used in the experiment.
 			#then
	 	   t(rbind( expset$causalgenos$genoM[rep(1:nrow(expset$causalgenos$genoM)  ,each=2), ], 
	 	   			## duplicating rows for pretend diploid microbe genotype matrix (lociin2rowsxindividuals), 
	 	   		   expset$neutralgenos$genoM[rep(1:nrow(expset$neutralgenos$genoM) ,each=2), ])+1) [expdat$IDMG,]  
	 	   		   #concatenating neutral loci treated the same way, transposing to ind x loci(2 cols per), adding 1 so bt 1 and 2, 
	 	   		   #and then grabbing individuals in rows as they are used in the experiment.
	 	   )
	 HOLOmap <- rbind( cbind( paste("scaffold_",c(expset$causalgenos$locusdatP$linkage,expset$neutralgenos$locusdatP$linkage),sep=""),  #      chromosome (1-22, X, Y or 0 if unplaced)
	#linkage group must be alphanumeric not numeric if > 22 chromosomes. so paste "scaffold" first "--allow-extra-chr" flag to make plink accept them.
			c(paste("cP",1:nrow(expset$causalgenos$locusdatP),sep=""),paste("nP",1:nrow(expset$neutralgenos$locusdatP),sep="") ),#      rs# or snp identifier
	 		rep(0,times=nrow(expset$causalgenos$locusdatP)+ nrow(expset$neutralgenos$locusdatP)), #      Genetic distance (morgans) # For basic association testing, the genetic distance column can be set at 0. 
	 		c(expset$causalgenos$locusdatP$location,expset$neutralgenos$locusdatP$location) ), #      Base-pair position (bp units)
	 		#
	   		cbind( paste("scaffold_",c(expset$causalgenos$locusdatM$linkage,expset$neutralgenos$locusdatM$linkage),sep=""),  #same for microbe
		 	c(paste("cM",1:nrow(expset$causalgenos$locusdatM),sep=""),paste("nM",1:nrow(expset$neutralgenos$locusdatM),sep="") ),
		 	rep(0,times=nrow(expset$causalgenos$locusdatM)+ nrow(expset$neutralgenos$locusdatM)), 
		 	c(expset$causalgenos$locusdatM$location,expset$neutralgenos$locusdatM$location) ) 
			)
	 PLANTgeno <- cbind( paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""), #plants again but alone
			paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),
			rep(0,times=nrow(expdat)) , rep(0,times=nrow(expdat)) ,rep(0,times=nrow(expdat)) ,
			expdat$traitvalue,
		    (RGP+1) [expdat$IDPG,]  # adding 1 so bt 1 and 2, and then grabbing rows as they are used in the experiment.
	  	   )
	 PLANTmap <-  cbind( paste("scaffold_",c(expset$causalgenos$locusdatP$linkage,expset$neutralgenos$locusdatP$linkage),sep=""),  #  now separately for plants
			c(paste("cP",1:nrow(expset$causalgenos$locusdatP),sep=""),paste("nP",1:nrow(expset$neutralgenos$locusdatP),sep="") ),
	 		rep(0,times=nrow(expset$causalgenos$locusdatP)+ nrow(expset$neutralgenos$locusdatP)),
	 		c(expset$causalgenos$locusdatP$location,expset$neutralgenos$locusdatP$location) ) 
	 MICRgeno <- cbind( paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),  #microbes again, but alone
			paste("p",expdat$IDPG,"m",expdat$IDMG,sep=""),
			rep(0,times=nrow(expdat)) , rep(0,times=nrow(expdat)) ,rep(0,times=nrow(expdat)) ,
			expdat$traitvalue,
		 	t(rbind( expset$causalgenos$genoM[rep(1:nrow(expset$causalgenos$genoM)  ,each=2), ], 
	 	   		   expset$neutralgenos$genoM[rep(1:nrow(expset$neutralgenos$genoM) ,each=2), ])+1) [expdat$IDMG,]  
	 	   )
	 MICRmap <- cbind( paste("scaffold_",c(expset$causalgenos$locusdatM$linkage,expset$neutralgenos$locusdatM$linkage),sep=""), #and separately for microbes
		 	c(paste("cM",1:nrow(expset$causalgenos$locusdatM),sep=""),paste("nM",1:nrow(expset$neutralgenos$locusdatM),sep="") ),
		 	rep(0,times=nrow(expset$causalgenos$locusdatM)+ nrow(expset$neutralgenos$locusdatM)), 
		 	c(expset$causalgenos$locusdatM$location,expset$neutralgenos$locusdatM$location) ) 
	 return(list(HOLOmap =HOLOmap,HOLOgeno=HOLOgeno,PLANTmap =PLANTmap,PLANTgeno =PLANTgeno,MICRmap= MICRmap,MICRgeno=MICRgeno))
	} else {print("unknown type")}
 }
 
 
  #
##PASS .ped and .map to PLINK to get binary formats .bim .fam .bed
##PASS to GEMMA to run models.





#make a matrix for geneflow, function. 
#this one currently breaks, when the conditions are less easy to fill -- more pops, more pops per pop involved in genflow
filltheta <- function(npopsource, geneflowrate,npops){
	popinbredrate <- 1 - geneflowrate*npopsource
	thetamat <- diag(popinbredrate,nrow=npops)
	for(i in 1:(ncol(thetamat)-1)){
		if(sum(sign(thetamat[,i])) < (npopsource+1)  ){
				whichcanbesampled <-  (1:nrow(thetamat)) [  (rowSums(sign(thetamat)) < (npopsource+1) ) & thetamat[(1:nrow(thetamat)),i] == 0 ]  
				if(length(whichcanbesampled)>1){
					sourcepops <- sample(whichcanbesampled ,1 + npopsource - sum(sign(thetamat[,i])) ,repl=F)
					} else{ sourcepops <- whichcanbesampled }
				thetamat[sourcepops,i] <- geneflowrate #populate the column
				thetamat[i,(i+1):ncol(thetamat)]  <- thetamat[(i+1):ncol(thetamat),i]
			} else{}
	}
	return(thetamat)
}



 ####GENERATE OUTPUT FILES FOR SIMULATIONS WITH GWAS
###simulate multiple populations; run experiment on result; create genotypes. 
#list parameters; keep everything the same except optima?

npopsource <- 4
geneflowrate <- 0.02
npops <- 200
randtheta <- 	thetamat <- diag(1-npopsource*geneflowrate,nrow=npops)
while(any(colSums(sign(randtheta)) < (npopsource+1) ) ) {
	randtheta <- filltheta(npopsource =npopsource, geneflowrate = geneflowrate, npops = npops)
}


#temporary sim for testing code
# thetamat <- matrix(c(0.98 , 0.00 , 0    , 0.01 , 0.01 ,
# 					 0    , 0.98 , 0.01 , 0    , 0.01 , 
# 					 0    , 0.01 , 0.98 , 0.01 , 0    ,
# 					 0.01 , 0    , 0.01 , 0.98 , 0    ,
# 					 0.01 , 0.01 , 0    , 0    , 0.98 ), ncol=5,byrow=T)
# 	npops <- 5
# popset <-sim.cotraitV(NP=rep(50,times=npops),NM=rep(50,times=npops),nlP=10,nlM=20,nlnP=100,nlnM=200,
# 					zoP=seq(from=1,to=5,length.out=npops),zoM=seq(from=2,to=6,length.out=npops),wP=rep(1,times=npops),wM=rep(1,times=npops),timesteps=5,
# 					Lambda=10,mutprb=0.1,prbHorz=0.1,
# 					pfP=0.7,pfM=0.7,ratemigr= 0.5,npops=npops,GFmat=thetamat) #note ratemigr doesn't matter if thetamat specified

#sim for running code and testing hypotheses
# popset <-sim.cotraitV(NP=rep(50,times=npops),NM=rep(50,times=npops),nlP=10,nlM=20,nlnP=100,nlnM=200,
# 					zoP=seq(from=1,to=5,length.out=npops),zoM=seq(from=2,to=6,length.out=npops),wP=rep(1,times=npops),wM=rep(1,times=npops),timesteps=500,
# 					Lambda=10,mutprb=0.001,prbHorz=0.1,
# 					pfP=0.7,pfM=0.7,ratemigr= 0.5,npops=npops,GFmat=randtheta) #note ratemigr doesn't matter if thetamat specified
#temporary commnt out
#expecting about 70 causal alleles ea per plnt and microbe, and maybe 150 each neutral ones.
#save(popset,file=paste(Sys.getenv("SCRATCH"),'/popset.RData',sep=""))
load(file=paste(Sys.getenv("SCRATCH"),'/popset.RData',sep=""))

reduceset <- sort(sample(1:npops, 40, repl=F)) #for below, needs to be a subsample of 60 pops
red_popset <- list( Plant = lapply(reduceset, function(pop) popset$Plant[[pop]]),       Microbe = lapply(reduceset, function(pop) popset$Microbe[[pop]]),
			    P_neutral = lapply(reduceset, function(pop) popset$P_neutral[[pop]]), M_neutral = lapply(reduceset, function(pop) popset$M_neutral[[pop]])) 

save(red_popset,file=paste(Sys.getenv("SCRATCH"),'/red_popset.RData',sep=""))


expsetabo <- run.exp.allbyone(popset, numperpop= 4,exp.err=0.05,nreps=4) #2 plant and 2  micr from 200 pops , 4*(800*1 + 1*800) , = 6400 exp; 3200 ea gwas
save(expsetabo,file=paste(Sys.getenv("SCRATCH"),'/expset_abo.RData',sep=""))
plinkabo <- makegwasfiles(expsetabo,expsetabo$expdat,type="HOLO") #<- function(expset,name_append){
	write.table(plinkabo$geno, file=paste(Sys.getenv("SCRATCH"), "/HOLOevosims_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
  	write.table(plinkabo$map,  file=paste(Sys.getenv("SCRATCH"), "/HOLOevosims_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboP <- makegwasfiles(expsetabo,expdat=expsetabo$expdatPexp,type="PLANT") #<- function(expset,name_append){
	write.table(plinkaboP$geno,file=paste(Sys.getenv("SCRATCH"),"/PLANTevosims_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboP$map, file=paste(Sys.getenv("SCRATCH"),"/PLANTevosims_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboM <- makegwasfiles(expsetabo,expdat=expsetabo$expdatMexp,type="MICR") #<- function(expset,name_append){
	write.table(plinkaboM$geno, file=paste(Sys.getenv("SCRATCH"), "/MICRevosims_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboM$map,  file=paste(Sys.getenv("SCRATCH"), "/MICRevosims_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)




expsetaba <- run.exp.allbyall(red_popset,numperpop= 1,exp.err=0.05,nreps=4) #1 from 40 pops, 40x40x4 = 6400 exp
save(expsetaba,file=paste(Sys.getenv("SCRATCH"),'/expset_aba.RData',sep=""))
plinkaba <- makegwasfiles(expsetaba,expsetaba$expdat,type="ALL") #<- function(expset,name_append){
 	write.table(plinkaba$HOLOgeno, file=paste(Sys.getenv("SCRATCH"), "/HOLOevosims_ABA.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaba$PLANTgeno,file=paste(Sys.getenv("SCRATCH"),"/PLANTevosims_ABA.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaba$MICRgeno, file=paste(Sys.getenv("SCRATCH"), "/MICRevosims_ABA.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
  	write.table(plinkaba$HOLOmap,  file=paste(Sys.getenv("SCRATCH"), "/HOLOevosims_ABA.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaba$PLANTmap, file=paste(Sys.getenv("SCRATCH"),"/PLANTevosims_ABA.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaba$MICRmap,  file=paste(Sys.getenv("SCRATCH"), "/MICRevosims_ABA.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

 
