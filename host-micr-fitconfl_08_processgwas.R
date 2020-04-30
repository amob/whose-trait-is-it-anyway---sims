

# holoout <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/tmp_HOLOgemma.assoc.txt",header=T,sep="\t")
# holooutK <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/tmp_HOLOgemmaK.assoc.txt",header=T,sep="\t")

# PARAMS WeRE
#(NP=rep(50,times=npops),NM=rep(50,times=npops),nlP=10,nlM=20,nlnP=100,nlnM=200,
# 					zoP=seq(from=1,to=5,length.out=npops),zoM=seq(from=2,to=6,length.out=npops),wP=rep(1,times=npops),wM=rep(1,times=npops),timesteps=500,
# 					Lambda=10,mutprb=0.001,prbHorz=0.1,
# 					pfP=0.7,pfM=0.7,ratemigr= 0.5,npops=npops,GFmat=randtheta) #note ratemigr doesn't matter if thetamat specified
holoout <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/HOLOgemmaABO.assoc.txt",header=T,sep="\t")
holooutK <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/HOLOgemmaKABO.assoc.txt",header=T,sep="\t")




holoout$known_pos <- sapply(holoout$rs, function(z) length(grep("c",z))==1)
holoout$known_neg <- !holoout$known_pos

holooutK$known_pos <- sapply(holooutK$rs, function(z) length(grep("c",z))==1)
holooutK$known_neg <- !holooutK$known_pos

p4noK <- findInterval(holoout$p_score,sort(holoout$p_score))/nrow(holoout)
p4K <- findInterval(holooutK$p_score,sort(holooutK$p_score))/nrow(holooutK)


holooutK$false_pos<- (holooutK$p_score < 0.05 & holooutK$known_pos == F)
holooutK$true_pos<- (holooutK$p_score < 0.05 & holooutK$known_pos == T)
holooutK$true_neg <- (holooutK$p_score > 0.05 & holooutK$known_pos == F)
holooutK$false_neg <- (holooutK$p_score > 0.05 & holooutK$known_pos == T)


holoout$false_pos<- (p4noK < 0.05 & holoout$known_pos == F)
holoout$true_pos<- (p4noK < 0.05 & holoout$known_pos == T)
holoout$true_neg <- (p4noK > 0.05 & holoout$known_pos == F)
holoout$false_neg <- (p4noK > 0.05 & holoout$known_pos == T)

#plant plantpos micr micrpos
known_PposK <- sapply(holooutK$rs, function(z) length(grep("cP",z))==1)
known_PnK <- sapply(holooutK$rs, function(z) length(grep("nP",z))==1)
known_MposK <- sapply(holooutK$rs, function(z) length(grep("cM",z))==1)
known_MnK <- sapply(holooutK$rs, function(z) length(grep("nM",z))==1)
catsnpK <- rep(NA,times=nrow(holooutK))
		catsnpK[known_PposK] <- "Pcausal"
		catsnpK[known_MposK] <- "Mcausal"
		catsnpK[known_PnK] <- "Pneut"
		catsnpK[known_MnK] <- "Mneut"
known_Ppos <- sapply(holoout$rs, function(z) length(grep("cP",z))==1)
known_Pn <- sapply(holoout$rs, function(z) length(grep("nP",z))==1)
known_Mpos <- sapply(holoout$rs, function(z) length(grep("cM",z))==1)
known_Mn <- sapply(holoout$rs, function(z) length(grep("nM",z))==1)
catsnp <- rep(NA,times=nrow(holoout))
		catsnp[known_Ppos] <- "Pcausal"
		catsnp[known_Mpos] <- "Mcausal"
		catsnp[known_Pn] <- "Pneut"
		catsnp[known_Mn] <- "Mneut"



startcol <- grep("false_pos",colnames(holoout))
startcolK <- grep("false_pos",colnames(holooutK))

conttabprp <- matrix(
	colSums(holoout[,startcol:(startcol+3)]) /
		rep(c(sum(holoout$known_neg),sum(holoout$known_pos)),times=2),  
		#c(f+ t+ t- f- ) /  c(allneg, allpos, allneg,allpos)
	byrow=T,ncol=2
		)
conttab <- 	colSums(holoout[,startcol:(startcol+3)])

#
conttabKprp <- matrix(
	colSums(holooutK[,startcolK:(startcolK+3)]) /
		rep(c(sum(holooutK$known_neg),sum(holooutK$known_pos)),times=2),
	byrow=T,ncol=2
		) # false pos, true pos, true neg, false neg /  c(knownneg, knownpos, knownneg, knownpos)
		
conttabK <- 	colSums(holooutK[,startcolK:(startcolK+3)])


cattots <- table(paste(catsnp,p4noK < 0.05,sep="")) [c(4,2,3,1,8,6,7,5)] 
# to match order above
#matching denominator to above matrix
catshouldbesums <- c(  rep( c(sum(cattots[c(1,3)]),sum(cattots[c(2,4)] )),times=2 )  ,  
	rep( c(sum(cattots[c(5,7)]),sum(cattots[c(6,8)])),times=2 ) ) 
cattotsK <- table(paste(catsnpK,holooutK$p_score < 0.05,sep="")) [c(4,2,3,1,8,6,7,5)] #for KINSHIP
catshouldbesumsK <- c(  rep( c(sum(cattotsK[c(1,3)]),sum(cattotsK[c(2,4)] )),times=2 )  ,  
	rep( c(sum(cattotsK[c(5,7)]),sum(cattotsK[c(6,8)])),times=2 ) ) 
			



pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_.pdf",height=8,width=8)
par(mfrow=c(2,3))

image(conttabprp,xaxt="n",yaxt="n",main="Together")
	mtext("w/o Kinship, interval p",side=2,line=2)
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		conttab,rep("/",times=4),
		rep(c(sum(holoout$known_neg),sum(holoout$known_pos)),times=2),
		rep("=",times=4),
		round(as.vector(t(conttabprp)),digits=2), sep="" ))
image(t(matrix( (cattots/catshouldbesums)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Microbe loci")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		cattots[1:4],rep("/",times=4), catshouldbesums[1:4], rep("=",times=4),
		round(  (cattots/catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (cattots/catshouldbesums)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Plant loci")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		cattots[5:8],rep("/",times=4), catshouldbesums[5:8], rep("=",times=4),
		round(  (cattots/catshouldbesums)[5:8],digits=2), sep="" ))


image(conttabKprp,xaxt="n",yaxt="n",main="")
	mtext("w/ Kinship, model p",side=2,line=2)
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		conttabK,rep("/",times=4),
		rep(c(sum(holooutK$known_neg),sum(holooutK$known_pos)),times=2),
		rep("=",times=4),
		round(as.vector(t(conttabKprp)),digits=2), sep="" ))
image(t(matrix( (cattotsK/catshouldbesumsK)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		cattotsK[1:4],rep("/",times=4), catshouldbesumsK[1:4], rep("=",times=4),
		round(  (cattotsK/catshouldbesumsK)[1:4],digits=2), sep="" ))
image(t(matrix( (cattotsK/catshouldbesumsK)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		cattotsK[5:8],rep("/",times=4), catshouldbesumsK[5:8], rep("=",times=4),
		round(  (cattotsK/catshouldbesumsK)[5:8],digits=2), sep="" ))

dev.off()


