
addcols <- function(outfile, k=NULL){
	outfile$known_pos <- sapply(outfile$rs, function(z) length(grep("c",z))==1)
	outfile$known_neg <- !outfile$known_pos
	if(is.null(k)){
p <- findInterval(outfile$p_score,sort(outfile$p_score))/nrow(outfile)
#  	p <- findInterval(outfile$beta,sort(outfile$beta))/nrow(outfile) #fewer total loci identified in both, regardless of how count. BUT makes the holo analysis more effective.

	} else if(k==TRUE){
		p <- outfile$p_score
	}
	outfile$p_used <- p 
	outfile$false_pos <- (p < 0.05 & outfile$known_pos == F)
	outfile$true_pos <-  (p < 0.05 & outfile$known_pos == T)
	outfile$true_neg <-  (p > 0.05 & outfile$known_pos == F)
	outfile$false_neg <- (p > 0.05 & outfile$known_pos == T)

	known_Ppos <- sapply(outfile$rs, function(z) length(grep("cP",z))==1)
	known_Pn <- sapply(outfile$rs, function(z) length(grep("nP",z))==1)
	known_Mpos <- sapply(outfile$rs, function(z) length(grep("cM",z))==1)
	known_Mn <- sapply(outfile$rs, function(z) length(grep("nM",z))==1)
	catsnp <- rep(NA,times=nrow(outfile))
		if(any(known_Ppos)==TRUE){catsnp[known_Ppos] <- "Pcausal"}
		if(any(known_Mpos)==TRUE){catsnp[known_Mpos] <- "Mcausal"}
		if(any(known_Pn)==TRUE){catsnp[known_Pn] <- "Pneut"}
		if(any(known_Mn)==TRUE){catsnp[known_Mn] <- "Mneut"}
	outfile$catsnp <- catsnp
	return(outfile)
}

summarizeSNPcalls <- function(outfilewcols){
	startcol <- grep("false_pos",colnames(outfilewcols))
	conttabprp <- matrix(
		colSums(outfilewcols[,startcol:(startcol+3)]) /
		rep(c(sum(outfilewcols$known_neg),sum(outfilewcols$known_pos)),times=2),  
		#c(f+ t+ t- f- ) /  c(allneg, allpos, allneg,allpos)
		byrow=T,ncol=2
		)
	conttab <- 	colSums(outfilewcols[,startcol:(startcol+3)])
	possout <- c("MneutTRUE",  "McausalTRUE",   "MneutFALSE",    "McausalFALSE", "PneutTRUE",  "PcausalTRUE",   "PneutFALSE",    "PcausalFALSE")
	cattots <- sapply(1:length(possout), function(o)  sum(paste(outfilewcols$catsnp,outfilewcols$p_used < 0.05,sep="")==possout[o] ) )
	names(cattots)<-possout
	catshouldbesums <- c(  rep( c(sum(cattots[c(1,3)]),sum(cattots[c(2,4)] )),times=2 )  ,  
	rep( c(sum(cattots[c(5,7)]),sum(cattots[c(6,8)])),times=2 ) ) 
	catshouldbesums[is.na(catshouldbesums)] <- 0
	return(list(conttabprp=conttabprp,conttab=conttab, cattots=cattots,catshouldbesums=catshouldbesums))
}

connectIG <- function(inputlociinfo, inputgwasinfo){
	matchloci <- sapply(1:nrow(inputlociinfo), function(z) which(as.character(inputgwasinfo$rs) == as.character(inputlociinfo$locusname[z])))
	is.inHABA <- sapply(1:length(matchloci), function(z) ifelse(length(matchloci[[z]])==1,matchloci[[z]],NA))
	connectedinfo <- inputgwasinfo[is.inHABA,]
	return(connectedinfo)
}


# holoout <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/tmp_HOLOgemma.assoc.txt",header=T,sep="\t")
# holooutK <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/tmp_HOLOgemmaK.assoc.txt",header=T,sep="\t")

# PARAMS WeRE
#(NP=rep(50,times=npops),NM=rep(50,times=npops),nlP=10,nlM=20,nlnP=100,nlnM=200,
# 					zoP=seq(from=1,to=5,length.out=npops),zoM=seq(from=2,to=6,length.out=npops),wP=rep(1,times=npops),wM=rep(1,times=npops),timesteps=500,
# 					Lambda=10,mutprb=0.001,prbHorz=0.1,
# 					pfP=0.7,pfM=0.7,ratemigr= 0.5,npops=npops,GFmat=randtheta) #note ratemigr doesn't matter if thetamat specified
holooutabo <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/HOLOgemmaABO.assoc.txt",header=T,sep="\t")
holooutKabo <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/HOLOgemmaKABO.assoc.txt",header=T,sep="\t")
holooutaba <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/HOLOgemmaABA.assoc.txt",header=T,sep="\t")
holooutKaba <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/HOLOgemmaKABA.assoc.txt",header=T,sep="\t")
microutabo <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO.assoc.txt",header=T,sep="\t")
microutaba <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABA.assoc.txt",header=T,sep="\t")
plantoutabo <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO.assoc.txt",header=T,sep="\t")
plantoutaba <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABA.assoc.txt",header=T,sep="\t")


holooutabo <- addcols(holooutabo)
holooutKabo <- addcols(holooutKabo,k=TRUE)
holooutaba <- addcols(holooutaba)
holooutKaba <- addcols(holooutKaba,k=TRUE)
plantoutabo <- addcols(plantoutabo)
microutabo <- addcols(microutabo)
plantoutaba <- addcols(plantoutaba)
microutaba <- addcols(microutaba)

summarizehabo  <- summarizeSNPcalls(holooutabo)
summarizehKabo <- summarizeSNPcalls(holooutKabo)
summarizehaba  <- summarizeSNPcalls(holooutaba)
summarizehKaba <- summarizeSNPcalls(holooutKaba)
summarizepabo  <- summarizeSNPcalls(plantoutabo)
summarizemabo  <- summarizeSNPcalls(microutabo)
summarizepaba  <- summarizeSNPcalls(plantoutaba)
summarizemaba  <- summarizeSNPcalls(microutaba)


plantLinfoaba <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABA.csv",header=T)
micrLinfoaba <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABA.csv",header=T)
plantLinfoabo <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO.csv",header=T)
micrLinfoabo <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO.csv",header=T)
holoLinfoaba <- rbind(plantLinfoaba,micrLinfoaba)
holoLinfoabo <- rbind(plantLinfoabo,micrLinfoabo)


# CIGph_aba <- connectIG(plantLinfoaba,holooutaba)
CIGp_abo <- connectIG(plantLinfoabo,plantoutabo)
CIGm_abo <- connectIG(micrLinfoabo,microutabo)
CIGh_abo <- connectIG(holoLinfoabo,holooutabo)
CIGp_aba <- connectIG(plantLinfoaba,plantoutaba)
CIGm_aba <- connectIG(micrLinfoaba,microutaba)
CIGh_aba <- connectIG(holoLinfoaba,holooutaba)

CIGh_aba$is.plant <- sapply(1:nrow(CIGh_aba), function(z) length(grep("P",as.character(CIGh_aba$rs[z])) )==1) 
CIGh_abo$is.plant <- sapply(1:nrow(CIGh_abo), function(z) length(grep("P",as.character(CIGh_abo$rs[z])) )==1) 



#visualize what is lost, it is the small effect size loci
# hist((CIGp_abo$beta),freq=F)
# hist(plantLinfoabo$reststate,freq=T,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,200))#ylim=c(0,2.5))
# hist(plantLinfoabo$reststate[!is.na(CIGp_abo$beta)],freq=T,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/dens_loci_tossed.pdf",height=4,width=6)
par(mfrow=c(2,3))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
hist(holoLinfoabo$reststate,freq=T,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,400),main="") #ylim=c(0,200))#
	hist(holoLinfoabo$reststate[!is.na(CIGh_abo$beta)],freq=T,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Holo-estimated",side=3, line=0.5)
	mtext("Frequency",side=2, line=3)
	mtext("ABO",side=2, line=2,adj = -0.5)
hist(plantLinfoabo$reststate,freq=T,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,400),main="") #ylim=c(0,200))#
	hist(plantLinfoabo$reststate[!is.na(CIGp_abo$beta)],freq=T,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Plant-estimated",side=3, line=0.5)
hist(micrLinfoabo$reststate,freq=T,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,400),main="") #ylim=c(0,200))#
	hist(micrLinfoabo$reststate[!is.na(CIGm_abo$beta)],freq=T,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Microbe-estimated",side=3, line=0.5)
hist(holoLinfoabo$reststate,freq=F,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,3.5),main="") #ylim=c(0,200))#
	hist(holoLinfoabo$reststate[!is.na(CIGh_abo$beta)],freq=F,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Density",side=2, line=3)
hist(plantLinfoabo$reststate,freq=F,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,3.5),main="") #ylim=c(0,200))#
	hist(plantLinfoabo$reststate[!is.na(CIGp_abo$beta)],freq=F,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Known effect",side=1, line=2)
hist(micrLinfoabo$reststate,freq=F,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,3.5),main="") #ylim=c(0,200))#
	hist(micrLinfoabo$reststate[!is.na(CIGm_abo$beta)],freq=F,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
#both distributions are biased positive, and  more are thrown out in the smaller effect size range.
# hist(holoLinfoaba$reststate,freq=F,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,3.5),main="") #ylim=c(0,200))#
# 	hist(holoLinfoaba$reststate[!is.na(CIGh_aba$beta)],freq=F,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
# 	mtext("ABA",side=2, line=3)
# hist(plantLinfoaba$reststate,freq=F,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,3.5),main="") #ylim=c(0,200))#
# 	hist(plantLinfoaba$reststate[!is.na(CIGp_aba$beta)],freq=F,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
# 	mtext("Known effect",side=1, line=2)
# hist(micrLinfoaba$sreststate,freq=F,breaks=seq(from=-2,to=2,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,3.5),main="") #ylim=c(0,200))#
# 	hist(micrLinfoaba$reststate[!is.na(CIGm_aba$beta)],freq=F,breaks=seq(from=-2,to=2,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
#nothing gets tossed
dev.off()

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/betas_maf_sig.pdf",height=4,width=6)
par(mfrow=c(2,3))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
plot(CIGh_abo$beta~(CIGh_abo$af), cex = ifelse(CIGh_abo$p_used < 0.05,1,0.1), #col = rgb(range01(abs(holoLinfoabo$reststate)),0,0), 
	col = c(rgb(0.5,0,0.5),rgb(0,0.5,0))[CIGh_abo$is.plant+1], xlim=c(0,0.5),ylim=c(-4,4)) #limits might need to change
	mtext("Holo-estimated",side=3, line=0.5)
	mtext("Estimated Beta",side=2, line=2,adj = -5)
	mtext("ABO",side=2, line=3)
plot(CIGp_abo$beta~(CIGp_abo$af), cex = ifelse(CIGp_abo$p_used < 0.05,1,0.1), xlim=c(0,0.5),ylim=c(-4,4)) #limits might need to change
	mtext("Plant-estimated",side=3, line=0.5)
plot(CIGm_abo$beta~(CIGm_abo$af), cex = ifelse(CIGm_abo$p_used < 0.05,1,0.1),  xlim=c(0,0.5),ylim=c(-4,4)) 
	mtext("Microbe-estimated",side=3, line=0.5)
plot(CIGh_aba$beta~(CIGh_aba$af), cex = ifelse(CIGh_aba$p_used < 0.05,1,0.1), 
	col = c(rgb(0.5,0,0.5),rgb(0,0.5,0))[CIGh_aba$is.plant+1], 	xlim=c(0,0.5),ylim=c(-4,4)) 
	mtext("ABA",side=2, line=3)
plot(CIGp_aba$beta~(CIGp_aba$af), cex = ifelse(CIGp_aba$p_used < 0.05,1,0.1), xlim=c(0,0.5),ylim=c(-4,4)) 
	mtext("Minor Allele Frequency",side=1, line=2)
plot(CIGm_aba$beta~(CIGm_aba$af), cex = ifelse(CIGm_aba$p_used < 0.05,1,0.1), xlim=c(0,0.5),ylim=c(-4,4)) 
dev.off()
#
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/knownEffs_maf_sig.pdf",height=4,width=6)
par(mfrow=c(2,3))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
plot(holoLinfoabo$reststate~(CIGh_abo$af), cex = ifelse(CIGh_abo$p_used < 0.05,1,0.1), 
	col = c(rgb(0.5,0,0.5),rgb(0,0.5,0))[CIGh_abo$is.plant+1], xlim=c(0,0.5),ylim=c(-1.5,2)) #limits might need to change
	mtext("Holo-estimated",side=3, line=0.5)
	mtext("Known effect",side=2, line=2,adj = -2)
	mtext("ABO",side=2, line=3)
plot(plantLinfoabo$reststate~(CIGp_abo$af), cex = ifelse(CIGp_abo$p_used < 0.05,1,0.1),xlim=c(0,0.5),ylim=c(-1.5,2)) #limits might need to change
	mtext("Plant-estimated",side=3, line=0.5)
plot(micrLinfoabo$reststate~(CIGm_abo$af), cex = ifelse(CIGm_abo$p_used < 0.05,1,0.1), xlim=c(0,0.5),ylim=c(-1.5,2)) 
	mtext("Microbe-estimated",side=3, line=0.5)
plot(holoLinfoaba$reststate~(CIGh_aba$af), cex = ifelse(CIGh_aba$p_used < 0.05,1,0.1), #pch=16,
	col = c(rgb(0.5,0,0.5),rgb(0,0.5,0))[CIGh_aba$is.plant+1], xlim=c(0,0.5),ylim=c(-1.5,2)) 
	mtext("ABA",side=2, line=3)
plot(plantLinfoaba$reststate~(CIGp_aba$af), cex = ifelse(CIGp_aba$p_used < 0.05,1,0.1), xlim=c(0,0.5),ylim=c(-1.5,2)) 
	mtext("Minor Allele Frequency",side=1, line=2)
plot(micrLinfoaba$reststate~(CIGm_aba$af), cex = ifelse(CIGm_aba$p_used < 0.05,1,0.1),  xlim=c(0,0.5),ylim=c(-1.5,2)) 
dev.off()
#
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/knownEffs_beta_sig.pdf",height=4,width=6)
par(mfrow=c(2,3))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
plot(CIGh_abo$beta~(holoLinfoabo$reststate), cex = ifelse(CIGh_abo$p_used < 0.05,1,0.1), 
	col = c(rgb(0.5,0,0.5),rgb(0,0.5,0))[CIGh_abo$is.plant+1], xlim=c(-1.5,2),ylim=c(-4,4)) #limits might need to change
	mtext("Holo-estimated",side=3, line=0.5)
	mtext("Estimated Beta",side=2, line=2,adj = -1)
	mtext("ABO",side=2, line=3)
plot(CIGp_abo$beta~(plantLinfoabo$reststate), cex = ifelse(CIGp_abo$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2),ylim=c(-4,4)) #limits might need to change
	mtext("Plant-estimated",side=3, line=0.5)
plot(CIGm_abo$beta~(micrLinfoabo$reststate), cex = ifelse(CIGm_abo$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2),ylim=c(-4,4)) 
	mtext("Microbe-estimated",side=3, line=0.5)
plot(CIGh_aba$beta~(holoLinfoaba$reststate), cex = ifelse(CIGh_aba$p_used < 0.05,1,0.1), 
	col = c(rgb(0.5,0,0.5),rgb(0,0.5,0))[CIGh_aba$is.plant+1], xlim=c(-1.5,2),ylim=c(-4,4)) 
# 	mtext("Holo-estimated",side=3, line=0.5)
	mtext("ABA",side=2, line=3)
plot(CIGp_aba$beta~(plantLinfoaba$reststate), cex = ifelse(CIGp_aba$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2),ylim=c(-4,4)) 
	mtext("Known effect",side=1, line=2)
# 	mtext("Plant-estimated",side=3, line=0.5)
plot(CIGm_aba$beta~(micrLinfoaba$reststate), cex = ifelse(CIGm_aba$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2),ylim=c(-4,4)) 
dev.off()


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_nokin.pdf",height=8,width=8)
par(mfrow=c(2,3))

image(summarizehabo$conttabprp,xaxt="n",yaxt="n",main="Together")
	mtext("ABO w/o Kinship, interval P",side=2,line=2)
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizehabo$conttab,rep("/",times=4),
		rep(c(sum(holooutabo$known_neg),sum(holooutabo$known_pos)),times=2),
		rep("=",times=4),
		round(as.vector(t(summarizehabo$conttabprp)),digits=2), sep="" ))
image(t(matrix( (summarizemabo$cattots/summarizemabo$catshouldbesums)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Microbe loci")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizemabo$cattots[1:4],rep("/",times=4), summarizemabo$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizemabo$cattots/summarizemabo$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizepabo$cattots/summarizepabo$catshouldbesums)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Plant loci")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizepabo$cattots[5:8],rep("/",times=4), summarizepabo$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizepabo$cattots/summarizepabo$catshouldbesums)[5:8],digits=2), sep="" ))
image(summarizehaba$conttabprp,xaxt="n",yaxt="n",main="Together")
	mtext("ABA w/o Kinship, interval P",side=2,line=2)
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizehaba$conttab,rep("/",times=4),
		rep(c(sum(holooutaba$known_neg),sum(holooutaba$known_pos)),times=2),
		rep("=",times=4),
		round(as.vector(t(summarizehaba$conttabprp)),digits=2), sep="" ))
image(t(matrix( (summarizemaba$cattots/summarizemaba$catshouldbesums)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Microbe loci")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizemaba$cattots[1:4],rep("/",times=4), summarizemaba$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizemaba$cattots/summarizemaba$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizepaba$cattots/summarizepaba$catshouldbesums)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Plant loci")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizepaba$cattots[5:8],rep("/",times=4), summarizepaba$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizepaba$cattots/summarizepaba$catshouldbesums)[5:8],digits=2), sep="" ))

# image(summarizehKabo$conttabprp,xaxt="n",yaxt="n",main="Together")
# 	mtext("w/o Kinship, interval p",side=2,line=2)
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehKabo$conttab,rep("/",times=4),
# 		rep(c(sum(holooutKabo$known_neg),sum(holooutKabo$known_pos)),times=2),
# 		rep("=",times=4),
# 		round(as.vector(t(summarizehKabo$conttabprp)),digits=2), sep="" ))
# image(t(matrix( (summarizehKabo$cattots/summarizehKabo$catshouldbesums)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Microbe loci")
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehKabo$cattots[1:4],rep("/",times=4), summarizehKabo$catshouldbesums[1:4], rep("=",times=4),
# 		round(  (summarizehKabo$cattots/summarizehKabo$catshouldbesums)[1:4],digits=2), sep="" ))
# image(t(matrix( (summarizehKabo$cattots/summarizehKabo$catshouldbesums)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Plant loci")
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehKabo$cattots[5:8],rep("/",times=4), summarizehKabo$catshouldbesums[5:8], rep("=",times=4),
# 		round(  (summarizehKabo$cattots/summarizehKabo$catshouldbesums)[5:8],digits=2), sep="" ))

dev.off()


# pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_ABA.pdf",height=8,width=8)
# par(mfrow=c(2,3))
# 
# image(summarizehaba$conttabprp,xaxt="n",yaxt="n",main="Together")
# 	mtext("w/o Kinship, interval p",side=2,line=2)
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehaba$conttab,rep("/",times=4),
# 		rep(c(sum(holooutaba$known_neg),sum(holooutaba$known_pos)),times=2),
# 		rep("=",times=4),
# 		round(as.vector(t(summarizehaba$conttabprp)),digits=2), sep="" ))
# image(t(matrix( (summarizehaba$cattots/summarizehaba$catshouldbesums)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Microbe loci")
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehaba$cattots[1:4],rep("/",times=4), summarizehaba$catshouldbesums[1:4], rep("=",times=4),
# 		round(  (summarizehaba$cattots/summarizehaba$catshouldbesums)[1:4],digits=2), sep="" ))
# image(t(matrix( (summarizehaba$cattots/summarizehaba$catshouldbesums)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Plant loci")
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehaba$cattots[5:8],rep("/",times=4), summarizehaba$catshouldbesums[5:8], rep("=",times=4),
# 		round(  (summarizehaba$cattots/summarizehaba$catshouldbesums)[5:8],digits=2), sep="" ))
# 
# 
# image(summarizehKaba$conttabprp,xaxt="n",yaxt="n",main="Together")
# 	mtext("w/o Kinship, interval p",side=2,line=2)
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehKaba$conttab,rep("/",times=4),
# 		rep(c(sum(holooutKaba$known_neg),sum(holooutKaba$known_pos)),times=2),
# 		rep("=",times=4),
# 		round(as.vector(t(summarizehKaba$conttabprp)),digits=2), sep="" ))
# image(t(matrix( (summarizehKaba$cattots/summarizehKaba$catshouldbesums)[1:4],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Microbe loci")
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehKaba$cattots[1:4],rep("/",times=4), summarizehKaba$catshouldbesums[1:4], rep("=",times=4),
# 		round(  (summarizehKaba$cattots/summarizehKaba$catshouldbesums)[1:4],digits=2), sep="" ))
# image(t(matrix( (summarizehKaba$cattots/summarizehKaba$catshouldbesums)[5:8],nrow=2,byrow=F)),xaxt="n",yaxt="n",main="Plant loci")
# 	axis(2,at=c(0,1),labels=c("neutral","causal"))
# 	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
# 	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
# 		summarizehKaba$cattots[5:8],rep("/",times=4), summarizehKaba$catshouldbesums[5:8], rep("=",times=4),
# 		round(  (summarizehKaba$cattots/summarizehKaba$catshouldbesums)[5:8],digits=2), sep="" ))
# 
# dev.off()
# 
# 
