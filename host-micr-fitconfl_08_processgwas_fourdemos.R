
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


plantout4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4b.assoc.txt",header=T,sep="\t")
plantout4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4bff.assoc.txt",header=T,sep="\t")
plantout4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4f.assoc.txt",header=T,sep="\t")
plantout4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4fff.assoc.txt",header=T,sep="\t")
microut4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4b.assoc.txt",header=T,sep="\t")
microut4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4bff.assoc.txt",header=T,sep="\t")
microut4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4f.assoc.txt",header=T,sep="\t")
microut4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4fff.assoc.txt",header=T,sep="\t")

plantout4b <- addcols(plantout4b)
plantout4bff <- addcols(plantout4bff)
plantout4f <- addcols(plantout4f)
plantout4fff <- addcols(plantout4fff)
microut4b <- addcols(microut4b)
microut4bff <- addcols(microut4bff)
microut4f <- addcols(microut4f)
microut4fff <- addcols(microut4fff)

summarizep4b  <- summarizeSNPcalls(plantout4b)
summarizep4bff  <- summarizeSNPcalls(plantout4bff)
summarizep4f  <- summarizeSNPcalls(plantout4f)
summarizep4fff  <- summarizeSNPcalls(plantout4fff)
summarizem4b  <- summarizeSNPcalls(microut4b)
summarizem4bff  <- summarizeSNPcalls(microut4bff)
summarizem4f  <- summarizeSNPcalls(microut4f)
summarizem4fff  <- summarizeSNPcalls(microut4fff)

plantLinfo4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4b.csv",header=T)
plantLinfo4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4bff.csv",header=T)
plantLinfo4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4f.csv",header=T)
plantLinfo4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4fff.csv",header=T)
micrLinfo4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4b.csv",header=T)
micrLinfo4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4bff.csv",header=T)
micrLinfo4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4f.csv",header=T)
micrLinfo4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4fff.csv",header=T)


CIGp_4b <- connectIG(plantLinfo4b,plantout4b)
CIGp_4bff <- connectIG(plantLinfo4bff,plantout4bff)
CIGp_4f <- connectIG(plantLinfo4f,plantout4f)
CIGp_4fff <- connectIG(plantLinfo4fff,plantout4fff)
CIGm_4b <- connectIG(micrLinfo4b,microut4b)
CIGm_4bff <- connectIG(micrLinfo4bff,microut4bff)
CIGm_4f <- connectIG(micrLinfo4f,microut4f)
CIGm_4fff <- connectIG(micrLinfo4fff,microut4fff)


#visualize what is lost, it is the small effect size loci
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/dens_loci_tossed_fourdemos.pdf",height=4,width=8)
par(mfrow=c(2,4))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
hist(plantLinfo4b$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(plantLinfo4b$reststate[!is.na(CIGp_4b$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Frequency",side=2, line=3)
	mtext("Plant estimated",side=2, line=2)
hist(plantLinfo4bff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(plantLinfo4bff$reststate[!is.na(CIGp_4bff$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(plantLinfo4f$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(plantLinfo4f$reststate[!is.na(CIGp_4f$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(plantLinfo4fff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(plantLinfo4fff$reststate[!is.na(CIGp_4fff$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4b$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(micrLinfo4b$reststate[!is.na(CIGm_4b$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Frequency",side=2, line=3)
	mtext("Plant estimated",side=2, line=2)
hist(micrLinfo4bff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(micrLinfo4bff$reststate[!is.na(CIGm_4bff$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4f$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(micrLinfo4f$reststate[!is.na(CIGm_4f$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4fff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,100),main="") #ylim=c(0,200))#
	hist(micrLinfo4fff$reststate[!is.na(CIGm_4fff$beta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
dev.off()

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/betas_maf_sig_fourdemos.pdf",height=4,width=8)
# par(mfrow=c(2,3))
par(mfrow=c(2,4))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
plot(CIGp_4b$beta~(CIGp_4b$af), cex = ifelse(CIGp_4b$p_used < 0.05,1,0.1), xlim=c(0,0.2),ylim=c(-4,2)) #limits might need to change
	mtext("Plant genome",side=2, line=3)
	mtext("Estimated Beta",side=2, line=2)
plot(CIGp_4bff$beta~(CIGp_4bff$af), cex = ifelse(CIGp_4bff$p_used < 0.05,1,0.1), xlim=c(0,0.2),ylim=c(-4,2))
plot(CIGp_4f$beta~(CIGp_4f$af), cex = ifelse(CIGp_4f$p_used < 0.05,1,0.1), xlim=c(0,0.2),ylim=c(-4,2)) 
plot(CIGp_4fff$beta~(CIGp_4fff$af), cex = ifelse(CIGp_4fff$p_used < 0.05,1,0.1), xlim=c(0,0.2),ylim=c(-4,2))

plot(CIGm_4b$beta~(CIGm_4b$af), cex = ifelse(CIGm_4b$p_used < 0.05,1,0.1),  xlim=c(0,0.2),ylim=c(-4,2)) 
	mtext("Microbe genome",side=2, line=2)
	mtext("Minor Allele Frequency",side=1, line=2, adj = 3)
plot(CIGm_4bff$beta~(CIGm_4bff$af), cex = ifelse(CIGm_4bff$p_used < 0.05,1,0.1),  xlim=c(0,0.2),ylim=c(-4,2)) 
plot(CIGm_4f$beta~(CIGm_4f$af), cex = ifelse(CIGm_4f$p_used < 0.05,1,0.1),  xlim=c(0,0.2),ylim=c(-4,2)) 
plot(CIGm_4fff$beta~(CIGm_4fff$af), cex = ifelse(CIGm_4fff$p_used < 0.05,1,0.1),  xlim=c(0,0.2),ylim=c(-4,2)) 
dev.off()
#
# pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/knownEffs_maf_sig_fourdemos.pdf",height=4,width=8)
# par(mfrow=c(1,2))
# par(oma=c(3,3,2,0))
# par(mar=c(2,2,1,1))
# plot(plantLinfo4b$reststate~(CIGp_4b$af), cex = ifelse(CIGp_4b$p_used < 0.05,1,0.1), 
# 	 xlim=c(0,0.5),ylim=c(-1.5,2.5)) #limits might need to change
# 	mtext("Plant-estimated",side=3, line=0.5)
# 	mtext("Known effect",side=2, line=2,adj = -2)
# plot(micrLinf4b$reststate~(CIGm_4b$af), cex = ifelse(CIGm_4b$p_used < 0.05,1,0.1),  xlim=c(0,0.5),ylim=c(-1.5,2.5)) 
# dev.off()
#
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/knownEffs_beta_sig_fourdemos.pdf",height=4,width=8)
par(mfrow=c(2,4))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
plot(CIGp_4b$beta~(plantLinfo4b$reststate), cex = ifelse(CIGp_4b$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) #limits might need to change
	mtext("Plant genome",side=3, line=0.5)
	mtext("Estimated Beta",side=2, line=2)
plot(CIGp_4bff$beta~(plantLinfo4bff$reststate), cex = ifelse(CIGp_4bff$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
plot(CIGp_4f$beta~(plantLinfo4f$reststate), cex = ifelse(CIGp_4f$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
plot(CIGp_4fff$beta~(plantLinfo4f$reststate), cex = ifelse(CIGp_4fff$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
plot(CIGm_4b$beta~(micrLinfo4b$reststate), cex = ifelse(CIGm_4b$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
	mtext("Microbe genome",side=3, line=0.5)
	mtext("Known Effect",side=1, line=2, adj = -1.5)
plot(CIGm_4bff$beta~(micrLinfo4bff$reststate), cex = ifelse(CIGm_4bff$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
plot(CIGm_4f$beta~(micrLinfo4f$reststate), cex = ifelse(CIGm_4f$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
plot(CIGm_4fff$beta~(micrLinfo4f$reststate), cex = ifelse(CIGm_4fff$p_used < 0.05,1,0.1), 
	 xlim=c(-1.5,2.5),ylim=c(-4,2)) 
dev.off()



pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_nokin_fourdemos_plant.pdf",height=4,width=4)
par(mfrow=c(2,2))
image(t(matrix( (summarizep4b$cattots/summarizep4b$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4b$cattots[1:4],rep("/",times=4), summarizep4b$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizep4b$cattots/summarizep4b$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizep4f$cattots/summarizep4f$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4f$cattots[1:4],rep("/",times=4), summarizep4f$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizep4f$cattots/summarizep4f$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizep4bff$cattots/summarizep4bff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4bff$cattots[1:4],rep("/",times=4), summarizep4bff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizep4bff$cattots/summarizep4bff$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizep4fff$cattots/summarizep4fff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4fff$cattots[1:4],rep("/",times=4), summarizep4fff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizep4fff$cattots/summarizep4fff$catshouldbesums)[1:4],digits=2), sep="" ))
dev.off()

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_nokin_fourdemos_micr.pdf",height=4,width=4)
par(mfrow=c(2,2))
image(t(matrix( (summarizem4b$cattots/summarizem4b$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4b$cattots[1:4],rep("/",times=4), summarizem4b$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4b$cattots/summarizem4b$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizem4f$cattots/summarizem4f$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4f$cattots[1:4],rep("/",times=4), summarizem4f$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4f$cattots/summarizem4f$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizem4bff$cattots/summarizem4bff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4bff$cattots[1:4],rep("/",times=4), summarizem4bff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4bff$cattots/summarizem4bff$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizem4fff$cattots/summarizem4fff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4fff$cattots[1:4],rep("/",times=4), summarizem4fff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4fff$cattots/summarizem4fff$catshouldbesums)[1:4],digits=2), sep="" ))
dev.off()

