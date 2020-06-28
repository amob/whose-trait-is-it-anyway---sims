
addcols <- function(outfile, k=NULL){
	outfile$recodebeta <- sapply(1:nrow(outfile), function(z) ifelse(outfile$allele1[z]==2, outfile$beta[z],-1*outfile$beta[z]) )#allele1 is the minor allele, but states 1 and 2 are derived and ancestral (0)
#	#the other thing that messes with beta estimation is the linkage to another negative allele, especially if the others are more negative, because then the effect of the focal allele must be *relatively* positive. 
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
plantout4a <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4a.assoc.txt",header=T,sep="\t")
plantout4aff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4aff.assoc.txt",header=T,sep="\t")
plantout4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4f.assoc.txt",header=T,sep="\t")
plantout4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/PLANTgemmaABO4fff.assoc.txt",header=T,sep="\t")
microut4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4b.assoc.txt",header=T,sep="\t")
microut4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4bff.assoc.txt",header=T,sep="\t")
microut4a <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4a.assoc.txt",header=T,sep="\t")
microut4aff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4aff.assoc.txt",header=T,sep="\t")
microut4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4f.assoc.txt",header=T,sep="\t")
microut4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/MICRgemmaABO4fff.assoc.txt",header=T,sep="\t")

plantout4b <- addcols(plantout4b)
plantout4bff <- addcols(plantout4bff)
plantout4a <- addcols(plantout4a)
plantout4aff <- addcols(plantout4aff)
plantout4f <- addcols(plantout4f)
plantout4fff <- addcols(plantout4fff)
microut4b <- addcols(microut4b)
microut4bff <- addcols(microut4bff)
microut4a <- addcols(microut4a)
microut4aff <- addcols(microut4aff)
microut4f <- addcols(microut4f)
microut4fff <- addcols(microut4fff)

summarizep4b  <- summarizeSNPcalls(plantout4b)
summarizep4bff  <- summarizeSNPcalls(plantout4bff)
summarizep4a  <- summarizeSNPcalls(plantout4a)
summarizep4aff  <- summarizeSNPcalls(plantout4aff)
summarizep4f  <- summarizeSNPcalls(plantout4f)
summarizep4fff  <- summarizeSNPcalls(plantout4fff)
summarizem4b  <- summarizeSNPcalls(microut4b)
summarizem4bff  <- summarizeSNPcalls(microut4bff)
summarizem4a  <- summarizeSNPcalls(microut4a)
summarizem4aff  <- summarizeSNPcalls(microut4aff)
summarizem4f  <- summarizeSNPcalls(microut4f)
summarizem4fff  <- summarizeSNPcalls(microut4fff)

#verifying that the zerostate is always zero. if not, then have to re-evaluate how deal with recoding beta!!!
plantLinfo4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4b.csv",header=T)
plantLinfo4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4bff.csv",header=T)
plantLinfo4a <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4a.csv",header=T)
plantLinfo4aff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4aff.csv",header=T)
plantLinfo4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4f.csv",header=T)
plantLinfo4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/plantlociABO4fff.csv",header=T)
micrLinfo4b <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4b.csv",header=T)
micrLinfo4bff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4bff.csv",header=T)
micrLinfo4a <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4a.csv",header=T)
micrLinfo4aff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4aff.csv",header=T)
micrLinfo4f <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4f.csv",header=T)
micrLinfo4fff <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/micrlociABO4fff.csv",header=T)

table(plantLinfo4bff$zerostate)
table(plantLinfo4a$zerostate)
table(plantLinfo4aff$zerostate)
table(plantLinfo4f$zerostate)
table(plantLinfo4fff$zerostate)
table(micrLinfo4b$zerostate)
table(micrLinfo4bff$zerostate)
table(micrLinfo4a$zerostate)
table(micrLinfo4aff$zerostate)
table(micrLinfo4f$zerostate)
table(micrLinfo4fff$zerostate)



CIGp_4b <- connectIG(plantLinfo4b,plantout4b)
CIGp_4bff <- connectIG(plantLinfo4bff,plantout4bff)
CIGp_4a <- connectIG(plantLinfo4a,plantout4a)
CIGp_4aff <- connectIG(plantLinfo4aff,plantout4aff)
CIGp_4f <- connectIG(plantLinfo4f,plantout4f)
CIGp_4fff <- connectIG(plantLinfo4fff,plantout4fff)
CIGm_4b <- connectIG(micrLinfo4b,microut4b)
CIGm_4bff <- connectIG(micrLinfo4bff,microut4bff)
CIGm_4a <- connectIG(micrLinfo4a,microut4a)
CIGm_4aff <- connectIG(micrLinfo4aff,microut4aff)
CIGm_4f <- connectIG(micrLinfo4f,microut4f)
CIGm_4fff <- connectIG(micrLinfo4fff,microut4fff)


#visualize what is lost, it is the small effect size loci
pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/dens_loci_tossed_sixdemos.pdf",height=4,width=12)
par(mfrow=c(2,6))
par(oma=c(3,3,2,0))
par(mar=c(2,2,1,1))
hist(plantLinfo4b$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(plantLinfo4b$reststate[!is.na(CIGp_4b$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Frequency",side=2, line=3)
	mtext("Host estimated",side=2, line=2)
hist(plantLinfo4bff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(plantLinfo4bff$reststate[!is.na(CIGp_4bff$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(plantLinfo4a$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(plantLinfo4a$reststate[!is.na(CIGp_4a$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(plantLinfo4aff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(plantLinfo4aff$reststate[!is.na(CIGp_4aff$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(plantLinfo4f$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(plantLinfo4f$reststate[!is.na(CIGp_4f$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(plantLinfo4fff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(plantLinfo4fff$reststate[!is.na(CIGp_4fff$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4b$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(micrLinfo4b$reststate[!is.na(CIGm_4b$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
	mtext("Frequency",side=2, line=3)
	mtext("Microbe estimated",side=2, line=2)
hist(micrLinfo4bff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(micrLinfo4bff$reststate[!is.na(CIGm_4bff$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4a$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(micrLinfo4a$reststate[!is.na(CIGm_4a$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4aff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(micrLinfo4aff$reststate[!is.na(CIGm_4aff$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4f$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(micrLinfo4f$reststate[!is.na(CIGm_4f$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
hist(micrLinfo4fff$reststate,freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),col=rgb(0.9,0.75,0,alpha=0.5),ylim=c(0,60),main="") #ylim=c(0,200))#
	hist(micrLinfo4fff$reststate[!is.na(CIGm_4fff$recodebeta)],freq=T,breaks=seq(from=-2.5,to=2.5,by=0.1),add=T,col=rgb(0,0,0,alpha=0.25))
dev.off()

#to adjust figure ranges if needed
range( c(plantLinfo4b$reststate,plantLinfo4bff$reststate,plantLinfo4a$reststate,plantLinfo4aff$reststate,plantLinfo4f$reststate,plantLinfo4fff$reststate))
range( c(micrLinfo4b$reststate,micrLinfo4bff$reststate,micrLinfo4a$reststate,micrLinfo4aff$reststate,micrLinfo4f$reststate,micrLinfo4fff$reststate))
range( c(CIGp_4b$beta,CIGp_4bff$beta,CIGp_4a$beta,CIGp_4aff$beta,CIGp_4f$beta,CIGp_4fff$beta),na.rm=T)
range( c(CIGm_4b$beta,CIGm_4bff$beta,CIGm_4a$beta,CIGm_4aff$beta,CIGm_4f$beta,CIGm_4fff$beta),na.rm=T)


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/betas_maf_sig_sixdemos.pdf",height=5,width=7)#4, 12
layout(matrix(1:6,ncol=3,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,3,0))
plot(	   CIGp_4b$recodebeta~(CIGp_4b$af), cex = ifelse(CIGp_4b$p_used < 0.05,1.5,0.5), col=rgb(0,0.5,0),  xlim=c(0,0.5),ylim=c(-3,1.5)) #limits might need to change
	points(CIGm_4b$recodebeta~(CIGm_4b$af), cex = ifelse(CIGm_4b$p_used < 0.05,1.5,0.5),  col=rgb(0.5,0,0.5))
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No direct link to microbe fitness",side=3,line=0.5)
	mtext("Estimated Beta",side=2, line=2,adj=-0.5)
	abline(h=0)
plot(      CIGp_4bff$recodebeta~(CIGp_4bff$af), cex = ifelse(CIGp_4bff$p_used < 0.05,1.5,0.5),, col=rgb(0,0.5,0), xlim=c(0,0.5),ylim=c(-3,1.5))
	points(CIGm_4bff$recodebeta~(CIGm_4bff$af), cex = ifelse(CIGm_4bff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	abline(h=0)
plot(	   CIGp_4a$recodebeta~(CIGp_4a$af), cex = ifelse(CIGp_4a$p_used < 0.05,1.5,0.5),col=rgb(0,0.5,0),  xlim=c(0,0.5),ylim=c(-3,1.5)) 
	points(CIGm_4a$recodebeta~(CIGm_4a$af), cex = ifelse(CIGm_4a$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("Equal links to fitness",side=3,line=0.5)
	mtext("Same optima",side=3,line=2)
	abline(h=0)
plot(	   CIGp_4aff$recodebeta~(CIGp_4aff$af), cex = ifelse(CIGp_4aff$p_used < 0.05,1.5,0.5), col=rgb(0,0.5,0),xlim=c(0,0.5),ylim=c(-3,1.5)) 
	points(CIGm_4aff$recodebeta~(CIGm_4aff$af), cex = ifelse(CIGm_4aff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	 mtext("Minor Allele Frequency",side=1, line=2)
	abline(h=0)
plot(	   CIGp_4f$recodebeta~(CIGp_4f$af), cex = ifelse(CIGp_4f$p_used < 0.05,1.5,0.5), col=rgb(0,0.5,0),xlim=c(0,0.5),ylim=c(-3,1.5)) 
	points(CIGm_4f$recodebeta~(CIGm_4f$af), cex = ifelse(CIGm_4f$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("Equal links to fitness",side=3,line=0.5)
	mtext("Different optima",side=3,line=2)
	abline(h=0)
plot(	   CIGp_4fff$recodebeta~(CIGp_4fff$af), cex = ifelse(CIGp_4fff$p_used < 0.05,1.5,0.5), col=rgb(0,0.5,0),xlim=c(0,0.5),ylim=c(-3,1.5))
	points(CIGm_4fff$recodebeta~(CIGm_4fff$af), cex = ifelse(CIGm_4fff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	abline(h=0)
legend(0.3,-1.7,c("Host","Microbe"),fill=c(rgb(0,0.5,0),rgb(0.5,0,0.5)),bty="n")
dev.off()
#

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/knownEffs_maf_sig_sixdemos.pdf",height=5,width=7)
layout(matrix(1:6,ncol=3,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,3,0))
plot(		plantLinfo4b$reststate~(CIGp_4b$af), cex = ifelse(CIGp_4b$p_used < 0.05,1.5,0.5), xlim=c(0,0.5),ylim=c(-0.5,0.7), col=rgb(0,0.5,0)) 
	points( micrLinfo4b$reststate~(CIGm_4b$af), cex = ifelse(CIGm_4b$p_used < 0.05,1.5,0.5), col=rgb(0.5,0,0.5) )
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No direct link to microbe fitness",side=3,line=0.5)
	mtext("Known effect",side=2, line=2,adj=-0.5)
	abline(h=0)
plot(		plantLinfo4bff$reststate~(CIGp_4bff$af), cex = ifelse(CIGp_4bff$p_used < 0.05,1.5,0.5), xlim=c(0,0.5),ylim=c(-0.5,0.7), col=rgb(0,0.5,0))
	points(	micrLinfo4bff$reststate~(CIGm_4bff$af), cex = ifelse(CIGm_4bff$p_used < 0.05,1.5,0.5), col=rgb(0.5,0,0.5) )
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	abline(h=0)
plot(		plantLinfo4a$reststate~(CIGp_4a$af), cex = ifelse(CIGp_4a$p_used < 0.05,1.5,0.5), xlim=c(0,0.5),ylim=c(-0.5,0.7), col=rgb(0,0.5,0))
	points(	micrLinfo4a$reststate~(CIGm_4a$af), cex = ifelse(CIGm_4a$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("Equal links to fitness",side=3,line=0.5)
	mtext("Same optima",side=3,line=2)
	abline(h=0)
plot(		plantLinfo4aff$reststate~(CIGp_4aff$af), cex = ifelse(CIGp_4aff$p_used < 0.05,1.5,0.5), xlim=c(0,0.5),ylim=c(-0.5,0.7), col=rgb(0,0.5,0))
	points(	micrLinfo4aff$reststate~(CIGm_4aff$af), cex = ifelse(CIGm_4aff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	 mtext("Minor Allele Frequency",side=1, line=2)
	abline(h=0)
plot(		plantLinfo4f$reststate~(CIGp_4f$af), cex = ifelse(CIGp_4f$p_used < 0.05,1.5,0.5), xlim=c(0,0.5),ylim=c(-0.5,0.7), col=rgb(0,0.5,0))
	points(	micrLinfo4f$reststate~(CIGm_4f$af), cex = ifelse(CIGm_4f$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("Equal links to fitness",side=3,line=0.5)
	mtext("Different optima",side=3,line=2)
	abline(h=0)
plot(		plantLinfo4fff$reststate~(CIGp_4fff$af), cex = ifelse(CIGp_4fff$p_used < 0.05,1.5,0.5), xlim=c(0,0.5),ylim=c(-0.5,0.7), col=rgb(0,0.5,0))
	points(	micrLinfo4fff$reststate~(CIGm_4fff$af), cex = ifelse(CIGm_4fff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	abline(h=0)
legend(0.3,-0.3,c("Host","Microbe"),fill=c(rgb(0,0.5,0),rgb(0.5,0,0.5)),bty="n")
dev.off()
#


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/knownEffs_beta_sig_sixdemos.pdf",height=5,width=7)
layout(matrix(1:6,ncol=3,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,3,0))
plot(		CIGp_4b$recodebeta~(plantLinfo4b$reststate), cex = ifelse(CIGp_4b$p_used < 0.05,1.5,0.5), xlim=c(-0.5,0.7),ylim=c(-3,1.5),col=rgb(0,0.5,0)) 
	points(	CIGm_4b$recodebeta~(micrLinfo4b$reststate), cex = ifelse(CIGm_4b$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No direct link to microbe fitness",side=3,line=0.5)
	mtext("Estimated Beta",side=2, line=2,adj=-0.5)
	text(0.5,y=-2, bquote(rho==.(round(cor(CIGp_4b$recodebeta,plantLinfo4b$reststate,use="complete.obs"),digits=2)) ),col=rgb(0,0.5,0)) 
	text(0.5,y=-2.3, bquote(rho==.(round(cor(CIGm_4b$recodebeta, micrLinfo4b$reststate,use="complete.obs"),digits=2)) ),col=rgb(0.5,0,0.5) )
plot(		CIGp_4bff$recodebeta~(plantLinfo4bff$reststate), cex = ifelse(CIGp_4bff$p_used < 0.05,1.5,0.5),  xlim=c(-0.5,0.7),ylim=c(-3,1.5),col=rgb(0,0.5,0)) 
	points(	CIGm_4bff$recodebeta~(micrLinfo4bff$reststate), cex = ifelse(CIGm_4bff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
	text(0.5,y=-2, bquote(rho==.(round(cor(CIGp_4bff$recodebeta,plantLinfo4bff$reststate,use="complete.obs"),digits=2)) ),col=rgb(0,0.5,0)) 
	text(0.5,y=-2.3, bquote(rho==.(round(cor(CIGm_4bff$recodebeta, micrLinfo4bff$reststate,use="complete.obs"),digits=2)) ),col=rgb(0.5,0,0.5) )
plot(		CIGp_4a$recodebeta~(plantLinfo4a$reststate), cex = ifelse(CIGp_4a$p_used < 0.05,1.5,0.5),  xlim=c(-0.5,0.7),ylim=c(-3,1.5),col=rgb(0,0.5,0)) 
	points(	CIGm_4a$recodebeta~(micrLinfo4a$reststate), cex = ifelse(CIGm_4a$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("Equal links to fitness",side=3,line=0.5)
	mtext("Same optima",side=3,line=2)
	text(0.5,y=-2, bquote(rho==.(round(cor(CIGp_4a$recodebeta,plantLinfo4a$reststate,use="complete.obs"),digits=2)) ),col=rgb(0,0.5,0)) 
	text(0.5,y=-2.3, bquote(rho==.(round(cor(CIGm_4a$recodebeta, micrLinfo4a$reststate,use="complete.obs"),digits=2)) ),col=rgb(0.5,0,0.5) )
plot(		CIGp_4aff$recodebeta~(plantLinfo4aff$reststate), cex = ifelse(CIGp_4aff$p_used < 0.05,1.5,0.5),  xlim=c(-0.5,0.7),ylim=c(-3,1.5),col=rgb(0,0.5,0)) 
	points(	CIGm_4aff$recodebeta~(micrLinfo4aff$reststate), cex = ifelse(CIGm_4aff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	 mtext("Known Effect",side=1, line=2)
	text(0.5,y=-2, bquote(rho==.(round(cor(CIGp_4aff$recodebeta,plantLinfo4aff$reststate,use="complete.obs"),digits=2)) ),col=rgb(0,0.5,0)) 
	text(0.5,y=-2.3, bquote(rho==.(round(cor(CIGm_4aff$recodebeta, micrLinfo4aff$reststate,use="complete.obs"),digits=2)) ),col=rgb(0.5,0,0.5) )
plot(		CIGp_4f$recodebeta~(plantLinfo4f$reststate), cex = ifelse(CIGp_4f$p_used < 0.05,1.5,0.5),  xlim=c(-0.5,0.7),ylim=c(-3,1.5),col=rgb(0,0.5,0)) 
	points(	CIGm_4f$recodebeta~(micrLinfo4f$reststate), cex = ifelse(CIGm_4f$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	mtext("Equal links to fitness",side=3,line=0.5)
	mtext("Different optima",side=3,line=2)
	text(0.5,y=-2, bquote(rho==.(round(cor(CIGp_4f$recodebeta,plantLinfo4f$reststate,use="complete.obs"),digits=2)) ),col=rgb(0,0.5,0)) 
	text(0.5,y=-2.3, bquote(rho==.(round(cor(CIGm_4f$recodebeta, micrLinfo4f$reststate,use="complete.obs"),digits=2)) ),col=rgb(0.5,0,0.5) )
plot(		CIGp_4fff$recodebeta~(plantLinfo4fff$reststate), cex = ifelse(CIGp_4fff$p_used < 0.05,1.5,0.5),  xlim=c(-0.5,0.7),ylim=c(-3,1.5),col=rgb(0,0.5,0)) 
	points(	CIGm_4fff$recodebeta~(micrLinfo4fff$reststate), cex = ifelse(CIGm_4fff$p_used < 0.05,1.5,0.5),col=rgb(0.5,0,0.5) )
	text(0.5,y=-2, bquote(rho==.(round(cor(CIGp_4fff$recodebeta,plantLinfo4fff$reststate,use="complete.obs"),digits=2)) ),col=rgb(0,0.5,0)) 
	text(0.5,y=-2.3, bquote(rho==.(round(cor(CIGm_4fff$recodebeta, micrLinfo4fff$reststate,use="complete.obs"),digits=2)) ),col=rgb(0.5,0,0.5) )
legend(-0.3,-1.8,c("Host","Microbe"),fill=c(rgb(0,0.5,0),rgb(0.5,0,0.5)),bty="n")
dev.off()


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_nokin_sixdemos_plant.pdf",height=5,width=7)
layout(matrix(1:6,ncol=3,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,3,0))
image(t(matrix( (summarizep4b$cattots/summarizep4b$catshouldbesums)[5:8],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4b$cattots[5:8],rep("/",times=4), summarizep4b$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizep4b$cattots/summarizep4b$catshouldbesums)[5:8],digits=2), sep="" ))
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No link to microbe fitness",side=3)
image(t(matrix( (summarizep4bff$cattots/summarizep4bff$catshouldbesums)[5:8],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4bff$cattots[5:8],rep("/",times=4), summarizep4bff$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizep4bff$cattots/summarizep4bff$catshouldbesums)[5:8],digits=2), sep="" ))
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
image(t(matrix( (summarizep4a$cattots/summarizep4a$catshouldbesums)[5:8],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4a$cattots[5:8],rep("/",times=4), summarizep4a$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizep4a$cattots/summarizep4a$catshouldbesums)[5:8],digits=2), sep="" ))
	mtext("Equal links to fitness",side=3)
	mtext("Same optima",side=3,line=2)
image(t(matrix( (summarizep4aff$cattots/summarizep4aff$catshouldbesums)[5:8],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4aff$cattots[5:8],rep("/",times=4), summarizep4aff$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizep4aff$cattots/summarizep4aff$catshouldbesums)[5:8],digits=2), sep="" ))
image(t(matrix( (summarizep4f$cattots/summarizep4f$catshouldbesums)[5:8],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4f$cattots[5:8],rep("/",times=4), summarizep4f$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizep4f$cattots/summarizep4f$catshouldbesums)[5:8],digits=2), sep="" ))
	mtext("Equal links to fitness",side=3)
	mtext("Different optima",side=3,line=2)
image(t(matrix( (summarizep4fff$cattots/summarizep4fff$catshouldbesums)[5:8],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizep4fff$cattots[5:8],rep("/",times=4), summarizep4fff$catshouldbesums[5:8], rep("=",times=4),
		round(  (summarizep4fff$cattots/summarizep4fff$catshouldbesums)[5:8],digits=2), sep="" ))
dev.off()

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest_nokin_sixdemos_micr.pdf",height=5,width=7)
layout(matrix(1:6,ncol=3,byrow=F))
par(mar=c(1.5,3,1,1))
par(oma=c(2,2,3,0))
image(t(matrix( (summarizem4b$cattots/summarizem4b$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4b$cattots[1:4],rep("/",times=4), summarizem4b$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4b$cattots/summarizem4b$catshouldbesums)[1:4],digits=2), sep="" ))
	mtext("No fitness feedback",side=2, line=3.5,cex=1.25)
	mtext("No link to microbe fitness",side=3)
image(t(matrix( (summarizem4bff$cattots/summarizem4bff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4bff$cattots[1:4],rep("/",times=4), summarizem4bff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4bff$cattots/summarizem4bff$catshouldbesums)[1:4],digits=2), sep="" ))
	mtext("+ fitness feedback",side=2, line=3.5,cex=1.25)
image(t(matrix( (summarizem4a$cattots/summarizem4a$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4a$cattots[1:4],rep("/",times=4), summarizem4a$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4a$cattots/summarizem4a$catshouldbesums)[1:4],digits=2), sep="" ))
	mtext("Equal links to fitness",side=3)
	mtext("Same optima",side=3,line=2)
image(t(matrix( (summarizem4aff$cattots/summarizem4aff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4aff$cattots[1:4],rep("/",times=4), summarizem4aff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4aff$cattots/summarizem4aff$catshouldbesums)[1:4],digits=2), sep="" ))
image(t(matrix( (summarizem4f$cattots/summarizem4f$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4f$cattots[1:4],rep("/",times=4), summarizem4f$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4f$cattots/summarizem4f$catshouldbesums)[1:4],digits=2), sep="" ))
	mtext("Equal links to fitness",side=3)
	mtext("Different optima",side=3,line=2)
image(t(matrix( (summarizem4fff$cattots/summarizem4fff$catshouldbesums)[1:4],nrow=2,byrow=F)),zlim=c(0,1),xaxt="n",yaxt="n",main="")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		summarizem4fff$cattots[1:4],rep("/",times=4), summarizem4fff$catshouldbesums[1:4], rep("=",times=4),
		round(  (summarizem4fff$cattots/summarizem4fff$catshouldbesums)[1:4],digits=2), sep="" ))
dev.off()

