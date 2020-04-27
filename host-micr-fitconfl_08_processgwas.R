

holoout <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/tmp_HOLOgemma.assoc.txt",header=T,sep="\t")
holooutK <- read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/tmp_HOLOgemmaK.assoc.txt",header=T,sep="\t")


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


startcol <- grep("false_pos",colnames(holoout))
startcolK <- grep("false_pos",colnames(holooutK))

conttabprp <- matrix(
	colSums(holoout[,startcol:(startcol+3)]) /
		rep(c(sum(holoout$known_neg),sum(holoout$known_pos)),times=2),
	byrow=T,ncol=2
		)
			
conttabKprp <- matrix(
	colSums(holooutK[,startcolK:(startcolK+3)]) /
		rep(c(sum(holooutK$known_neg),sum(holooutK$known_pos)),times=2),
	byrow=T,ncol=2
		) # false pos, true pos, true neg, false neg
		
conttabK <- 	colSums(holooutK[,startcolK:(startcolK+3)])

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/quickGWAStest.pdf")
image(conttabKprp,xaxt="n",yaxt="n")
	axis(2,at=c(0,1),labels=c("neutral","causal"))
	axis(1,at=c(0,1),labels=c("p < 0.05","p > 0.05"))
	text(x=c(0,0,1,1),y=c(0,1,0,1),labels=paste(
		conttabK,rep("/",times=4),
		rep(c(sum(holoout$known_neg),sum(holoout$known_pos)),times=2),
		rep("=",times=4),
		round(as.vector(t(conttabKprp)),digits=2), sep="" ))
dev.off()