
library(SDMTools)#legend.gradient


bufferX <- function(x,p) { 
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	}

range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

# simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_feedbackparameters.csv",header=T)

simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_feedbackparameters_holdplant.csv",header=T)

pars <- simsens[,1:17] #
resps <- simsens[,18:26]
mpars <- c(7,8,9,10,15,16)#now its only  -- ws, zs, and pfs

#1 - 2800
# however, each 10 rows together are the same simulation parameters.
parby10 <- as.data.frame(t(sapply(1:280, function(p) colMeans(  pars[c(1:10+10*(p-1)),] ) )))
respsby10 <- as.data.frame(t(sapply(1:280, function(p) colMeans(  resps[c(1:10+10*(p-1)),] ) )))
respsby10V <- as.data.frame(t(sapply(1:280, function(p) sapply(1:ncol(resps), function(colmn) var(resps[c(1:10+10*(p-1)),colmn] ) )  )))

respsby10$ltVp <- log(respsby10$tVp)
respsby10$ltVm <- log(respsby10$tVm)

rb <- colorRampPalette(c( rgb(1,0,0), rgb(1,1,1), rgb(0,0,1) ))
wb <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,1) ))


ylabv <- c( bquote("Fitness"~rho),"pVp",bquote("V"[P]),"pVm",bquote("V"[M]),
			expression(paste("|","D","|",sep="")~from~Z[opt][P]),expression(paste("|","D","|",sep="")~from~Z[opt][P]),
			bquote("Slope Z"[P]),"tcoefM", 
			bquote("log V "[P]),bquote("log V"[M]))
			
		# c("Fitness correlation",expression(Plant~V[A]~(V[P])),
# 			expression(Mean~paste("|","D","|",sep="")~from~Z[opt][P]),"Microbe avg |D| from Zopt","Plant terminal slope", "Microbe terminal slope")


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/feedbacks_mostres.pdf",height=6,width=6)
#rofint <- c(1,3,5,6,8,9)
rofint <- c(10,11,1,8,6)
par(mfrow=c(length(rofint),length(unique(parby10$zoM)))) 
par(mar =c(0,0,0,0))
par(oma = c(4,10,4,1))
for(i in rofint){
	for(j in unique(parby10$zoM)){
		if(i ==1){
			ranges <- c(-1.1,1.1)
			} else if(i%in%c(6,7)){
			ranges <- c(-0.3,3.3)
			} else if(i == 5){ #high variance in this one, I think its a "real" phenomena, just messes with visualization!
			ranges <- c(0,0.05) #currently plot log instead to avoid issue
			} else {ranges <- bufferX(range(respsby10[,i]),0.1)} #[parby10$zoM==j]; keeping on same scale means only 1 y axis
		if(i == rofint[1]){mains <- bquote("Z"[opt][M]==.(j))} else{mains <- ""}
		if(i == rofint[1]){mains2 <- bquote("Z"[opt][P]==5)} else{mains2 <- ""}
		if(j==  unique(parby10$zoM)[1]){ylabs <- ylabv[i]} else{ ylabs <- ""}
		plot(respsby10[,i][parby10$zoM==j]~I(1-parby10$pfM[parby10$zoM==j]),ylim=ranges,
	  		col=rgb(range01(parby10$wM[parby10$zoM==j]),0,0,alpha=0.5),pch=16,
	  		ylab="",xlab="",main ="", 
	  		yaxt = ifelse(j== unique(parby10$zoM)[1],"s","n"), xaxt = ifelse(i==rofint[length(rofint)], "s","n") )
	  	mtext(ylabs,side=2,line=2,las=1)
	  	mtext(mains,side=3,line=1.5)
	  	mtext(mains2,side=3,line=0.25)
		if(i%in%c(1,8,9)) {abline(h=0,lty=2)}
		abline(v=0.6,lty=3,col=rgb(0,0,0,alpha=0.5))
	}
}
legend.gradient(cbind(c(0.6,0.7,0.7,0.6),c(0.5,0.5,2.5,2.5)),cols=rgb(range01(unique(parby10$wM)),0,0,alpha=0.75),title=bquote(omega[M]),limits=c("0.25","1.75"))
mtext(bquote(1-alpha[M]),side=1,line=3,adj=2)
dev.off()

