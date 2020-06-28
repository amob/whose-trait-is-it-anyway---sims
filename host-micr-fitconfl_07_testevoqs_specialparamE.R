
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

reps <- 5 ## MUST SET NUMBER OF REPS.

# simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_feedbackparameters.csv",header=T)

simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_feedbackparameters_holdplant.csv",header=T)

pars <- simsens[,1:17] #
resps <- simsens[,18:26]
mpars <- c(7,8,9,10,15,16)#now its only  -- ws, zs, and pfs

#1 - 2800
# however, each 10 rows together are the same simulation parameters.
parbyreps <- as.data.frame(t(sapply(1:280, function(p) colMeans(  pars[c(1:reps+reps*(p-1)),] ) )))
respsbyreps <- as.data.frame(t(sapply(1:280, function(p) colMeans(  resps[c(1:reps+reps*(p-1)),] ) )))
respsbyrepsV <- as.data.frame(t(sapply(1:280, function(p) sapply(1:ncol(resps), function(colmn) var(resps[c(1:reps+reps*(p-1)),colmn] ) )  )))

respsbyreps$ltVp <- log(respsbyreps$tVp)
respsbyreps$ltVm <- log(respsbyreps$tVm)

rb <- colorRampPalette(c( rgb(1,0,0), rgb(1,1,1), rgb(0,0,1) ))
wb <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,1) ))


ylabv <- c( bquote("Fitness"~rho),"pVh",bquote("V"[H]),"pVm",bquote("V"[M]),
			expression(paste("|","D","|",sep="")~from~Z[opt][H]),expression(paste("|","D","|",sep="")~from~Z[opt][H]),
			bquote("Slope Z"[H]),"tcoefM", 
			bquote("log V "[H]),bquote("log V"[M]))
			
		# c("Fitness correlation",expression(Plant~V[A]~(V[P])),
# 			expression(Mean~paste("|","D","|",sep="")~from~Z[opt][P]),"Microbe avg |D| from Zopt","Plant terminal slope", "Microbe terminal slope")


pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/feedbacks_mostres.pdf",height=8,width=7)
#rofint <- c(1,3,5,6,8,9)
rofint <- c(10,11,1,8,6)
par(mfrow=c(length(rofint),length(unique(parbyreps$zoP)))) 
par(mar =c(0,0,0,0))
par(oma = c(4,10,4,1))
for(i in rofint){
	for(j in unique(parbyreps$zoP)){
		if(i ==1){
			ranges <- c(-1.1,1.1)
			} else if(i%in%c(6,7)){
			ranges <- c(-0.3,3.3)
			} else if(i == 5){ #high variance in this one, I think its a "real" phenomena, just messes with visualization!
			ranges <- c(0,0.05) #currently plot log instead to avoid issue
			} else {ranges <- bufferX(range(respsbyreps[,i]),0.1)} #[parby10$zoM==j]; keeping on same scale means only 1 y axis
		if(i == rofint[1]){mains <- bquote("Z"[opt][H]==.(j))} else{mains <- ""}
		if(i == rofint[1]){mains2 <- bquote("Z"[opt][M]==2)} else{mains2 <- ""}
		if(j==  unique(parbyreps$zoP)[1]){ylabs <- ylabv[i]} else{ ylabs <- ""}
		plot(respsbyreps[,i][parbyreps$zoP==j]~I(1-parbyreps$pfM[parbyreps$zoP==j]),ylim=ranges,
# 		plot(respsby10[,i][parby10$zoM==j]~I(parby10$pfM[parby10$zoM==j]),ylim=ranges,
	  		col=rgb(range01(parbyreps$wM[parbyreps$zoP==j]),0,0,alpha=0.5),pch=16,
	  		ylab="",xlab="",main ="", 
	  		yaxt = ifelse(j== unique(parbyreps$zoP)[1],"s","n"), xaxt = ifelse(i==rofint[length(rofint)], "s","n") )
	  	mtext(ylabs,side=2,line=2,las=1)
	  	mtext(mains,side=3,line=1.5)
	  	if(j == 4 & i==tail(rofint,1)) {mtext(bquote("1-alpha"[M]),side=1,line=2.5,adj=-0.3)}
# 	  	if(j == 4 & i==tail(rofint,1)) {mtext(bquote("alpha"[M]),side=1,line=1.5,adj=-0.5)}
	  	mtext(mains2,side=3,line=0.25)
		if(i%in%c(1,8,9)) {abline(h=0,lty=2)}
		abline(v=0.6,lty=3,col=rgb(0,0,0,alpha=0.5))
	}
}
legend.gradient(cbind(c(0.6,0.7,0.7,0.6),c(0.7,0.7,2.7,2.7)),cols=rgb(range01(unique(parbyreps$wM)),0,0,alpha=0.75),title=bquote(omega[M]),limits=c("0.25","1.75"))
mtext(bquote(1-alpha[M]),side=1,line=3,adj=2)
dev.off()

