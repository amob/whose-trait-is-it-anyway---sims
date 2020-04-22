
simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_feedbackparameters.csv",header=T)

pars <- simsens[,1:19] #
resps <- simsens[,20:28]
mpars <- c(7,8,9,10,18,19)#now its only 7,8,9,10, 18,19 -- ws, zs, and pfs

#1 - 3240
# however, each 10 rows together are the same simulation parameters.




parby10 <- as.data.frame(t(sapply(1:324, function(p) colMeans(  pars[c(1:10+10*(p-1)),] ) )))
respsby10 <- as.data.frame(t(sapply(1:324, function(p) colMeans(  resps[c(1:10+10*(p-1)),] ) )))
respsby10V <- as.data.frame(t(sapply(1:324, function(p) sapply(1:ncol(resps), function(colmn) var(resps[c(1:10+10*(p-1)),colmn] ) )  )))

#going to aim for 3 panel figures, one for when zoM>zoP, when zoP>zoM, and when zoP=zoM

par10noconf <- parby10[which(parby10$zoP==parby10$zoM),]
par10pm <- parby10[which(parby10$zoP > parby10$zoM),]
par10mp <- parby10[which(parby10$zoP < parby10$zoM),]
res10noconf <- respsby10[which(parby10$zoP==parby10$zoM),]
res10pm <- respsby10[which(parby10$zoP > parby10$zoM),]
res10mp <- respsby10[which(parby10$zoP < parby10$zoM),]


rb <- colorRampPalette(c( rgb(1,0,0), rgb(1,1,1), rgb(0,0,1) ))
wb <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,1) ))
pg <- colorRampPalette(c( rgb(0.5,0,0.5), rgb(1,1,1), rgb(0,0.5,0) ))
wg <- colorRampPalette(c( rgb(1,1,1), rgb(0,0.5,0) ))
wbk <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,0) ))


image(t(respsby10[,c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(respsby10[ which(parby10$wM==1 & parby10$wP==1)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))


par(mfrow=c(2,2))
image(t(res10noconf[ which(par10noconf$wM==1 & par10noconf$wP==1)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(res10noconf[ which(par10noconf$wM==1 & par10noconf$wP==0.25)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(res10noconf[ which(par10noconf$wM==0.25 & par10noconf$wP==1)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(res10noconf[ which(par10noconf$wM==0.25 & par10noconf$wP==0.25)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))

par(mfrow=c(2,2))
image(t(res10pm[ which(par10pm$wM==1 & par10pm$wP==1)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(res10pm[ which(par10pm$wM==1 & par10pm$wP==0.25)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(res10pm[ which(par10pm$wM==0.25 & par10pm$wP==1)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))
image(t(res10pm[ which(par10pm$wM==0.25 & par10pm$wP==0.25)  ,  c(1,2,4)]),zlim=c(-1,1),col=rb(100))


###taking variables 1 at a time

ylabs <- c("fitness correlation","Vp/(Vp+Vm)","Plant Va (Vp)","Vm/(Vp+Vm)","Microbe Va (Vm)",
			"Plant avg |D| from Zopt","Microbe avg |D| from Zopt","plant terminal slope", "microbe terminal slope")
tops <- c("zP=3, zM=1","zP=5, zM=1","zP=5, zM=3",expression(alpha[M]~and~alpha[P]))
#check.
mptops <- c("zM=3, zP=1","zM=5, zP=1","zM=5, zP=3",expression(alpha[M]~and~alpha[P]))
nctops <- c("zP=5, zM=5","zP=3, zM=3","zP=1, zM=1",expression(alpha[M]~and~alpha[P]) )

xlabs <- c(expression(omega[P]==1~omega[M]==1),expression(omega[P]==1~omega[M]==0.25),
			expression(omega[P]==0.25~omega[M]==1),expression(omega[P]==0.25~omega[M]==0.25))

xlabs2 <- c("1,1","1,0.25","0.25,1","0.25,0.25" )

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/feedbacks_ZpgreaterthanZm.pdf",width=4,height=6)
par(mfrow=c(5,2))
par(mar=c(0,0,2,0))
par(oma=c(5,3,2,5))
for(i in c(1:3,6,8)){
		if(i == 1){
				zlims <- c(-1,1)
				cols <- pg(50)
			} else if(i%in%c(2,4)){
				zlims <- c(0,1)
				cols <- wg(50)
			} else if(i%in%(8:9)){
				intend <- max(abs(range(res10pm[,8:9])))
				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
				cols <- pg(50)
			} else{
				zlims <- range(res10pm[,i])
				cols <- wg(50)
			}
	tmp <- cbind(res10pm[which(par10pm$wP==1 & par10pm$wM==1),i],
		res10pm[which(par10pm$wP==1 & par10pm$wM==0.25),i],
		res10pm[which(par10pm$wP==0.25 & par10pm$wM==1),i],
		res10pm[which(par10pm$wP==0.25 & par10pm$wM==0.25),i])
	partmp <- (par10pm[which(par10pm$wP==1 & par10pm$wM==1),7:8])
	image(t(tmp[1:9,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
#		if(i==8){axis(side=1,at = seq(from=0,to=1,length.out=4),labels =xlabs,las=2 ) }
		if(i==8){
			mtext(expression(omega[P]~omega[M]),side =1, padj=-7,las=2,line=0)
			axis(side=1,at = c(seq(from=0,to=1,length.out=4)),labels =xlabs2,las=2 ) }
		mtext(ylabs[i],side=3,line=0.5)
		if(i==1){mtext(tops[1],side=3,line=2)}
	image(t(tmp[10:18,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
#		if(i==8){axis(side=1,at = seq(from=0,to=1,length.out=4),labels =xlabs,las=2 ) }
		if(i==8){axis(side=1,at = seq(from=0,to=1,length.out=4),labels =xlabs2,las=2 ) }
		mtext(paste(round(zlims,digits=2),collapse=" to "),side=3,line=0.5)
		if(i==1){
			mtext(tops[2],side=3,line=2)
			mtext(expression(alpha[P]~alpha[M]),side =3, adj=1.5,line=0)}
		axis(side=4,at=seq(from=0,to=1,length.out=9), las=2,
			labels=c("0.9,0.9","0.4,0.9","0.1,0.9","0.9,0.4","0.4,0.4","0.1,0.4","0.1,0.9","0.1,0.4","0.1,0.1"))
# 	image(t(tmp[19:27,]), zlim=zlims,col=cols)
# 	if(i==1){mtext(tops[3],side=3,line=2)}
# 	image(t(cbind(par10pm[which(par10pm$wP==1 & par10pm$wM==1)[1:9],17],
# 		par10pm[which(par10pm$wP==1 & par10pm$wM==1)[1:9],18])),
# 		zlim=c(0,1),col=wbk(100),xaxt="n",yaxt="n")
# 		if(i==1){mtext(tops[4],side=3,line=2)}
# 		if(i==1){mtext("0 to 1",side=3,line=0.5)}
}
dev.off()

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/feedbacks_ZpequalsZm.pdf",width=4,height=6)
par(mfrow=c(5,2))
par(mar=c(0,0,2,0))
par(oma=c(5,3,2,5))
for(i in c(1:3,6,8)){
		if(i == 1){
				zlims <- c(-1,1)
				cols <- pg(50)
			} else if(i%in%c(2,4)){
				zlims <- c(0,1)
				cols <- wg(50)
			} else if(i%in%(8:9)){
				intend <- max(abs(range(res10noconf[,8:9])))
				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
				cols <- pg(50)
			} else{
				zlims <- range(res10noconf[,i])
				cols <- wg(50)
			}
	tmp <- cbind(res10noconf[which(par10noconf$wP==1 & par10noconf$wM==1),i],
		res10noconf[which(par10noconf$wP==1 & par10noconf$wM==0.25),i],
		res10noconf[which(par10noconf$wP==0.25 & par10noconf$wM==1),i],
		res10noconf[which(par10noconf$wP==0.25 & par10noconf$wM==0.25),i])
	image(t(tmp[1:9,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
		mtext(ylabs[i],side=3,line=0.5)
		if(i==1){mtext(nctops[1],side=3,line=2)}
		if(i==8){
			mtext(expression(omega[P]~omega[M]),side =1, padj=-7,las=2,line=0)
			axis(side=1,at = c(seq(from=0,to=1,length.out=4)),labels =xlabs2,las=2 ) 
		}
	image(t(tmp[10:18,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
		if(i==1){
			mtext(nctops[2],side=3,line=2)
			mtext(expression(alpha[P]~alpha[M]),side =3, adj=1.5,line=0)}
		if(i==8){axis(side=1,at = c(seq(from=0,to=1,length.out=4)),labels =xlabs2,las=2 ) }
		axis(side=4,at=seq(from=0,to=1,length.out=9), las=2,
			labels=c("0.9,0.9","0.4,0.9","0.1,0.9","0.9,0.4","0.4,0.4","0.1,0.4","0.1,0.9","0.1,0.4","0.1,0.1"))
		mtext(paste(round(zlims,digits=2),collapse=" to "),side=3,line=0.5)
}
dev.off()



pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/feedbacks_tog.pdf",width=5,height=6)
par(mfrow=c(5,4))
par(mar=c(0,0,2,0))
par(oma=c(5,3,2,5))
for(i in c(1:3,6,8)){
		if(i == 1){
				zlims <- c(-1,1)
				cols <- pg(50)
			} else if(i%in%c(2,4)){
				zlims <- c(0,1)
				cols <- wg(50)
			} else if(i%in%(8:9)){
				intend <- max(abs(range(rbind(res10pm[,8:9],res10noconf[,8:9]))))
				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
				cols <- pg(50)
			} else{
				zlims <- range(c(res10pm[,i],res10noconf[,i]))
				cols <- wg(50)
			}
	tmp <- cbind(res10pm[which(par10pm$wP==1 & par10pm$wM==1),i],
		res10pm[which(par10pm$wP==1 & par10pm$wM==0.25),i],
		res10pm[which(par10pm$wP==0.25 & par10pm$wM==1),i],
		res10pm[which(par10pm$wP==0.25 & par10pm$wM==0.25),i])
	tmp2 <- cbind(res10noconf[which(par10noconf$wP==1 & par10noconf$wM==1),i],
		res10noconf[which(par10noconf$wP==1 & par10noconf$wM==0.25),i],
		res10noconf[which(par10noconf$wP==0.25 & par10noconf$wM==1),i],
		res10noconf[which(par10noconf$wP==0.25 & par10noconf$wM==0.25),i])
	image(t(tmp[1:9,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
#		if(i==8){axis(side=1,at = seq(from=0,to=1,length.out=4),labels =xlabs,las=2 ) }
		if(i==8){
			mtext(expression(omega[P]~omega[M]),side =1, padj=-5,las=2,line=0)
			axis(side=1,at = c(seq(from=0,to=1,length.out=4)),labels =xlabs2,las=2 ) }
	#	mtext(ylabs[i],side=3,line=0.5)
		if(i==1){mtext(tops[1],side=3,line=2)}
	image(t(tmp[10:18,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
#		if(i==8){axis(side=1,at = seq(from=0,to=1,length.out=4),labels =xlabs,las=2 ) }
		if(i==8){axis(side=1,at = seq(from=0,to=1,length.out=4),labels =xlabs2,las=2 ) }
		mtext(paste(ylabs[i],", range: ",paste(round(zlims,digits=2),collapse=" to "),sep=""),side=3,line=0.5)
		if(i==1){
			mtext(tops[2],side=3,line=2)
		}
			
	image(t(tmp2[1:9,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
		if(i==1){mtext(nctops[1],side=3,line=2)}
		if(i==8){axis(side=1,at = c(seq(from=0,to=1,length.out=4)),labels =xlabs2,las=2 ) }
	image(t(tmp2[10:18,]), zlim=zlims,col=cols,xaxt="n",yaxt="n")
		if(i==1){
			mtext(nctops[2],side=3,line=2)
			mtext(expression(alpha[P]~alpha[M]),side =3, adj=1.75,line=0)}
		if(i==8){axis(side=1,at = c(seq(from=0,to=1,length.out=4)),labels =xlabs2,las=2 ) }
		axis(side=4,at=seq(from=0,to=1,length.out=9), las=2,
			labels=c("0.9,0.9","0.4,0.9","0.1,0.9","0.9,0.4","0.4,0.4","0.1,0.4","0.1,0.9","0.1,0.4","0.1,0.1"))
}
dev.off()

####mostly reversed for zm >zp
# par(mfrow=c(5,4))
# par(mar=c(1,1,4,1))
# par(oma=c(2,2,0,0))
# for(i in c(1:3,6,8)){
# 		if(i == 1){
# 				zlims <- c(-1,1)
# 				cols <- pg(50)
# 			} else if(i%in%c(2,4)){
# 				zlims <- c(0,1)
# 				cols <- wg(50)
# 			} else if(i%in%(8:9)){
# 				intend <- max(abs(range(res10mp[,8:9])))
# 				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
# 				cols <- pg(50)
# 			} else{
# 				zlims <- range(res10mp[,i])
# 				cols <- wg(50)
# 			}
# 	tmp <- cbind(res10mp[which(par10mp$wP==1 & par10mp$wM==1),i],
# 		res10mp[which(par10mp$wP==1 & par10mp$wM==0.25),i],
# 		res10mp[which(par10mp$wP==0.25 & par10mp$wM==1),i],
# 		res10mp[which(par10mp$wP==0.25 & par10mp$wM==0.25),i])
# 	image(t(tmp[1:9,]), zlim=zlims,col=cols)
# 	mtext(ylabs[i],side=3,line=0.5)
# 	if(i==1){mtext(mptops[1],side=3,line=2)}
# 	image(t(tmp[10:18,]), zlim=zlims,col=cols)
# 	if(i==1){mtext(mptops[2],side=3,line=2)}
# 	image(t(tmp[19:27,]), zlim=zlims,col=cols)
# 	if(i==1){mtext(mptops[3],side=3,line=2)}
# 	image(t(cbind(par10mp[which(par10mp$wP==1 & par10mp$wM==1)[1:9],17],
# 		par10mp[which(par10mp$wP==1 & par10mp$wM==1)[1:9],18])),
# 		zlim=c(0,1),col=wg(100))
# 		if(i==1){mtext(mptops[4],side=3,line=2)}
# 
# }

#####EDITED FOR now up to here

# ylabs <- c(expression(alpha[M]),expression(alpha[P]),expression(P[hrz]),
# 			expression(theta[M]),expression(theta[P]),expression(P[mu]),
# 			expression(lambda),expression(omega[M]),expression(omega[P]),
# 			expression(L),expression(N))
# y2labs <- c( paste(pfP.v,collapse=","), paste(pfP.v,collapse=","), paste(round(prbHorz.v,digits=2),collapse=","),
# 			paste(fiterrM.v,collapse=","), paste(fiterrP.v,collapse=","), paste(mutprb.v,collapse=","),  
# 			paste(Lambda.v,collapse=","), paste(wM.v,collapse=","), paste(wP.v,collapse=","),
# 			paste(nloc.v,collapse=","), paste(popsz.v,collapse=",") )
# 
