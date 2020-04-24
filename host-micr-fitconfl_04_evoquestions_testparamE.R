
simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters.csv",header=T)

pars <- simsens[,1:17]
resps <- simsens[,18:26]
#not all pars were modified
mpars <- c(1,3,9,10,12:16)#which then is manipulated one after the other, all but zoM here.

#1:90, 91, 92:181, 182
# 1 thru 900 are part of first set of sims, then 901-910 are base, 911-1810 are second set and 1811-1820 are second base with matched zopt
# however, each 10 rows together are the same simulation parameters.
91:100; 101:110  1101
tmp <- lapply( 1:ncol(resps), function(r) sapply( 1:length(mpars), function(p) resps[ c( 1:100+(100*(p-1)) , 911:1010+(100*(p-1)) ) , r]) )
#each result of the sapply goes into a column

tmp2<- lapply(1:length(tmp), function(r) t(sapply(seq(from=1, to =nrow(tmp[[r]]),by=10), function(startrow) colMeans(tmp[[r]][startrow:(startrow+9),]))  ) )
#each resul of the sapply goes into a column, but we want it to go into a row

tmp2var<- lapply(1:length(tmp), function(r) t(sapply(seq(from=1, to =nrow(tmp[[r]]),by=10), function(startrow) sapply(1:length(mpars), function(clmn)  var(tmp[[r]][startrow:(startrow+9),clmn])  ))  ) )
#each resul of the sapply goes into a column, but we want it to go into a row




#####PASTED FROM 02 series. CHECK ACCURACY
		wP.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
		wM.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75?
		popsz.v <- c(50, 75, 125, 150, 175, 200, 225, 250, 275, 300) #note turning up M without N is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
		nloc.v <- c(50, 75, 125, 150, 175, 200, 225, 250, 275, 300)# multiply by 2 for microbes, base set to 100
		prbHorz.v <- seq(from =0, to =1, length.out=10)
		Lambda.v <- seq(from = 41, to = 23, by =-2) #base 30
		mutprb.v <- seq( from = 0.0001,to= 0.001 , length.out=10) #base 0.0005
		pfP.v <- seq(from = 0.0, to =0.95, by =0.1) #base 0.6	



par(mfrow=c(3,3))
par(mar=c(1,2,3,1))
par(oma=c(4,4,0,1))
for(i in 1:9){
	image(tmp[[i]],main="",xaxt="n",yaxt="n")
	abline(v=0.5)
	mtext(colnames(resps)[i],side=3,line=1)
	if(i%in%c(7:9)){axis(side=1,at = c(0.225,0.725),lab=c("zoP>zoM","zoP=zoM"))}
	if(i%in%c(1,4,7)){axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=colnames(pars)[mpars],las=2)}
}

mains <- c("fitness correlation","Vp/(Vp+Vm)","Plant Va (Vp)","Vm/(Vp+Vm)","Microbe Va (Vm)",
			"Plant avg |D| from Zopt","Microbe avg |D| from Zopt","plant terminal slope", "microbe terminal slope")

ylabs <- c(expression(alpha[M]),expression(alpha[P]),expression(P[hrz]),
			expression(theta[M]),expression(theta[P]),expression(P[mu]),
			expression(lambda),expression(omega[M]),expression(omega[P]),
			expression(L),expression(N))
y2labs <- c( paste(pfP.v,collapse=","), paste(pfP.v,collapse=","), paste(round(prbHorz.v,digits=2),collapse=","),
			paste(mutprb.v,collapse=","),  
			paste(Lambda.v,collapse=","), paste(wM.v,collapse=","), paste(wP.v,collapse=","),
			paste(nloc.v,collapse=","), paste(popsz.v,collapse=",") )


rb <- colorRampPalette(c( rgb(1,0,0), rgb(1,1,1), rgb(0,0,1) ))
wb <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,1) ))
pg <- colorRampPalette(c( rgb(0.5,0,0.5), rgb(1,1,1), rgb(0,0.5,0) ))
wg <- colorRampPalette(c( rgb(1,1,1), rgb(0,0.5,0) ))

pdf("~/Dropbox/host microbe trait evo and gwas/sens_reps_finalparameters.pdf",width=10,height=6)
par(mfrow=c(3,3))
par(mar=c(1,3,3,1))
par(oma=c(2,1,0,26))
for(i in 1:9){
	if(i == 1){
				zlims <- c(-1,1)
				cols <- pg(50)
			} else if(i%in%c(2,4,6,7)){
				zlims <- c(0,1)
				cols <- wg(50)
			} else if(i%in%(8:9)){
				intend <- max(abs(range(tmp2[8:9])))
				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
# 				zlims <- c(-0.03,0.03)#artificial range. if change sims this may need to change
				cols <- pg(50)
			} else{
				zlims <- range(tmp2[[i]])
				cols <- wg(50)
			}
	image(tmp2[[i]],main="",xaxt="n",yaxt="n",zlim=zlims,col=cols)
	abline(v=0.5)
	mtext(mains[i],side=3,line=1,cex=0.75)
	if(i%in%c(7:9)){axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) }
	axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)
	if(i%in%c(3,6,9)){axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)}
	mtext(paste(round(zlims,digits=3),collapse=" to "),side=3,line=0.2,cex=0.5)
}
dev.off()


pdf("~/Dropbox/host microbe trait evo and gwas/sens_reps_finalparameters_vars.pdf",width=10,height=6)
par(mfrow=c(3,3))
par(mar=c(1,3,3,1))
par(oma=c(2,1,0,26))
for(i in 1:9){
# 	if(i == 1){
# 				zlims <- c(-1,1)
# 				cols <- pg(50)
# 			} else if(i%in%c(2,4,6,7)){
# 				zlims <- c(0,1)
# 				cols <- wg(50)
# 			} else if(i%in%(8:9)){
# 				zlims <- c(-0.03,0.03)#artificial range. if change sims this may need to change
# 				cols <- pg(50)
# 			} else{
# 				zlims <- range(tmp2var[[i]])
# 				cols <- wg(50)
# 			}
	zlims <- range(tmp2var[[i]])
	cols <- wg(50)
	image(tmp2var[[i]],main="",xaxt="n",yaxt="n",zlim=zlims,col=cols)
# 	image(tmp2var[[i]],main="",xaxt="n",yaxt="n",col=wg(50))
	abline(v=0.5)
	mtext(mains[i],side=3,line=1,cex=0.75)
	if(i%in%c(7:9)){axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) }
	axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)
	if(i%in%c(3,6,9)){axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)}
	mtext(paste(round(zlims,digits=5),collapse=" to "),side=3,line=0.2,cex=0.5)
}
dev.off()
