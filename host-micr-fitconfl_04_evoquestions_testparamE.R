
simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters.csv",header=T)

pars <- simsens[,1:17]
resps <- simsens[,18:26]
#not all pars were modified
mpars <- c(1,3,9,10,12:16)#which then is manipulated one after the other, all but zoM here.

#1:90, 91, 92:181, 182
# 1 thru 900 are part of first set of sims, then 901-910 are base, 911-1810 are second set and 1811-1820 are second base with matched zopt
# however, each 10 rows together are the same simulation parameters.
#91:100; 101:110  1101
tmp <- lapply( 1:ncol(resps), function(r) sapply( 1:length(mpars), function(p) resps[ c( 1:100+(100*(p-1)) , 911:1010+(100*(p-1)) ) , r]) )
#each result of the sapply goes into a column

tmp2<- lapply(1:length(tmp), function(r) t(sapply(seq(from=1, to =nrow(tmp[[r]]),by=10), function(startrow) colMeans(tmp[[r]][startrow:(startrow+9),]))  ) )
#each resul of the sapply goes into a column, but we want it to go into a row

tmp2var<- lapply(1:length(tmp), function(r) t(sapply(seq(from=1, to =nrow(tmp[[r]]),by=10), function(startrow) sapply(1:length(mpars), function(clmn)  var(tmp[[r]][startrow:(startrow+9),clmn])  ))  ) )
#each resul of the sapply goes into a column, but we want it to go into a row




#####PASTED FROM 02 series. CHECK ACCURACY
params <- data.frame(
		popsz.v = c(50, 75, 125, 150, 175, 200, 225, 250, 275, 300), #note turning up M without N is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
		nloc.v = c(50, 75, 125, 150, 175, 200, 225, 250, 275, 300),# multiply by 2 for microbes, base set to 100
		wP.v = c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5),#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75?
		wM.v = c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5, 5),#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
		Lambda.v = seq(from = 41, to = 23, by =-2), #base 30
		mutprb.v = seq( from = 0.0001,to= 0.001 , length.out=10), #base 0.0005
		prbHorz.v = seq(from =0, to =1, length.out=10),
		alphaP.v = seq(from = 0.0, to =0.95, by =0.1), #base 0.6	
		alphaM.v = seq(from = 0.0, to =0.95, by =0.1) #base 0.6	
)


mains <- c("Fitness correlation","Vp/(Vp+Vm)",expression(Plant~V[A]~(V[P])),"Vm/(Vp+Vm)",expression(Microbe~V[A]~(V[M])),
			expression(Mean~paste("|","D","|",sep="")~from~Z[opt][P]),"Microbe avg |D| from Zopt","Plant terminal slope", "Microbe terminal slope")
# mains <- c("Plant Va (Vp)","Microbe Va (Vm)",
# 			"Plant avg |D| from Zopt","Microbe avg |D| from Zopt",
# 			"plant terminal slope", "microbe terminal slope")

ylabs <- c(expression(alpha[M]),expression(alpha[P]),expression(P[hrz]),
			expression(P[mu]),
			expression(lambda),expression(omega[M]),expression(omega[P]),
			expression(L[P]),expression(N))
# y2labs <- c( paste(pfP.v,collapse=","), paste(pfP.v,collapse=","), paste(round(prbHorz.v,digits=2),collapse=","),
# 			paste(mutprb.v,collapse=","),  
# 			paste(Lambda.v,collapse=","), paste(wM.v,collapse=","), paste(wP.v,collapse=","),
# 			paste(nloc.v,collapse=","), paste(popsz.v,collapse=",") )
y2labs <- c(paste(params[c(1,10),9],collapse="-"), paste(params[c(1,10),8],collapse="-"), paste(round(params[c(1,10),7],digits=2),collapse="-"),
			paste(params[c(1,10),6],collapse="-"),  
			paste(params[c(1,10),5],collapse="-"), paste(params[c(1,10),4],collapse="-"), paste(params[c(1,10),3],collapse="-"),
			paste(params[c(1,10),2],collapse="-"), paste(params[c(1,10),1],collapse="-") )


rb <- colorRampPalette(c( rgb(1,0,0), rgb(1,1,1), rgb(0,0,1) ))
wb <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,1) ))

indices <- c(1,3,5,6) #parts of tmp2 or responses that we are including.

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters.pdf",width=9,height=2)
layout(matrix(c(1:5),ncol=5,byrow=T),widths=c(2,2,2,2,1.1))
par(mar=c(1,0,3,1))
par(oma=c(2,3,0,6))
for(i in indices){
	if(i%in%(8:9)){
				intend <- max(abs(range(tmp2[8:9])))
				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
				cols <- rb(50)
			} else if(i == 1){
				zlims <- c(-1,1)
				cols <- rb(50)
			} else{
				zlims <- range(tmp2[[i]])
				cols <- wb(50)
			}
	image(tmp2[[i]],main="",xaxt="n",yaxt="n",zlim=zlims,col=cols)
	abline(v=0.5)
	mtext(mains[i],side=3,line=1,cex=0.75)
	axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) 
#	if(i%in%c(7:9)){axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) }
	if(i==1){axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)}
#	if(i%in%c(6,9)){axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)}
	mtext(paste(round(zlims,digits=3),collapse=" to "),side=3,line=0.2,cex=0.5)
}
image(t(as.matrix(t(params)/colSums(params))),main="",xaxt="n",yaxt="n")
	axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) 
	#axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)
	axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)
	mtext("Parameter Values",side=3,line=0.5,cex=0.75)
dev.off()



pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters_vars.pdf",width=9,height=2)
layout(matrix(c(1:5),ncol=5,byrow=T),widths=c(2,2,2,2,1.1))
par(mar=c(1,0,3,1))
par(oma=c(2,3,0,6))
for(i in indices){
	zlims <- range(tmp2var[[i]])
	cols <- wb(50)
	image(tmp2var[[i]],main="",xaxt="n",yaxt="n",zlim=zlims,col=cols)
	abline(v=0.5)
	mtext(mains[i],side=3,line=1,cex=0.75)
	axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) 
#	if(i%in%c(7:9)){axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) }
	if(i==1){axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)}
#	if(i%in%c(6,9)){axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)}
	mtext(paste(round(zlims,digits=4),collapse=" to "),side=3,line=0.2,cex=0.5)
}
image(t(as.matrix(t(params)/colSums(params))),main="",xaxt="n",yaxt="n")
	axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) 
	#axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)
	axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)
	mtext("Parameter Values",side=3,line=0.5,cex=0.75)
dev.off()

