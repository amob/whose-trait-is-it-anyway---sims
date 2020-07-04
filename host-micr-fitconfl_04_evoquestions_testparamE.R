
simsens <-read.csv("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters.csv",header=T)
reps = 5 #MUST SET REPLICATES

pars <- simsens[,1:17]
resps <- simsens[,18:26]
#not all parameters were modified
mpars <- c(1,3,4,4, 9,10,12:16)#which then is manipulated one after the other, all but zoM here.
#4 is twice because loci were manipulated separately then together.

#####PASTED FROM 02 series. CHECK ACCURACY if you change things
#copy paste design of 02 series that generated the files
basevals <- c(2000,2000, 20,40, 3,3, 3,2,      0.75,0.75, 1000,       25, 0.0005,      0.2,    0.6,0.6,   0.1)
#sim.cotrait(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,
# 		popsz.v <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000,10000)  #other interesting ranges here, commented out because do not match current settings.
		popsz.v <- c(100, 500, 900, 1300, 1700, 2100, 2500, 2900, 3300, 4100) 
# 		nlocP.v = c(2, 2, 4, 8, 16, 32, 64, 128, 256, 516),#  base set to 20 --
		nloc.v = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)# 
		w.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)#
		Lambda.v <- seq(from = 35, to = 17, by =-2) #base 25
		mutprb.v <- c( 0.0001, 0.0002, 0.0003,0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001 )
# 		mutprb.v <- c(0.0000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.1) #base 0.0001
		prbHorz.v <- seq(from =0, to =1, length.out=10)  
		alpha.v <- seq(from = 0.0, to =0.9, by =0.1) #base 0.6	


parm <- data.frame(matrix(rep(basevals,times=111),nrow=111,byrow=T)) #81 is the base case
parm[1:10,1] <- popsz.v
parm[1:10,2] <- popsz.v
parm[11:20,3] <- nloc.v #P
parm[21:30,4] <- nloc.v #M
parm[31:40,3] <- nloc.v #tog
parm[31:40,4] <- nloc.v*2 #tog
parm[41:50,9] <- w.v #P
parm[51:60,10] <- w.v #M
parm[61:70,12] <- Lambda.v
parm[71:80,13] <- mutprb.v
parm[81:90,14] <- prbHorz.v
parm[91:100,15] <- alpha.v #repeated 2x! once for plants, once for microbes
parm[101:110,16] <- alpha.v
#111st and 222nd rows are the base state
parm2 <- rbind(parm,parm)
parm2[112:222,8] <- 3 #change to fitness agreement; now both have optima at 3
#


params <- data.frame( # for axis variables, thus repeated rows
		popsz.v = popsz.v,
		nloc.vP = nloc.v,
		nloc.vM = nloc.v,
		nloc.v = nloc.v,
		w.vP = w.v,
		w.vM = w.v,
		Lambda.v =Lambda.v,
		mutprb.v = mutprb.v,
		prbHorz.v = prbHorz.v,
		alpha.vP= alpha.v,
		alpha.vM= alpha.v
)


scenarios <- rep(1:nrow(parm2),each=reps)
#remove rows are subtracted
#there are 6 of 1110 that hav NA values in 3 columns. it seems to be because there is 0 additive genetic variance, so it can't calculate a proportion or a fitness correlation
means.s <- lapply( 1:ncol(resps), function(r) matrix(
					sapply((1:nrow(parm2) )[-c(111,222)], function(scenario) mean(resps[ which(scenarios==scenario) , r],na.rm=T) ), 
					ncol= 10, byrow = T)   )
vars.s <- lapply( 1:ncol(resps), function(r) matrix(
					sapply((1:nrow(parm2) )[-c(111,222)], function(scenario) var(resps[ which(scenarios==scenario) , r],na.rm=T) ), 
					ncol= 10, byrow = T)   )

means <- lapply( means.s, function(r) cbind(r[1:length(mpars),],r[(length(mpars)+1):(length(mpars)*2),])  )
vars <- lapply( vars.s, function(r) cbind(r[1:length(mpars),],r[(length(mpars)+1):(length(mpars)*2),])  )




mains <- c("Fitness correlation","Vh/(Vh+Vm)",expression(Plant~V[A]~(V[H])),"Vm/(Vh+Vm)",expression(Microbe~V[A]~(V[H])),
			expression(Mean~paste("|","D","|",sep="")~from~z[opt][H]),"Microbe avg |D| from zopt","Host terminal slope", "Microbe terminal slope")
ylabs <- c(expression(alpha[M]),expression(alpha[H]),expression(P[hrz]),
			expression(P[mu]),
			expression(lambda),expression(omega[M]),expression(omega[H]),
			expression(L),expression(L[M]),expression(L[H]),expression(N))
y2labs <- c(paste(params[c(1,10),11],collapse="-"), paste(params[c(1,10),10],collapse="-"), 
			paste(params[c(1,10),9],collapse="-"), paste(params[c(1,10),8],collapse="-"), paste(round(params[c(1,10),7],digits=2),collapse="-"),
			paste(params[c(1,10),6],collapse="-"),  
			paste(params[c(1,10),5],collapse="-"), paste(params[c(1,10),4],collapse="-"), paste(params[c(1,10),3],collapse="-"),
			paste(params[c(1,10),2],collapse="-"), paste(params[c(1,10),1],collapse="-") )

rb <- colorRampPalette(c( rgb(1,0,0), rgb(1,1,1), rgb(0,0,1) ))
wb <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,1) ))
wb2 <- colorRampPalette(c( rgb(1,1,1), rgb(0,0,0.75,alpha=0.75), rgb(0,0,1) ))

indices <- c(1,3,5,6) #parts of tmp2 or responses that we are including.

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters.pdf",width=9,height=2)
layout(matrix(c(1:5),ncol=5,byrow=T),widths=c(2,2,2,2,1.1))
par(mar=c(1,0,3,1))
par(oma=c(2,3,0,6))
for(i in indices){
	if(i%in%(8:9)){
				intend <- max(abs(range(means[8:9])))
				zlims <- c(-1*intend, intend)#artificial range. if change sims this may need to change
				cols <- rb(50)
			} else if(i == 1){
				zlims <- c(-1,1)
				cols <- rb(50)
			} else{
				zlims <- range(means[[i]])
				cols <- wb(50)
			}
	image(t(means[[i]]),main="",xaxt="n",yaxt="n",zlim=zlims,col=cols)
	abline(v=0.5)
	mtext(mains[i],side=3,line=1,cex=0.75)
	axis(side=1,at = c(0.22,0.73),lab=c(expression("Z"[opt][H]>"Z"[opt][M]),expression("Z"[opt][H]*'='*"Z"[opt][M]) ) ) 
	if(i==1){axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)}
	mtext(paste(round(zlims,digits=3),collapse=" to "),side=3,line=0.2,cex=0.5)
}
image(t(as.matrix(t(params)/colSums(params))),main="",xaxt="n",yaxt="n")
# 	axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][H]>"z"[opt][M]),expression("z"[opt][H]*'='*"z"[opt][M]) ) ) 
# 	axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)
	axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)
	mtext("Parameter Values",side=3,line=0.5,cex=0.75)
dev.off()



pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/sens_reps_finalparameters_vars.pdf",width=9,height=2)
layout(matrix(c(1:5),ncol=5,byrow=T),widths=c(2,2,2,2,1.1))
par(mar=c(1,0,3,1))
par(oma=c(2,3,0,6))
for(i in indices){
	zlims <- range(vars[[i]])
	cols <- wb(50)
	image(t(vars[[i]]),main="",xaxt="n",yaxt="n",zlim=zlims,col=cols)
	abline(v=0.5)
	mtext(mains[i],side=3,line=1,cex=0.75)
	axis(side=1,at = c(0.22,0.73),lab=c(expression("Z"[opt][H]>"Z"[opt][M]),expression("Z"[opt][H]*'='*"Z"[opt][M]) ) ) 
#	if(i%in%c(7:9)){axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][P]>"z"[opt][M]),expression("z"[opt][P]*'='*"z"[opt][M]) ) ) }
	if(i==1){axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)}
#	if(i%in%c(6,9)){axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)}
	mtext(paste(round(zlims,digits=4),collapse=" to "),side=3,line=0.2,cex=0.5)
}
image(t(as.matrix(t(params)/colSums(params))),main="",xaxt="n",yaxt="n")
# 	axis(side=1,at = c(0.22,0.73),lab=c(expression("z"[opt][H]>"z"[opt][M]),expression("z"[opt][H]*'='*"z"[opt][M]) ) ) 
	#axis(side=2,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(ylabs),las=2)
	axis(side=4,at = seq(from=0,to=1,length.out=length(mpars)),lab=rev(y2labs),las=2)
	mtext("Parameter Values",side=3,line=0.5,cex=0.75)
dev.off()

