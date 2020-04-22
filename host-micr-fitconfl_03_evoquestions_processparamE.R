#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to evaluate simulation outputs parameter ranges.
##


source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 


#how do the following parameters change ans to above:
	#strength of link between trait and fitness; wP and wM, lower is stronger link -- with lower links, PFF becomes more important
		wP.v <- c(0.1, 0.15, 0.25, 0.5, 0.75, 1.25, 1.5, 2, 2.5, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
		wM.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75?
	#drift and selection within each species; # uhhh. wP and wM vs fiterrP and fiterrM
		#N and s
#		fiterrP.v <- c(0.0001, 0.0005, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
#		fiterrM.v <- c(0.0001, 0.0005, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		fiterrP.v <- c( 0.001, 0.005, 0.01, 0.05, 0.01, 0.05, 0.1,0.5,1,1.5) #base set both low, to 0.001
		fiterrM.v <- c( 0.001, 0.005, 0.01, 0.05, 0.01, 0.05, 0.1,0.5,1,1.5) #base set both low, to 0.001
		popsz.v <- c(50, 75, 125, 150, 175, 200, 225, 250, 275, 300) #note turning up M without N is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
		nloc.v <- c(50, 75, 125, 150, 175, 200, 225, 250, 275, 300)# multiply by 2 for microbes, base set to 100
#		nloc.v <- c(10, 25, 50, 75, 125, 150, 200, 250, 300, 500) # multiply by 2 for microbes, base set to 100
		prbHorz.v <- seq(from =0, to =0.70, length.out=10)
	# mutational inputs;  Lambda
		Lambda.v <- seq(from = 41, to = 23, by =-2) #base 30
#		mutprb.v <- c(0.00001, 0.000025, 0.00005, 0.000075, 0.0001, 0.00025, 0.00075, 0.001 ,0.0015, 0.002)
		mutprb.v <- seq( from = 0.0001,to= 0.001 , length.out=10) #base 0.0005
	#symmetries in conflict between partners (are these poss? I think just change pff for one only)
		pfP.v <- seq(from = 0.0, to =0.95, by =0.1) #base 0.6	
	#fitness optima?
		#maybe 2 sets, one for conflict and one for shared optima?
	
##since one simulation generates a datafile of about 5MB on disk, then 200 would be 1000 MB, or about 1 GB. seems totally reasonable amount of space.

basevals <- c(100,100, 100,200, 3,3, 3,2,      1,1, 1000,       30, 0.0005, 0.001,0.001,      0.2,    0.6,0.6,   0.1)
#sim.cotrait(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,
													#timesteps,Lambda,mutprb,fiterrP,fiterrM ,prbHorz, pfP, pfM,FLFC,startmats = "n",zoptvects = "n")


parm <- data.frame(matrix(rep(basevals,times=111),nrow=111,byrow=T)) #81 is the base case
parm[1:10,1] <- popsz.v
parm[1:10,2] <- popsz.v
parm[11:20,3] <- nloc.v
parm[11:20,4] <- nloc.v*2
parm[21:30,9] <- wP.v
parm[31:40,10] <- wM.v
parm[41:50,12] <- Lambda.v
parm[51:60,13] <- mutprb.v
parm[61:70,14] <- fiterrP.v
parm[71:80,15] <- fiterrM.v
parm[81:90,16] <- prbHorz.v
parm[91:100,17] <- pfP.v
parm[101:110,18] <- pfP.v
#111th row is the base state
parm2 <- rbind(parm,parm)
parm2[112:222,8] <- 3 #change to fitness agreement; now both have optima at 3
#then this run as array job so that each row of parameters went into sim.cotrait
#simres outputs are saved with names: where jn is 1:nrow(parm2)
#file=paste(Sys.getenv("SCRATCH"),'/sensitivity_',jn,'rep',repnum,'.RData',sep=""))
#repnum goes from 1 to 5 

#they have the following stats output to a file
# FC <- getfitcon(10, pv[11]+1, 1, simres,zoP=pv[7],zoM=pv[8], wP=pv[9], wM=pv[10],pfP=pv[17],pfM=pv[18])$fitnesscorrelation
# VmVp <- extractVmVp(simres, 1,pv[11]+1,1)
# pVp <- VmVp$PVp
# tVp <- VmVp$Vp#currently pVx is a ratio of each to the sum, but not to the breeding value variance.
# pVm <- VmVp$PVm
# tVm <- VmVp$Vm
# win  <- extractwinning(simres,first=1,last=pv[11]+1,1,zoP=pv[7],zoM= pv[8])
# dP <- win$dP
# dM <- win$dM
# dyn <- extractDyn(simres,first=1,last=pv[11]+1,10)
# tcoefP <- dyn$tcoefP
# tcoefM <- dyn$tcoefM
#stats <- data.frame( FC=c(rep(0,times=9),FC), pVp=pVp, tVp=tVp, pVm=pVm, tVm=tVm, dP=dP, dM=dM, 
#			tcoefP=c(0,rep(tcoefP,each=20)), tcoefM= c(0,rep(tcoefM,each=20) ))#a cheater move to make tcoefP the same length....
#named:
#file=paste(Sys.getenv("SCRATCH"),'/sensitivity_stats',jn,'rep',repnum,'.csv',sep=""),row.names=F)

#get last row of each simulation data output and append to matrix

examplefinal <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/sens_stats/sensitivity_stats',1,'rep',1,'.csv',sep=""),header=T)
timecol <- which(formalArgs(sim.cotrait)=="timesteps")

finalstate <-  matrix(NA,nrow=0,ncol=ncol(examplefinal))#ncol in stats file
# finalvar <- matrix(NA,nrow=0,ncol=ncol(examplefinal))#ncol in stats file

# for(i in 1:nrow(parm2)){
# 	simoutreps <- matrix(NA,nrow=0,ncol=ncol(examplefinal))#ncol in stats file
# 	for(r in 1:5){
# 		simout <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/sensitivity_stats',i,'rep','r','.csv',sep=""),header=T)
# 		simoutreps <- rbind(simoutreps, simout[ parm2[i,timecol]+1 ,] )
# 	}
# 		finalstate <- cbind(finalstate,colMeans(simoutreps))#column #needs to be timesteps parameter
# 		finalvar <- cbind(finalvar,sapply(1:ncol(simoutreps), function(column) var(simoutreps[,column]) ))#column #needs to be timesteps parameter
# }


for(i in 1:nrow(parm2)){
	for(r in 1:10){
		simout <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/sens_stats/sensitivity_stats',i,'rep',r,'.csv',sep=""),header=T)
		finalstate <- rbind(finalstate, simout[ parm2[i,timecol]+1 ,] )
	}
}

finalstate$repl <- rep(1:10,times=nrow(parm2))


numargs <- length(formalArgs(sim.cotrait))
colnames(parm2)<- formalArgs(sim.cotrait)[-c(numargs-1,numargs)]

parmWfs <- cbind(parm2[rep(1:nrow(parm2),each=10),],finalstate)
# parmWfs <- cbind(parm2,finalstate)
# parmVfs <- cbind(parm2,finalstate)
# 
# write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_av_parameters.csv",sep=""),row.names=F)
# write.csv(parmVfs,file=paste(Sys.getenv("SCRATCH"),"/sens_var_parameters.csv",sep=""),row.names=F)

write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_reps_finalparameters.csv",sep=""),row.names=F)
