#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to evaluate simulation outputs parameter ranges.
##

Reps <- 5

source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

#how do the following parameters change ans to above:

		popsz.v <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000,10000) #note turning up NM without NP is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
# 		nlocP.v = c(2, 2, 4, 8, 16, 32, 64, 128, 256, 516),#  base set to 20 --
		nloc.v = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)# this has been lowered. since last run
		w.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75. #remains unchanged
		Lambda.v <- seq(from = 35, to = 17, by =-2) #base 25
		mutprb.v <- c(0.0000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.1) #base 0.0001
		prbHorz.v <- seq(from =0, to =1, length.out=10)  #unchanged
		alpha.v <- seq(from = 0.0, to =0.9, by =0.1) #base 0.6	 #unchanged

	
##since one simulation generates a datafile of about 5MB on disk, then 200 would be 1000 MB, or about 1 GB. seems totally reasonable amount of space.

basevals <- c(2000,2000, 20,40, 3,3, 3,2,      0.75,0.75, 1000,       25, 0.0005,      0.2,    0.6,0.6,   0.1)
#sim.cotrait(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,
													#timesteps,Lambda,mutprb,prbHorz, pfP, pfM,FLFC,startmats = "n",zoptvects = "n")


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
# parm[61:70,14] <- fiterrP.v
# parm[71:80,15] <- fiterrM.v
parm[81:90,14] <- prbHorz.v
parm[91:100,15] <- alpha.v #repeated 2x! once for plants, once for microbes
parm[101:110,16] <- alpha.v
#111st and 222nd rows are the base state
parm2 <- rbind(parm,parm)
parm2[112:222,8] <- 3 #change to fitness agreement; now both have optima at 3
#

# repnum from 1-5?, as each gets repeated 5 times.


#they have the following stats output to a file
# FC <- getfitcon(10, pv[11]+1, 1, simres,zoP=pv[7],zoM=pv[8], wP=pv[9], wM=pv[10],pfP=pv[15],pfM=pv[16])$fitnesscorrelation
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
	for(r in 1:Reps){
		simout <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/sens_stats/sensitivity_stats',i,'rep',r,'.csv',sep=""),header=T)
		finalstate <- rbind(finalstate, simout[ parm2[i,timecol]+1 ,] )
	}
}

finalstate$repl <- rep(1:Reps,times=nrow(parm2))


numargs <- length(formalArgs(sim.cotrait))
colnames(parm2)<- formalArgs(sim.cotrait)[-c(numargs-1,numargs)]

parmWfs <- cbind(parm2[rep(1:nrow(parm2),each=Reps),],finalstate)
# parmWfs <- cbind(parm2,finalstate)
# parmVfs <- cbind(parm2,finalstate)
# 
# write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_av_parameters.csv",sep=""),row.names=F)
# write.csv(parmVfs,file=paste(Sys.getenv("SCRATCH"),"/sens_var_parameters.csv",sep=""),row.names=F)

write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_reps_finalparameters.csv",sep=""),row.names=F)
