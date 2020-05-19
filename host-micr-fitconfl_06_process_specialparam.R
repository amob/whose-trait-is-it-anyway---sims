#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to evaluate simulation outputs parameter ranges.
##


source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

# 
# ###COPIED FROM 05 .R script
		pf.v <- c(1:10)/10 # can have 10 now
		w.v <- c( 0.25,0.5,0.75, 1,1.25,1.5,1.75)#7
		zopt.v <- c(2:5)#4
#full factorial combination for each of these as microbe and plant parameters.	
##since one simulation generates a datafile of about 5MB on disk, then 200 would be 1000 MB, or about 1 GB. seems totally reasonable amount of space.
basevals <- c(100,100, 100,200, 3,3, 5,5,      1,1, 1000,       25, 0.0005,     0.2,    0.6,0.6,   0.1)
params <- data.frame(	pfm =   rep( rep( pf.v, times=4),    times=7 ), #pfp
						zoptM = rep( rep(zopt.v, each=10), times = 7  ), #zoptp
						wM =    rep(w.v, each=40) #wm
			)
parm <- data.frame(matrix(rep(basevals,times=nrow(params)),nrow=nrow(params),byrow=T)) #
parm[,8] <- params$zoptM
parm[,10] <- params$wM
parm[,16] <- params$pfm

# #how do the following parameters change ans to above:
# 		w.v <- c( 0.25, 1)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
# 		pf.v <- c(0,0.4,0.9) #
# 		zopt.v <- c(1,3,5)
# #full factorial combination for each of these as microbe and plant parameters.	
# ##since one simulation generates a datafile of about 5MB on disk, then 200 would be 1000 MB, or about 1 GB. seems totally reasonable amount of space.
# basevals <- c(100,100, 100,200, 3,3, 3,2,      1,1, 1000,       25, 0.0005,     0.2,    0.6,0.6,   0.1)
# #sim.cotrait(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,
# 													#timesteps,Lambda,mutprb,prbHorz, pfP, pfM,FLFC,startmats = "n",zoptvects = "n")
# 
# params <- data.frame(	pfp = rep( rep( rep( rep(pf.v, times=3),    times=3 ), times=3), times=4), #pfp
# 						pfm = rep( rep( rep( rep(pf.v, each =3),    times=3 ), times=3), times=4), #pfm
# 						zoptP = rep( rep( rep( rep(zopt.v, each=3), each = 3  ), times=3), times=4), #zoptp
# 						zoptM = rep( rep( rep( rep(zopt.v, each =3), each = 3 ), each =3), times=4), #zoptm			
# 						wP = rep( rep(w.v, each=81), times=2 ), #wp
# 						wM = rep(w.v, each=162) #wm
# 			)
# parm <- data.frame(matrix(rep(basevals,times=nrow(params)),nrow=nrow(params),byrow=T)) #
# parm[,7] <- params$zoptP
# parm[,8] <- params$zoptM
# parm[,9] <- params$wP
# parm[,10] <- params$wM
# parm[,15] <- params$pfp
# parm[,16] <- params$pfm
# 

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
#file=paste(Sys.getenv("SCRATCH"),'/feedback_stats/sensitivity_stats',jn,'rep',repnum,'.csv',sep=""),row.names=F)

#get last row of each simulation data output and append to matrix

examplefinal <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/feedback_stats/sensitivityholdplant_stats',1,'rep',1,'.csv',sep=""),header=T)
timecol <- which(formalArgs(sim.cotrait)=="timesteps")

finalstate <-  matrix(NA,nrow=0,ncol=ncol(examplefinal))#ncol in stats file

###ASSUMES 10 replicates run from here on out!

for(i in 1:nrow(parm)){
	for(r in 1:10){
		simout <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/feedback_stats/sensitivityholdplant_stats',i,'rep',r,'.csv',sep=""),header=T)
		finalstate <- rbind(finalstate, simout[ parm[i,timecol]+1 ,] )
	}
}

finalstate$repl <- rep(1:10,times=nrow(parm))


numargs <- length(formalArgs(sim.cotrait))
colnames(parm)<- formalArgs(sim.cotrait)[-c(numargs-1,numargs)]

parmWfs <- cbind(parm[rep(1:nrow(parm),each=10),],finalstate)

write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_reps_feedbackparameters_holdplant.csv",sep=""),row.names=F)
