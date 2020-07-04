#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to evaluate simulation outputs parameter ranges.
##


source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 
reps <- 5 #SET REPLICATES

# 
# ###COPIED FROM 05 series .R script
		pf.v <- c(1:10)/10 # can have 10 now
		w.v <- c( 0.25,0.5,0.75, 1,1.25,1.5,1.75)#7
		zopt.v <- c(2:5)#4
#full factorial combination for each of these as microbe and plant parameters.	
##note one simulation generates a datafile of at least 5MB on disk, then 200 would be 1000 MB, or about 1 GB. 
basevals <- c(2000,2000, 20,40, 3,3, 2,2,      0.75,0.75, 1000,       25, 0.0001,      0.2,    0.6,0.6,   0.1)
params <- data.frame(	pfm =   rep( rep( pf.v, times=4),    times=7 ), #pfp
						zoptP = rep( rep(zopt.v, each=10), times = 7  ), #zoptp
						wM =    rep(w.v, each=40) #wm
			)
parm <- data.frame(matrix(rep(basevals,times=nrow(params)),nrow=nrow(params),byrow=T)) #
parm[,7] <- params$zoptP
parm[,10] <- params$wM
parm[,16] <- params$pfm


#they have the following stats output to a file like this (copied over from 05 series R script)
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


for(i in 1:nrow(parm)){
	for(r in 1:reps){
		simout <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/feedback_stats/sensitivityholdplant_stats',i,'rep',r,'.csv',sep=""),header=T)
		finalstate <- rbind(finalstate, simout[ parm[i,timecol]+1 ,] )
	}
}

finalstate$repl <- rep(1:reps,times=nrow(parm))


numargs <- length(formalArgs(sim.cotrait))
colnames(parm)<- formalArgs(sim.cotrait)[-c(numargs-1,numargs)]

parmWfs <- cbind(parm[rep(1:nrow(parm),each=reps),],finalstate)

write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_reps_feedbackparameters_holdplant.csv",sep=""),row.names=F)
