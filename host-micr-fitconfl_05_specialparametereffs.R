#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to run simulations across different parameter ranges.
##
options(scipen = 999)

repnum <- Sys.getenv("REP")

source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

#how do the following parameters change ans to above:
		pf.v <- c(1:10)/10 # 10 levels
		w.v <- c( 0.25,0.5,0.75, 1,1.25,1.5,1.75)#7
		zopt.v <- c(2:5)#4
#full factorial combination for each of these as microbe and plant parameters.	
##since one simulation generates a datafile of about 5MB on disk, then 200 would be 1000 MB, or about 1 GB. seems totally reasonable amount of space.

basevals <- c(2000,2000, 20,40, 3,3, 2,2,      0.75,0.75, 1000,       25, 0.0001,      0.2,    0.6,0.6)
#sim.cotrait(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,
													#timesteps,Lambda,mutprb,prbHorz, pfP, pfM,startmats = "n",zoptvects = "n")


params <- data.frame(	pfm =   rep( rep( pf.v, times=4),    times=7 ), #pfp
						zoptP = rep( rep(zopt.v, each=10), times = 7  ), #zoptp
						wM =    rep(w.v, each=40) #wm
			)


#check
nrow(params)
length(unique(sapply(1:nrow(params), function(z) paste(params[z,],collapse=".") )))
#all same length!

parm <- data.frame(matrix(rep(basevals,times=nrow(params)),nrow=nrow(params),byrow=T)) #
parm[,7] <- params$zoptP
# parm[,8] <- params$zoptM # other variables that one might vary are commented out, not varied in this script
#parm[,9] <- params$wP
parm[,10] <- params$wM
#parm[,15] <- params$pfp
parm[,16] <- params$pfm


#####THIS SECTION ASSUMES RUNNING AS ARRAY JOB FROM PAIRED BASH SCRIPT
#pull out a row corresponding to a system variable that is the number in the arrayjob
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# numeric job number
jn <- as.numeric(slurm_arrayid)

pv <- as.numeric(parm[jn,])

jn

pv

#run sim cotrait on those variables
simres <- sim.cotrait(NP=pv[1],NM=pv[2],nlP=pv[3],nlM=pv[4],nlnP=pv[5],nlnM=pv[6],
#sim.cotrait(         NP,NM,            nlP,nlM,            nlnP,nlnM,
			zoP=pv[7],zoM=pv[8],wP=pv[9],wM=pv[10],
#			zoP,zoM,           wP,wM,
			timesteps = pv[11], Lambda = pv[12],mutprb =pv[13],
			#timesteps,Lambda,mutprb,
			prbHorz = pv[14], pfP =pv[15], pfM=pv[16])
		# ,prbHorz, pfP, pfM,startmats = "n",zoptvects = "n")


#write sim with jobbum in the name
save(simres,file=paste(Sys.getenv("SCRATCH"),'/feedback_rdata/sensitivityholdplant_',jn,'rep',repnum,'.RData',sep=""))

#calculate diagnostic stats
FC <- getfitcon(10, pv[11]+1, 1, simres,zoP=pv[7],zoM=pv[8], wP=pv[9], wM=pv[10],pfP=pv[15],pfM=pv[16])$fitnesscorrelation
VmVp <- extractVmVp(simres, 1,pv[11]+1,1)
pVp <- VmVp$PVp
tVp <- VmVp$Vp#
pVm <- VmVp$PVm
tVm <- VmVp$Vm
win  <- extractwinning(simres,first=1,last=pv[11]+1,1,zoP=pv[7],zoM= pv[8])
dP <- win$dP
dM <- win$dM
dyn <- extractDyn(simres,first=1,last=pv[11]+1,20)
tcoefP <- dyn$tcoefP
tcoefM <- dyn$tcoefM

stats <- data.frame( FC=c(rep(0,times=9),FC), pVp=pVp, tVp=tVp, pVm=pVm, tVm=tVm, dP=dP, dM=dM, 
			tcoefP=c(0,rep(tcoefP,each=20)), tcoefM=c(0,rep(tcoefM,each=20) ))#the last argument makes tcoefP/M the same length as the other variables through repeating the value
#write stats with jobnum in name

write.csv(stats,file=paste(Sys.getenv("SCRATCH"),'/feedback_stats/sensitivityholdplant_stats',jn,'rep',repnum,'.csv',sep=""),row.names=F)