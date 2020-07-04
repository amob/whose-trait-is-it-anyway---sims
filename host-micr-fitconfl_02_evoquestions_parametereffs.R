#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to run simulations across different parameter ranges.
##
options(scipen = 999)

repnum <- Sys.getenv("REP")

source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

		popsz.v <- c(100, 500, 900, 1300, 1700, 2100, 2500, 2900, 3300, 4100) #
		nloc.v = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)# 
		w.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)#
		Lambda.v <- seq(from = 35, to = 17, by =-2) #base 25
		mutprb.v <- c(0.00002, 0.00004, 0.00006, 0.00008, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001 )
# 		mutprb.v <- c(0.0000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.1) #base 0.0001
		prbHorz.v <- seq(from =0, to =1, length.out=10)  #
		alpha.v <- seq(from = 0.0, to =0.9, by =0.1) #base 0.6	


	

basevals <- c(2000,2000, 20,40, 3,3, 3,2,      0.75,0.75, 1000,       25, 0.0001,      0.2,    0.6,0.6,   0.1)
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
parm[81:90,14] <- prbHorz.v
parm[91:100,15] <- alpha.v #repeated 2x! once for plants, once for microbes
parm[101:110,16] <- alpha.v
#111st and 222nd rows are the base state
parm2 <- rbind(parm,parm)
parm2[112:222,8] <- 3 #change to fitness agreement; now both have optima at 3
#

#####THIS SECTION ASSUMES RUNNING AS ARRAY JOB FROM PAIRED BASH SCRIPT
#pull out a row corresponding to a system variable that is the number in the arrayjob
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# numeric job number
jn <- as.numeric(slurm_arrayid)

pv <- as.numeric(parm2[jn,])

jn

pv

#run sim cotrait on those variables
simres <- sim.cotrait(NP=pv[1],NM=pv[2],nlP=pv[3],nlM=pv[4],nlnP=pv[5],nlnM=pv[6],
#sim.cotrait(         NP,NM,            nlP,nlM,            nlnP,nlnM,
			zoP=pv[7],zoM=pv[8],wP=pv[9],wM=pv[10],
#			zoP,zoM,           wP,wM,
			timesteps = pv[11], Lambda = pv[12],mutprb =pv[13],
			#timesteps,Lambda,mutprb,
			prbHorz = pv[14], pfP =pv[15], pfM=pv[16],FLFC=0.1)
		# ,prbHorz, pfP, pfM,FLFC,startmats = "n",zoptvects = "n")


#write sim with jobbum in the name
save(simres,file=paste(Sys.getenv("SCRATCH"),'/sens_rdata/sensitivity_',jn,'rep',repnum,'.RData',sep=""))

#calculate diagnostic stats
FC <- getfitcon(10, pv[11]+1, 1, simres,zoP=pv[7],zoM=pv[8], wP=pv[9], wM=pv[10],pfP=pv[15],pfM=pv[16])$fitnesscorrelation
VmVp <- extractVmVp(simres, 1,pv[11]+1,1)
pVp <- VmVp$PVp
tVp <- VmVp$Vp#currently pVx is a ratio of each to the sum, but not to the breeding value variance.
pVm <- VmVp$PVm
tVm <- VmVp$Vm
win  <- extractwinning(simres,first=1,last=pv[11]+1,1,zoP=pv[7],zoM= pv[8])
dP <- win$dP
dM <- win$dM
dyn <- extractDyn(simres,first=1,last=pv[11]+1,20)
tcoefP <- dyn$tcoefP
tcoefM <- dyn$tcoefM

stats <- data.frame( FC=c(rep(0,times=9),FC), pVp=pVp, tVp=tVp, pVm=pVm, tVm=tVm, dP=dP, dM=dM, 
			tcoefP=c(0,rep(tcoefP,each=20)), tcoefM=c(0,rep(tcoefM,each=20) ))#replicate to make tcoefP the same length....as is a window based stat.
#write stats with jobnum in name

write.csv(stats,file=paste(Sys.getenv("SCRATCH"),'/sens_stats/sensitivity_stats',jn,'rep',repnum,'.csv',sep=""),row.names=F)
