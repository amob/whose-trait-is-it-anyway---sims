
##
#goal of this script is to evaluate simulation outputs parameter ranges.
##

Reps <- 5 #this number should reflect the number of times the corresponding 02 series shell script was run -- with updated replicate numbers from 1:Reps

source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 

# ###THIS SECTION COPIED FROM 02 series .R script, should be updated to remain identical to 02 series if any changes are made
		popsz.v <- c(100, 500, 900, 1300, 1700, 2100, 2500, 2900, 3300, 4100) 
		nloc.v = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
		w.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)
		Lambda.v <- seq(from = 35, to = 17, by =-2) 
		mutprb.v <- c( 0.0001, 0.0002, 0.0003,0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001 )
		prbHorz.v <- seq(from =0, to =1, length.out=10)  
		alpha.v <- seq(from = 0.0, to =0.9, by =0.1) 

basevals <- c(2000,2000, 20,40, 3,3, 3,2,      0.75,0.75, 1000,       25, 0.0005,      0.2,    0.6,0.6)
#sim.cotrait(NP,NM,nlP,nlM,nlnP,nlnM,zoP,zoM,wP,wM,
													#timesteps,Lambda,mutprb,prbHorz, pfP, pfM,startmats = "n",zoptvects = "n")
parm <- data.frame(matrix(rep(basevals,times=111),nrow=111,byrow=T)) #81 is the base case
parm[1:10,1] <- popsz.v
parm[1:10,2] <- popsz.v
parm[11:20,3] <- nloc.v #L_P
parm[21:30,4] <- nloc.v #L_M
parm[31:40,3] <- nloc.v #L_P when varying loci for M and P tog
parm[31:40,4] <- nloc.v*2 ##L_P when varying loci for M and P tog
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

#get last row of each simulation data output and append to matrix
examplefinal <- read.csv(file=paste(Sys.getenv("SCRATCH"),'/sens_stats/sensitivity_stats',1,'rep',1,'.csv',sep=""),header=T)
timecol <- which(formalArgs(sim.cotrait)=="timesteps")

finalstate <-  matrix(NA,nrow=0,ncol=ncol(examplefinal))#ncol in stats file

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

write.csv(parmWfs,file=paste(Sys.getenv("SCRATCH"),"/sens_reps_finalparameters.csv",sep=""),row.names=F)
