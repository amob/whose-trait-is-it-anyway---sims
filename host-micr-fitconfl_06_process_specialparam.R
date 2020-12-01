##
#goal of this script is to evaluate simulation outputs parameter ranges.
##


source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_01_simfunction.R',sep="")) 
reps <- 5 #SET REPLICATES 
 #this number should reflect the number of times the corresponding 02 series shell script was run -- with updated replicate numbers from 1:Reps

# 
# ###THIS SECTION COPIED FROM 05 series .R script, should be updated to remain identical to 05 series if any changes are made
		pf.v <- c(1:10)/10 
		w.v <- c( 0.25,0.5,0.75, 1,1.25,1.5,1.75)#
		zopt.v <- c(2:5)#
basevals <- c(2000,2000, 20,40, 3,3, 2,2,      0.75,0.75, 1000,       25, 0.0001,      0.2,    0.6,0.6,   0.1)
params <- data.frame(	pfm =   rep( rep( pf.v, times=4),    times=7 ), #pfp
						zoptP = rep( rep(zopt.v, each=10), times = 7  ), #zoptp
						wM =    rep(w.v, each=40) #wm
			)
parm <- data.frame(matrix(rep(basevals,times=nrow(params)),nrow=nrow(params),byrow=T)) #
parm[,7] <- params$zoptP
parm[,10] <- params$wM
parm[,16] <- params$pfm






#get last row of each simulation data output and append to matrix, and output
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
