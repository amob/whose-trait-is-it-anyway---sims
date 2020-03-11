#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is to run simulations across different parameter ranges.
#might want to make parameter ranges shell/system variables in some way
##

source('/home/f/freder19/annamob/whosetrait/host-micr-fitconfl_01_simfunction.R') 

#how do the following parameters change ans to above:
	#strength of link between trait and fitness; wP and wM, lower is stronger link -- with lower links, PFF becomes more important
		wP.v <- c(0.1, 0.15, 0.25, 0.5, 0.75, 1.25, 1.5, 2, 2.5, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
		wM.v <- c(0.1, 0.15, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.5 ,5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75?
	#drift and selection within each species; # uhhh. wP and wM vs fiterrP and fiterrM
		#N and s
		fiterrP.v <- c(0.0001, 0.0005, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		fiterrM.v <- c(0.0001, 0.0005, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		popsz.v <- c(10, 25, 50, 75, 125, 150, 200, 250, 300, 500) #note turning up M without N is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
		nloc.v <- c(10, 25, 50, 75, 125, 150, 200, 250, 300, 500) # multiply by 2 for microbes, base set to 100
		prbHorz.v <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.25, 0.5)
	# mutational inputs;  Lambda
		Lambda.v <- seq(from = 41, to = 23, by =-2) #base 30
		mutprb.v <- c(0.00001, 0.000025, 0.00005, 0.000075, 0.0001, 0.00025, 0.00075, 0.001 ,0.0015, 0.002)
	#symmetries in conflict between partners (are these poss? I think just change pff for one only)
		pfP.v <- seq(from = 0.0, to =0.95, by =0.1) #base 0.6	
	#fitness optima?
		#maybe 2 sets, one for conflict and one for shared optima?
	
###needs code to run sims :)
##since one simulation generates a datafile of about 5MB on disk, then 200 would be 1000 MB, or about 1 GB. seems totally reasonable amount of space.

basevals <- c(100,100, 100,200, 3,3, 3,2,      1,1, 500,       30, 0.0005, 0.001,0.001,      0.2,    0.6,0.6,   0.1)
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

parm2 <- cbind(parm,parm)
parm[112:222,8] <- 3 #fitness agreement on 3

#pull out a row corresponding to a system variable that is the number in the arrayjob
jobnum <- 
pv <- parm[jobnum,]

#run sim cotrait on those variables
simres <- sim.cotrait(NP=pv[1],NM=pv[2],nlP=pv[3],nlM=pv[4],nlnP=pv[5],nlnM=pv[6],
			zoP=pv[7],zoM=pv[8],wP=pv[9],wM=pv[10],
			timesteps = pv[11], Lambda = pv[12],mutprb =pv[13],fiterrP=pv[14],fiterrM =pv[15],
			prbHorz = pv[16], pfP =pv[17], pfM=pv[18])

#write sim with jobbum in the name

#calculate diagnostic stats?

#write stats with jobnum in name
