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
		wP.v <- c(0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 1?  
		wM.v <- c(0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 5)#seq(from = 0.25, to = 5,lenght.out=10) # set base at 0.75?
	#drift and selection within each species; # uhhh. wP and wM vs fiterrP and fiterrM
		#N and s
		fiterrP.v <- c(0.0001, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		fiterrM.v <- c(0.0001, 0.001, 0.0015, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.05)#base set both low, to 0.001
		popsz.v <- c(10, 25, 50, 75, 100, 150, 200, 500) #note turning up M without N is similar to increasing fiterr. increasing hosts without microbes makes no sense and is not possible.
		nloc.v <- c(10, 25, 50, 75, 100, 150, 200,500) # multiply by 2 for microbes, base set to 100
	# mutational inputs;  Lambda
		Lambda.v <- seq(from = 40, to = 20, by =-2) #base 30
	#symmetries in conflict between partners (are these poss? I think just change pff for one only)
		pfP.v <- seq(from = 0.05, to =0.95, by =0.1) #base 0.6	


###needs code to run sims :)
