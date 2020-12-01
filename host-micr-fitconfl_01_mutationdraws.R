#This script produces the mutation distribution figure

#mutate function copied from simulation function
mutate.exp <- function(nL,N,lambda, prbmut) { sapply( 1:N , function(z)
			abs(rexp(nL,rate=lambda)) * ifelse(rbinom(nL,size=1,prob=0.5), 1 , -1 ) * rbinom(nL, size=1, prob = prbmut) ) 
			  }# DistofAbsValueofTraitEff * ProbofPosvNegMutation * ProbMutOccurs(u*loci) 
#high values of lambda give LOWER effect sizes of mutations on average.


vect <- mutate.exp(1,100000,lambda=25,prbmut=1)

sum(-0.002 < vect & vect < 0.002) /length(vect)
sum(-0.02 < vect & vect < 0.02) /length(vect)

pdf("~/Dropbox/host microbe trait evo and gwas/whose-trait-is-it-anyway---sims/mutationdistribution.pdf",height=3,width=3)
par(mar=c(3,3,1.5,1))
plot(density(vect),xlab="",main ="",ylab="",xlim=c(-0.4,0.4))
abline(v=0,lty=3)
mtext("mutation effect",side=1,line=2)
mtext("density",side=2,line=2)
mtext(expression(lambda~"= 25"),side=3,line=0)
dev.off()