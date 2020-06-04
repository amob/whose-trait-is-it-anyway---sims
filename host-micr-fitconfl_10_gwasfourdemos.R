#######
#Quantitative trait evolution and fitness conflict or not in plant-microbe interactions
#######

##
#goal of this script is
#for 4 scenarios, run in silico simulations and prep for gwas
##

source(paste(Sys.getenv("HOME"),'/whosetrait/host-micr-fitconfl_10_gwasfunctionsonly.R',sep="")) 
#source('~/Dropbox/host microbe trait evo and gwas/host-micr-fitconfl_10_gwasfunctionsonly.R') 

load(file=paste(Sys.getenv("SCRATCH"),'/Simulation_Results_fourdemos.RData',sep=""))

test.4b <- fourdemosims[[1]]
test.4bff <- fourdemosims[[2]]
test.4f <- fourdemosims[[3]]
test.4fff <- fourdemosims[[4]]

#.005 seems very low for experimental error, but there for now so that the range of effects is less than larger allele effect sizes
expsetabo4b <- 		run.exp.allbyone(test.4b, numperpop= 100,numsites=1, exp.err=0.005,nreps=4) 
expsetabo4bff <- 	run.exp.allbyone(test.4bff,numperpop= 100,numsites=1, exp.err=0.005,nreps=4) 
expsetabo4f <- 		run.exp.allbyone(test.4f,numperpop= 100,numsites=1, exp.err=0.005,nreps=4) 
expsetabo4fff <- 	run.exp.allbyone(test.4fff, numperpop= 100,numsites=1, exp.err=0.005,nreps=4) 

fourdemoexps <- list(expsetabo4b,expsetabo4bff,expsetabo4f,expsetabo4fff)
save(fourdemoexps,file=paste(Sys.getenv("SCRATCH"),'/fourdemoexps_abo.RData',sep=""))

# plinkaboH <- makegwasfiles(expsetabo,expsetabo$expdat,type="HOLO") 

plinkaboP4b <-   makegwasfiles(expsetabo4b,  expdat=expsetabo4b$expdatPexp,  type="PLANT") 
	write.table(plinkaboP4b$geno,file=paste(Sys.getenv("SCRATCH"),"/PLANT4b_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboP4b$map, file=paste(Sys.getenv("SCRATCH"),"/PLANT4b_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboP4bff <- makegwasfiles(expsetabo4bff,expdat=expsetabo4bff$expdatPexp,type="PLANT") 
	write.table(plinkaboP4bff$geno,file=paste(Sys.getenv("SCRATCH"),"/PLANT4bff_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboP4bff$map, file=paste(Sys.getenv("SCRATCH"),"/PLANT4bff_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboP4f <-   makegwasfiles(expsetabo4f,  expdat=expsetabo4f$expdatPexp,  type="PLANT") 
	write.table(plinkaboP4f$geno,file=paste(Sys.getenv("SCRATCH"),"/PLANT4f_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboP4f$map, file=paste(Sys.getenv("SCRATCH"),"/PLANT4f_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboP4fff <- makegwasfiles(expsetabo4fff,expdat=expsetabo4fff$expdatPexp,type="PLANT") 
	write.table(plinkaboP4fff$geno,file=paste(Sys.getenv("SCRATCH"),"/PLANT4fff_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboP4fff$map, file=paste(Sys.getenv("SCRATCH"),"/PLANT4fff_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

plinkaboM4b <- makegwasfiles(expsetabo4b,expdat=expsetabo4b$expdatMexp,type="MICR") 
	write.table(plinkaboM4b$geno, file=paste(Sys.getenv("SCRATCH"), "/MICR4b_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboM4b$map,  file=paste(Sys.getenv("SCRATCH"), "/MICR4b_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboM4bff <- makegwasfiles(expsetabo4b,expdat=expsetabo4bff$expdatMexp,type="MICR")
	write.table(plinkaboM4bff$geno, file=paste(Sys.getenv("SCRATCH"), "/MICR4bff_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboM4bff$map,  file=paste(Sys.getenv("SCRATCH"), "/MICR4bff_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboM4f <- makegwasfiles(expsetabo4f,expdat=expsetabo4f$expdatMexp,type="MICR") 
	write.table(plinkaboM4f$geno, file=paste(Sys.getenv("SCRATCH"), "/MICR4f_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboM4f$map,  file=paste(Sys.getenv("SCRATCH"), "/MICR4f_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
plinkaboM4fff <- makegwasfiles(expsetabo4fff,expdat=expsetabo4fff$expdatMexp,type="MICR")
	write.table(plinkaboM4fff$geno, file=paste(Sys.getenv("SCRATCH"), "/MICR4fff_ABO.ped",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
	write.table(plinkaboM4fff$map,  file=paste(Sys.getenv("SCRATCH"), "/MICR4fff_ABO.map",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

 #store allele affect size information for causal genos
plantlociABO4b <- expsetabo4b$causalgenos$locusdatP #alleles 1:
	plantlociABO4b$locusname <- paste("cP",1:nrow(plantlociABO4b),sep="") #names as in gwas
	write.csv(plantlociABO4b,file=paste(Sys.getenv("SCRATCH"),'/plantlociABO4b.csv',sep=""),row.names=F)
plantlociABO4bff <- expsetabo4bff$causalgenos$locusdatP 
	plantlociABO4bff$locusname <- paste("cP",1:nrow(plantlociABO4bff),sep="") 
	write.csv(plantlociABO4bff,file=paste(Sys.getenv("SCRATCH"),'/plantlociABO4bff.csv',sep=""),row.names=F)
plantlociABO4f <- expsetabo4f$causalgenos$locusdatP 
	plantlociABO4f$locusname <- paste("cP",1:nrow(plantlociABO4f),sep="") 
	write.csv(plantlociABO4f,file=paste(Sys.getenv("SCRATCH"),'/plantlociABO4f.csv',sep=""),row.names=F)
plantlociABO4fff <- expsetabo4fff$causalgenos$locusdatP 
	plantlociABO4fff$locusname <- paste("cP",1:nrow(plantlociABO4fff),sep="") 
	write.csv(plantlociABO4fff,file=paste(Sys.getenv("SCRATCH"),'/plantlociABO4fff.csv',sep=""),row.names=F)


micrlociABO4b <- expsetabo4b$causalgenos$locusdatM #alleles 1:
	micrlociABO4b$locusname <- paste("cM",1:nrow(micrlociABO4b),sep="") #names as in gwas
	write.csv(micrlociABO4b,file=paste(Sys.getenv("SCRATCH"),'/micrlociABO4b.csv',sep=""),row.names=F)
micrlociABO4bff <- expsetabo4bff$causalgenos$locusdatM #alleles 1:
	micrlociABO4bff$locusname <- paste("cM",1:nrow(micrlociABO4bff),sep="") #names as in gwas
	write.csv(micrlociABO4bff,file=paste(Sys.getenv("SCRATCH"),'/micrlociABO4bff.csv',sep=""),row.names=F)
micrlociABO4f <- expsetabo4f$causalgenos$locusdatM #alleles 1:
	micrlociABO4f$locusname <- paste("cM",1:nrow(micrlociABO4f),sep="") #names as in gwas
	write.csv(micrlociABO4f,file=paste(Sys.getenv("SCRATCH"),'/micrlociABO4f.csv',sep=""),row.names=F)
micrlociABO4fff <- expsetabo4fff$causalgenos$locusdatM #alleles 1:
	micrlociABO4fff$locusname <- paste("cM",1:nrow(micrlociABO4fff),sep="") #names as in gwas
	write.csv(micrlociABO4fff,file=paste(Sys.getenv("SCRATCH"),'/micrlociABO4fff.csv',sep=""),row.names=F)
