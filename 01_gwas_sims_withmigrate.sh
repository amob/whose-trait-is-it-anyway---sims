#!/bin/bash
#SBATCH --time=13:00:00
#SBATCH --job-name=simandgwastest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca


#temporarily comment out next two lines if no need to rerun
module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_01_gwas_sims_withmigrate.R
#outputs .ped and .map for holo, plant, and microbe genomes
#according to simulation parameters and experiment parameters as set up in script
#output is always to scratch

module load intel/2018.2
module load plink

##PASS .ped and .map to PLINK to get binary formats .bim .fam .bed
#plink --file mydata --out mydata --make-bed 
plink --file $SCRATCH/HOLOevosims_ABA --out HOLOevosimsABA --make-bed
plink --file $SCRATCH/PLANTevosims_ABA --out PLANTevosimsABA --make-bed
plink --file $SCRATCH/MICRevosims_ABA --out MICRevosimsABA --make-bed
plink --file $SCRATCH/HOLOevosims_ABO --out HOLOevosimsABO --make-bed
plink --file $SCRATCH/PLANTevosims_ABO --out PLANTevosimsABO --make-bed
plink --file $SCRATCH/MICRevosims_ABO --out MICRevosimsABO --make-bed
#note that if split up from R script, the input file location may need to be changed

##PASS to GEMMA to run models.

#get kinship
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABA  -gk 2 -o holokinABA
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABA  -gk 2 -o plantkinABA
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABA  -gk 2 -o micrkinABA
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABO  -gk 2 -o holokinABO
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO  -gk 2 -o plantkinABO
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO  -gk 2 -o micrkinABO
#given manual text, might actually prefer -gk 2, since the larger effect alleles will be low frequency (depends somewhat on simulation parameters)
#however, they say -gk 1 usually performs better.



#run models
#without K
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABA -lm 4 -o HOLOgemmaABA -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABA -lm 4 -o PLANTgemmaABA -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABA -lm 4 -o MICRgemmaABA -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABO -lm 4 -o HOLOgemmaABO -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO -lm 4 -o PLANTgemmaABO -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO -lm 4 -o MICRgemmaABO -maf 0.0125

#when ABA experiment is 6400 big as 4 reps of 1600 unique p-m combos, 
	#then there are 40 of each p and m in a sq design and the min af possible is 0.0125
# the identical maf in the abo experiment is when there are 20 copies of the 800 (1600 chrm.) plant (or 10 micr in other half) have an allele

#mixed models


$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABA -k $SCRATCH/output/holokinABA.cXX.txt -lmm 4 -o HOLOgemmaKABA  -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABA -k $SCRATCH/ouput/micrkinABA.cXX.txt -lmm 4 -o MICRgemmaKABA -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABA -k $SCRATCH/ouput/plantkinABA.cXX.txt -lmm 4 -o PLANTgemmaKABA -maf 0.0125
#
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABO -k $SCRATCH/output/holokinABO.cXX.txt -lmm 4 -o HOLOgemmaKABO -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO -k $SCRATCH/ouput/micrkinABO.cXX.txt -lmm 4 -o MICRgemmaKABO -maf 0.0125
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO -k $SCRATCH/ouput/plantkinABO.cXX.txt -lmm 4 -o PLANTgemmaKABO -maf 0.0125
#mixed models currently fail for plant and microbes 

