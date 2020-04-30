#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --job-name=simandgwastest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca


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
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABA  -gk 1 -o holokinABA
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABO  -gk 1 -o holokinABO
#given manual text, might actually prefer -gk 2, since the larger effect alleles will be low frequency.
#because of different optima across pops...etc
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABA  -gk 1 -o plantkinABA
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABA  -gk 1 -o micrkinABA
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO  -gk 1 -o plantkinABO
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO  -gk 1 -o micrkinABO

#run models
#without K
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABA -lm 4 -o HOLOgemmaABA
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABA -lm 4 -o PLANTgemmaABA
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABA -lm 4 -o MICRgemmaABA
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABO -lm 4 -o HOLOgemmaABA
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO -lm 4 -o PLANTgemmaABA
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO -lm 4 -o MICRgemmaABA

#mixed models


$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABA -k $SCRATCH/output/holokinABA.cXX.txt -lmm 4 -o HOLOgemmaKABA
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABA -k $SCRATCH/ouput/micrkinABA.cXX.txt -lmm 4 -o MICRgemmaKABA
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABA -k $SCRATCH/ouput/plantkinABA.cXX.txt -lmm 4 -o PLANTgemmaKABA
#
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosimsABO -k $SCRATCH/output/holokinABO.cXX.txt -lmm 4 -o HOLOgemmaKABO
$HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO -k $SCRATCH/ouput/micrkinABO.cXX.txt -lmm 4 -o MICRgemmaKABO
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO -k $SCRATCH/ouput/plantkinABO.cXX.txt -lmm 4 -o PLANTgemmaKABO
#mixed models currently fail for plant and microbes 

