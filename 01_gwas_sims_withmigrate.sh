#!/bin/bash
#SBATCH --time=40:00
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
plink --file $SCRATCH/HOLOevosims --out HOLOevosims --make-bed
plink --file $SCRATCH/PLANTevosims --out PLANTevosims --make-bed
plink --file $SCRATCH/MICRevosims --out MICRevosims --make-bed
#note that if split up from R script, the input file location may need to be changed

##PASS to GEMMA to run models.

#get kinship
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosims  -gk 1 -o holokin
#given manual text, might actually prefer -gk 2, since the larger effect alleles will be low frequency.
#because of different optima across pops...etc
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosims  -gk 1 -o plantkin
$HOME/gemma-0.98.1-linux-static -bfile MICRevosims  -gk 1 -o micrkin

#run models
#without K
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosims -lm 4 -o HOLOgemma
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosims -lm 4 -o PLANTgemma
$HOME/gemma-0.98.1-linux-static -bfile MICRevosims -lm 4 -o MICRgemma

#mixed models
$HOME/gemma-0.98.1-linux-static -bfile HOLOevosims -k $SCRATCH/output/holokin.cXX.txt -lmm 4 -o HOLOgemmaK
$HOME/gemma-0.98.1-linux-static -bfile MICRevosims -k $SCRATCH/ouput/micrkin.cXX.txt -lmm 4 -o MICRgemmaK
$HOME/gemma-0.98.1-linux-static -bfile PLANTevosims -k $SCRATCH/ouput/plantkin.cXX.txt -lmm 4 -o PLANTgemmaK
#mixed models currently fail for plant and microbes but not holo, NO IDEA WHY. otherwise everything runs!
