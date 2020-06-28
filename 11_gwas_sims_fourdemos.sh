#!/bin/bash
#SBATCH --time=0:12:00
#SBATCH --job-name=gwassim4demo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca
#SBATCH --partition=debug


#temporarily comment out next two lines if no need to rerun
module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_10_gwasfourdemos.R 
#outputs .ped and .map for plant and microbe and each of the 4 scenarios
#according to simulation parameters and experiment parameters as set up in script
#output is always to scratch

module load intel/2018.2
module load plink

plink --file $SCRATCH/PLANT4b_ABO --out PLANT4bABO --make-bed
plink --file $SCRATCH/PLANT4bff_ABO --out PLANT4bffABO --make-bed
plink --file $SCRATCH/PLANT4a_ABO --out PLANT4aABO --make-bed
plink --file $SCRATCH/PLANT4aff_ABO --out PLANT4affABO --make-bed
plink --file $SCRATCH/PLANT4f_ABO --out PLANT4fABO --make-bed
plink --file $SCRATCH/PLANT4fff_ABO --out PLANT4fffABO --make-bed

plink --file $SCRATCH/MICR4b_ABO --out MICR4bABO --make-bed
plink --file $SCRATCH/MICR4bff_ABO --out MICR4bffABO --make-bed
plink --file $SCRATCH/MICR4a_ABO --out MICR4aABO --make-bed
plink --file $SCRATCH/MICR4aff_ABO --out MICR4affABO --make-bed
plink --file $SCRATCH/MICR4f_ABO --out MICR4fABO --make-bed
plink --file $SCRATCH/MICR4fff_ABO --out MICR4fffABO --make-bed

##PASS to GEMMA to run models.
#get kinship; currently not using
# $HOME/gemma-0.98.1-linux-static -bfile PLANTevosimsABO  -gk 2 -o plantkinABO
# $HOME/gemma-0.98.1-linux-static -bfile MICRevosimsABO  -gk 2 -o micrkinABO
#given manual text, might actually prefer -gk 2, since the larger effect alleles will be low frequency (depends somewhat on simulation parameters)
#however, they say -gk 1 usually performs better.

#run models without K
$HOME/gemma-0.98.1-linux-static -bfile PLANT4bABO -lm 4 -o PLANTgemmaABO4b -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile PLANT4bffABO -lm 4 -o PLANTgemmaABO4bff -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile PLANT4aABO -lm 4 -o PLANTgemmaABO4a -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile PLANT4affABO -lm 4 -o PLANTgemmaABO4aff -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile PLANT4fABO -lm 4 -o PLANTgemmaABO4f -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile PLANT4fffABO -lm 4 -o PLANTgemmaABO4fff -maf 0.00375

$HOME/gemma-0.98.1-linux-static -bfile MICR4bABO -lm 4 -o MICRgemmaABO4b -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile MICR4bffABO -lm 4 -o MICRgemmaABO4bff -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile MICR4aABO -lm 4 -o MICRgemmaABO4a -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile MICR4affABO -lm 4 -o MICRgemmaABO4aff -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile MICR4fABO -lm 4 -o MICRgemmaABO4f -maf 0.00375
$HOME/gemma-0.98.1-linux-static -bfile MICR4fffABO -lm 4 -o MICRgemmaABO4fff -maf 0.00375

