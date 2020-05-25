#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=xpphenos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca


#temporarily comment out next two lines if no need to rerun
module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl-09-ExperimentalPhenos.R 
