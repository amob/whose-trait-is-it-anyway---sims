#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --job-name=procfitconevo_feedbacks
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca


module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_06_process_specialparam.R

#output is always to scratch
