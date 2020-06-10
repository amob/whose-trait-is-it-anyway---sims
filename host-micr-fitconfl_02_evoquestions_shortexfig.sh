#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --job-name=simdemofigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca


module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_02_evoquestions_shortexfig.R

#output is always to scratch
