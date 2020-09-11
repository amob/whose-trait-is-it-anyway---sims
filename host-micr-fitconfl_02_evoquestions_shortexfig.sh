#!/bin/bash
#SBATCH --time=01:40:00
#SBATCH --job-name=simdemofigs2
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca


module load NiaEnv/2018a # allows loading the correct version of R on the computer on which this was run
module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_02_evoquestions_shortexfig.R

#output is always to scratch
