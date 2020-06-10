#!/bin/bash
#SBATCH --time=0:60:00
#SBATCH --job-name=R1fitconevo_sensitivity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca
#SBATCH --array=1-20

export REP=1

module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_02_evoquestions_parametereffs.R $SCRATCH/sens_routs/parametereffs_${SLURM_ARRAY_TASK_ID}_${REP}.Rout

#output is always to scratch
