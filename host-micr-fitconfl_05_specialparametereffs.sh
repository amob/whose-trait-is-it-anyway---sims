#!/bin/bash
#SBATCH --time=0:20:00
#SBATCH --job-name=R6fitconevo_feedbacks_hp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca
#SBATCH --array=1-280

#280 in redo.
#ultimately we want array to be 1-324, since that is how many rows in parameters
export REP=6

module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_05_specialparametereffs.R $SCRATCH/feedback_routs/specialparams_${SLURM_ARRAY_TASK_ID}_${REP}.Rout

#output is always to scratch
