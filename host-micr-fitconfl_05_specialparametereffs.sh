#!/bin/bash
#SBATCH --time=0:25:00
#SBATCH --job-name=R5fitconevo_feedbacks_hp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.obrien@utoronto.ca
#SBATCH --array=1-280

# we want array to be 1-280, since that is how many rows in parameters
export REP=5

module load r/3.4.3-anaconda5.1.0
R CMD BATCH --no-restore --no-save $HOME/whosetrait/host-micr-fitconfl_05_specialparametereffs.R $SCRATCH/feedback_routs/specialparams_${SLURM_ARRAY_TASK_ID}_${REP}.Rout

#output is always to scratch
