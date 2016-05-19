#!/bin/bash
#SBATCH -p standard
#SBATCH -t 1-01:00:00
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem-per-cpu=1gb
#SBATCH -J RP_CIRC3
#SBATCH -o RP_CIRC3_%j
#SBATCH --mail-type=all
module load matlab
srun matlab -nodesktop -r circ3
