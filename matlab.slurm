#!/bin/bash
#SBATCH -p standard
#SBATCH -t 3-01:00:00
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem-per-cpu=1gb
#SBATCH -J RP_CIRC1
#SBATCH -o RP_CIRC1_%j
#SBATCH --mail-type=all
module load matlab
srun matlab -nodesktop -r circ1
