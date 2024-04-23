#!/bin/bash
#SBATCH -C v100
#SBATCH -A vermaaslab
##SBATCH --array=0-80%1
#SBATCH --gres=gpu:2
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:0:0 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=barkarar@msu.edu
#SBATCH --job-name=5pw

cd $SLURM_SUBMIT_DIR
#modules are loaded automatically by the NAMD module.
module use /mnt/home/vermaasj/modules
module load NAMD/3.0b5-gpu
NUM=`ls run*dcd | wc -l` 
PRINTNUM=`printf "%03d" $NUM`
srun namd3 ++ppn 2 +devices 0,1 run.namd > run-${PRINTNUM}.log

