#! /bin/bash

#=== part 1 ===

#SBATCH --partition=bldesign
#SBATCH --account=bldesign
#SBATCH --qos=blregular
#SBATCH --job-name=B1_Undulator
#SBATCH --ntasks=100
#SBATCH --mem-per-cpu=2048
#SBATCH --output=/scratchfs/heps/yangfg/job.log

#=== part 2 ===
date
echo $SLURM_JOB_NODELIST
source /afs/ihep.ac.cn/users/y/yangfg/.bash_profile
time srun python /scratchfs/heps/yangfg/B1_ME.py
date
