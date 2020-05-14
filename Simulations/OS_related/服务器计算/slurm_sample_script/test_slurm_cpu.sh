#! /bin/bash

#=== part 1 ===

#SBATCH --partition=bldesign
#SBATCH --account=bldesign
#SBATCH --qos=blregular
#SBATCH --job-name=test_cpu
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=2048
#SBATCH --output=/scratchfs/heps/yangym/test_slurm/job-%j.out

#=== part 2 ===
source /afs/ihep.ac.cn/users/z/zhaohf/.bash_profile
python3.7 /scratchfs/heps/yangym/test_oasys/MICRON_ME.py
