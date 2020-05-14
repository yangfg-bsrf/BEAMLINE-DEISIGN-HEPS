#! /bin/bash

#=====================================================
#===== Modify the following options for your job =====
#=====    DON'T remove the #! /bin/bash lines    =====
#=====      DON'T comment #SBATCH lines          =====
#=====        of partition,account and           =====
#=====                qos                        =====
#=====================================================

# Specify the partition name from which resources will be allocated  
#SBATCH --partition=mbh

# Specify which expriment group you belong to.
# This is for the accounting, so if you belong to many experiments,
# write the experiment which will pay for your resource consumption
#SBATCH --account=mbh

# Specify which qos(job queue) the job is submitted to.
#SBATCH --qos=regular

#=====================================================
#===== Modify the following options for your job =====
#=====   The following options are not mandatory =====
#=====================================================

# Specify your job name, e.g.: "test_job" is my job name
#SBATCH --job-name=test_job

# Specify how many cores you will need, e.g.: 16
#SBATCH --ntasks=16

# Specify the output file path of your job 
# Attention!! It's a output file name, not a directory
# Also you must have write access to this file
# An output file name example is : job-%j.out, where %j will be replace with job id
#SBATCH --output=/path/to/your/output/file

#=================================
#===== Run your job commands =====
#=================================

#=== You can define your variables following the normal bash style
VAR1="value1"
VAR2="value2"

#=== Run you job commands in the following lines
# srun is necessary, with srun, Slurm will allocate the resources you are asking for

# You can run an executable script with srun
# modify /path/to/your/script to a real file path name
srun /path/to/your/script

# Or if your program is written with MPI, you can run it with mpiexec 
# First, run a simple command with srun
srun -l hostname

# later, you can run your MPI program with mpiexec
# The output will be written under the path specified by the --output option
# modify /path/to/your/mpi_program to your real program file path
mpiexec /path/to/your/mpi_program
