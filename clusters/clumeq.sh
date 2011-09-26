# / Bin / bash

#
# Sample script for SGE job submission
#
# - The options for EMS is provided with
# Lines starting with # $
#
# For more information, see the qsub manual
#

# You must first load the necessary modules
module load compilers/intel/11.1.059
module load mpi/openmpi/1.4.3_intel

# The name of the task. It is advisable not to put spaces or
# Special characters. Otherwise, put between "".
# $-N JOB

# Project (Rap ID) to which this task is assigned. Use the command
Colossus-info # if you do not know this information.
# $-P xxx-yyy-zz 

# Parallel environment and number of cores
# $-Pe default 16 # 8 cores per node. Must be a multiple of 8.

# Queue Implementation
# $-Q short
# # $ Q-med
# # $-Q test

# Estimated execution time (in seconds)
# $-L = 300 h_rt

# File in which to save the stdout of the script.
# If not specified, will be recorded in <NomDeLaTache>. <IdTache> O
# # $-O $ HOME / job.out

# File in which to save the stderr output of the script.
# If not specified, will be recorded in <NomDeLaTache>. <IdTache> E
# # $-E $ HOME / job.err

# Shell to use for the task
# $-S / bin / bash

# Run the task from the current directory
# $-Cwd

# List of users who email
# # $-M user@domain.com

# When a mail is sent there?
# "B" = when the task begins
# "E" = when the job ends
# "A" = when the task looks
# "S" = when the task is suspended
# "N" = do not send email
# # $-M bea

# Send signal SIGUSR2 several seconds
# Before killing the task.
# # $-Notify
# Trap "" USR2

# Task batch. It is specified start-end: SKIP
The variable $ # contains the number SGE_TASK_ID
# Assigned to a given process.
# # $-T 1-10:1

# Number of threads for OpenMP tasks. It is
# Normally not necessary, since by default,
# 8 threads will be executed, and you will not want
# Run in a different number.
# Export OMP_NUM_THREADS = 8

# Command to execute
# / CLUMEQ / example / hello-openmp
mpirun / CLUMEQ / example / hello-mpi
