#PBS -q qfbb@mp2 
#PBS -l walltime=08:00:00 -l nodes=1:ppn=1
#PBS -N pgab_analysis

trap 'rm -rf /dev/shm/*; exit 0' EXIT SIGINT INT TERM KILL

cd $PBS_O_WORKDIR

mkdir /dev/shm/grace

python ~/labwork/abeta_analysis/trjcat/trjcat_project.py pgab_glucosamine 1 13 -n nonsolvent -o pgab_glucosamine_nonsolvent

