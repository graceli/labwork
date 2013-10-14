#!/bin/sh
#PBS -q qfbb@mp2 
#PBS -l walltime=24:00:00 -l nodes=1:ppn=1
#PBS -N pgab_analysis


# This script runs in /mnt/scratch_mp2/pomes/ligrace1/pgab

cd $PBS_O_WORKDIR

trap "clean; exit 0" EXIT SIGINT INT KILL

function clean {
	cp -rp /dev/shm/grace/* $PBS_O_WORKDIR/volmap/
	rm -rf /dev/shm/grace
}

python volmap.py

# vmd volmap command
mkdir -p /dev/shm/grace/volmap

XTC=volmap/pgab_37-50_all_dt10_fit.xtc
GRO=pgab/pgab.gro
vmd -dispdev text -e ~/labwork/pgab/volmap.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_full.log 2>&1 

