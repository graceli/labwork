#!/bin/sh
##PBS -l nodes=1:ppn=8,walltime=5:00:00
#PBS -N volmap

#cd $PBS_O_WORKDIR

XTC=scyllo_64_volmap_all_fit.xtc
GRO=scyllo_64_nosol.gro

if [ "$is_cer" == "1" ]; then
	vmd -dispdev text -e volmap_cer.tcl -args $XTC $GRO
else
	vmd -dispdev text -e volmap.tcl -args $XTC $GRO
fi

