#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=5:00:00
#PBS -N volmap

trap 'clean; exit $?' TERM KILL EXIT SIGINT

function clean {
    cp /dev/shm/grace/* .
    rm -rf /dev/shm/grace
}

mkdir /dev/shm/grace

cd $PBS_O_WORKDIR

if [ "$is_cer" == "1" ]; then
    vmd -dispdev text -e volmap_cer.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap.log 2>&1
else
    vmd -dispdev text -e volmap.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap.log 2>&1 
fi
