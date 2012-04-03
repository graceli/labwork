#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=8:00:00
#PBS -N volmap_full

set -x

trap 'clean "${iso}_${r}"; exit $?' TERM KILL EXIT SIGINT

function clean {
    cd /dev/shm/grace
    tar cvfz $PBS_O_WORKDIR/${1}_volmap_full_dt10.tar.gz *
    rm -rf /dev/shm/grace
}

mkdir /dev/shm/grace

cd $PBS_O_WORKDIR

XTC=${iso}_${r}_volmap_all_fit2.xtc
GRO=${iso}_${r}_nosol.gro

if [ "$iso" == "glycerol" ]; then
    vmd -dispdev text -e volmap_cer.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_full.log 2>&1 &
else
    vmd -dispdev text -e volmap.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_full.log 2>&1 &
fi

wait
