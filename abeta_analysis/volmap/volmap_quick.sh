#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=8:00:00
#PBS -N volmap_quick

set -x

trap 'clean; exit $?' TERM KILL EXIT SIGINT

function clean {
    cd /dev/shm/grace
    tar cvfz $PBS_O_WORKDIR/volmap.tar.gz *
    rm -rf /dev/shm/grace
}

mkdir /dev/shm/grace

cd $PBS_O_WORKDIR

for r in 15 64; do
    for iso in scyllo chiro glycerol; do
        XTC=${iso}_${r}_volmap_all_fit_smaller.xtc
        GRO=${iso}_${r}_nosol.gro

        if [ "$iso" == "glycerol" ]; then
            vmd -dispdev text -e volmap_cer.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_quick.log 2>&1 &
        else
            vmd -dispdev text -e volmap.tcl -args $XTC $GRO > /dev/shm/grace/${XTC}_volmap_quick.log 2>&1 &
        fi
    done
    wait
done
