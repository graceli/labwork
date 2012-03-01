#!/bin/sh

mkdir runs_mrsd
mkdir runs_q

for i in `seq 1 1`; do
    mkdir runs_mrsd/sys${i}
    cp max_mrsd_ions_final.gro max_mrsd.top runs_mrsd/sys${i}

    mkdir runs_q/sys${i}
    cp max_q_ions_final.gro max_q.top runs_q/sys${i}
done