#!/bin/sh

NFILES=1024
RUN_DATE=`date "+%F"`

bonnie -p2 -u root
nohup bonnie++ -d /3ware/c0u0 -s 4000 -n $NFILES -m "t2:/3ware/c0u0" -r 2000 -ys -u root >> /3ware/u0_sync_0_${NFILES}_${RUN_DATE}.csv &
nohup bonnie++ -d /3ware/c0u1 -s 4000 -n $NFILES -m "t2:/3ware/c0u1" -r 2000 -ys -u root >> /3ware/u1_sync_1_${NFILES}_${RUN_DATE}.csv &

