#!/bin/bash
 
# This script is named: data-offload.sh
#PBS -l walltime=72:00:00
#PBS -q archive
#PBS -N offload
#PBS -j oe
#PBS -me
 
trap "echo 'Job script not completed';exit 129" TERM INT
# individual tarballs already exist
 
/usr/local/bin/hsi  -v <<EOF1
mkdir put-away
cd put-away
cput $SCRATCH/workarea/finished-job1.tar.gz : finished-job1.tar.gz
end
EOF1
status=$?
if [ ! $status == 0 ];then
   echo 'HSI returned non-zero code.'
   /scinet/gpc/bin/exit2msg $status
   exit $status
else
   echo 'TRANSFER SUCCESSFUL'
fi
 
/usr/local/bin/hsi  -v <<EOF2
mkdir put-away
cd put-away
cput $SCRATCH/workarea/finished-job2.tar.gz : finished-job2.tar.gz
end
EOF2
status=$?
if [ ! $status == 0 ];then
   echo 'HSI returned non-zero code.'
   /scinet/gpc/bin/exit2msg $status
   exit $status
else
   echo 'TRANSFER SUCCESSFUL'
fi
 
trap - TERM INT