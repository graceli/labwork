#!/bin/bash 
#PBS -l nodes=1:compute-eth:ppn=8,walltime=48:00:00,os=centos53computeA
#PBS -N SRV_vt0.0x0

SCRIPT_NAME=00
EXED=/project/pomes/cneale/GPC/exe/DR/DR_version2.1.8_extended/bin

###########################################################
### Things below this line don't usually need to be changed


# If not an interactive job (i.e. -I), then cd into the directory where I typed qsub.
if [ "$PBS_ENVIRONMENT" != "PBS_INTERACTIVE" ]; then
  if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
  fi
fi


rm -f ./port ./serverip
port=$(grep PORT ./${SCRIPT_NAME}.script | awk '{print $2}')
echo "$port" > ./port
serverip=$(/sbin/ifconfig eth0 | grep 'inet addr:' | cut -d: -f2 | awk '{print $1}')
echo "$serverip" > ./serverip
rm -f ./package.tar.gz
tar cfv ./package.tar ./port ./serverip ./DR_client_body
gzip ./package.tar

mkdir -p output/cleanup
mkdir -p output/data
mkdir -p output/outerr
mkdir `whoami`

## Don't put it in the background if it is submitted
nohup ${EXED}/DR_server ${SCRIPT_NAME}.script 


