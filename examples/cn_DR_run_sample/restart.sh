#!/bin/bash 
#PBS -l nodes=1:compute-eth:ppn=8,walltime=48:00:00,os=centos53computeA
#PBS -N SRV_vt0.0x0

### server=$(cat serverip); port=$(cat port); /project/pomes/cneale/GPC/exe/DR/DR_version2.1.8_extended/bin/DR_commander $server $port EXIT

SCRIPT_NAME=00
EXED=/project/pomes/cneale/GPC/exe/DR/DR_version2.1.8_extended/bin
SNAPSHOT=00.1277252602.snapshot

###########################################################
### Things below this line don't usually need to be changed


# If not an interactive job (i.e. -I), then cd into the directory where I typed qsub.
if [ "$PBS_ENVIRONMENT" != "PBS_INTERACTIVE" ]; then
  if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
  fi
fi

date=`date`
CPD=`echo PRIOR_TO_RESTART_${date}|sed "s/ /_/g"`
mkdir ${CPD}
mv output outerr `whoami` ${SCRIPT_NAME}.log package.tar.gz ${CPD}
cp ${SCRIPT_NAME}.script ${SCRIPT_NAME}.forcedatabase ${SNAPSHOT} ${CPD}

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
nohup ${EXED}/DR_server ${SCRIPT_NAME}.script -s ${SNAPSHOT}


