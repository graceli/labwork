# Automated submission script
# Chris Neale November 2007
# chris.neale@utoronto.ca
# candrewn@gmail.com

# this is a grace modified version of head.sh
# NFRAMES_IN_XTC should be computed using formula, not manually edited

#!/bin/bash
PATH=$PATH:.

##########################################################
# Basic setup options:
REPLICA=1
MYMOL=sheet_slab_water
MD=/work/grace
CLUSTER=requin
NUM_TO_SUBMIT=1
MYNP=4
# Time length of simulation (ps)
TJUMP=1
# Number of segments of the simulation
NJUMPS=1
SAVE_FREQUENCY=500
NJUMPS_POSRE=2
#Used to check whether finished xtc is incomplete
NFRAMES_IN_XTC=$(echo "$TJUMP/0.002+1" | bc)


##########################################################
# More advanced setup options:

RUN_AS_A_TEST=0
NEVER_USE_SORT_SHUFFLE=0
RUNTIME_LIMIT=1w
APPLY_RUNTIME_LIMIT=0
DOUBLE_CHECK_DESORT=0

##########################################################
# Things below this line do not usually need to be changed

cd ${MD}

if((RUN_AS_A_TEST)); then
  # Override the previously setup options
  MYNP=2
  TJUMP=0.2
  NJUMPS=1
  NFRAMES_IN_XTC=2
  SAVE_FREQUENCY=100
  NJUMPS_POSRE=0
  testFlag="--test"
else
  testFlag=""
fi

case "$CLUSTER" in
narwhal)
  timeFlag="-r $RUNTIME_LIMIT"
  ;;
*)
  # Not sure if this is correct for ccb or aix
  if((APPLY_RUNTIME_LIMIT==1)); then
    timeFlag="-r $RUNTIME_LIMIT"
  else
    timeFlag=""
  fi
  ;;
esac

case "$CLUSTER" in
ccb*)
  subCommand="qsub"
  ;;
*)
  # Assume sharcnet
  subCommand="sqsub"
  ;;
esac

if((MYNP>1)); then
  case "$CLUSTER" in
  ccb)
    if((MYNP==4)); then
      queueIdent="-pe ompi4p-* 4"
    else
      # note that a <4cpu job can be run on ompi4p-* queue but the other cpus will be wasted
      echo "unsure how to qsub for more than 4cpus"
      exit
    fi
    ;;
  ccb_aix)
    if((MYNP<=16)); then
      queueIdent="-l power -pe power16 ${MYNP}"
    else
      echo "Can not run more than 16 cpu on ccb_aix"
      exit
    fi
    ;;
  greatwhite)
    queueIdent="-q parallel --nompirun"
    ;;
  *)
    # Assume sharcnet
    queueIdent="-q mpi --nompirun"
    ;;
  esac
else
  case "$CLUSTER" in
  ccb_aix)
    queueIdent="-l power";;
  *)
    queueIdent="";;
  esac
fi

for((i=0;i<NUM_TO_SUBMIT; i++)); do
  if [ -e DO_NOT_RUN ]; then
    echo "ERROR: file DO_NOT_RUN exists... exiting"
    exit
  fi
  
  if [ -e wait_for_this_pid ]; then
    pid=`cat wait_for_this_pid`
    case "$CLUSTER" in
    ccb*)
      waitFlag="-hold_jid $pid"
      ;;
    *)
      waitFlag="-w $pid"
      ;;
    esac
  else
    waitFlag=""
  fi

  case "$CLUSTER" in
  ccb*)
    nameFlag="-N ${MYMOL}_${REPLICA}"
    mynpFlag=""
    subName=".${CLUSTER}"
    ;;
  *)
    nameFlag=""
    mynpFlag="-n ${MYNP}"
    subname=""
    ;;
  esac

  ${subCommand} ${testFlag} ${waitFlag} ${timeFlag} ${queueIdent} ${nameFlag} ${mynpFlag} -o `pwd`/out.${i} -e `pwd`/err.${i} /home/cneale/scripts/submission/sub${subName}.sh ${MD} ${MYNP} ${MYMOL} ${TJUMP} ${NJUMPS} ${NFRAMES_IN_XTC} ${SAVE_FREQUENCY} ${NJUMPS_POSRE} ${NEVER_USE_SORT_SHUFFLE} ${CLUSTER} ${DOUBLE_CHECK_DESORT} > tmp.sub 2>&1
  case "$CLUSTER" in
  ccb*)
    pid=`tail -n 1 tmp.sub | awk -F '"' '{print $2}' | awk '{print $1}'`
    ;;
  silky|greatwhite|gulper)
    # Silky return value is slightly different (although this may change in the near future)
    # looks like this: Job <852684> is submitted to queue <mpi>.
    pid=`tail -n 1 tmp.sub | awk -F '<' '{print $2}' | awk -F '>' '{print $1}'`
    ;;
  *)
    # Regular Sharcnet looks like this: submitted as jobid 222166
    pid=`grep jobid tmp.sub | tail -n 1 | awk '{print $4}'`
    ;;
  esac
  
  echo "Submitting ${subCommand} ${testFlag} ${waitFlag} ${timeFlag} ${queueIdent} ${nameFlag} ${mynpFlag} -o `pwd`/out.${i} -e `pwd`/err.${i} /home/cneale/scripts/submission/sub${subName}.sh ${MD} ${MYNP} ${MYMOL} ${TJUMP} ${NJUMPS} ${NFRAMES_IN_XTC} ${SAVE_FREQUENCY} ${NJUMPS_POSRE} ${NEVER_USE_SORT_SHUFFLE} ${CLUSTER} ${DOUBLE_CHECK_DESORT} -- received pid = $pid"
  echo "$pid" > wait_for_this_pid
  sleep 1

done
