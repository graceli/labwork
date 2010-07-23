# Automated submission script
# Chris Neale November 2007
# chris.neale@utoronto.ca
# candrewn@gmail.com

echo -n "TIMING TEST (start): "
date
#!/bin/bash
MD="$1"
MYNP="$2"
MYMOL="$3"
TJUMP="$4"
NJUMPS="$5"
NFRAMES_IN_XTC="$6"
SAVE_FREQUENCY="$7"
NJUMPS_POSRE="$8"
NEVER_USE_SORT_SHUFFLE="$9"
CLUSTER="${10}"
DOUBLE_CHECK_DESORT="${11}"

case "$CLUSTER" in
bala|silky|coral)
  PATH=$PATH:/work/cneale/exe/gromacs-3.3.1/exec/bin
  ED=/work/cneale/exe/gromacs-3.3.1/exec/bin
  PED=${ED}
  ;;
greatwhite)
  PATH=$PATH:/home/cneale/exe/gromacs-3.3.1/exec/bin
  ED=/home/cneale/exe/gromacs-3.3.1/exec/bin
  PED=${ED}
  ;;
ccb)
  ED=/hpf/data/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/bin
  PED=${ED}
  mpiLocation="/tools/openmpi/1.2.1"
  mpirunProg="${mpiLocation}/bin/mpirun"
  mdrun_mpiProg="mdrun_openmpi_v1.2.1"

  ###############################################
  # For LAM mpi                                 #
  # mpiLocation="/tools/lam/lam-7.1.2"          #
  # mpirunProg="${mpiLocation}/bin/mpirun C"    #
  # mdrun_mpiProg="mdrun_mpi"                   #
  ###############################################

  ;;
ccb_aix)
  ED=/hpf/data/pomes/cneale/exe/gromacs-3.3.1_aix/exec/bin
  PED=${ED}
  mpiLocation="/usr"
  mpirunProg="${mpiLocation}/bin/poe"
  mdrun_mpiProg="mdrun_mpi"
 ;;
*)
  PATH=$PATH:/opt/sharcnet/gromacs/current/serial/bin:/opt/sharcnet/gromacs/current/mpi/bin
  ED=/opt/sharcnet/gromacs/current/serial/bin
  PED=/opt/sharcnet/gromacs/current/mpi/bin
  ;;
esac

case "$CLUSTER" in
bala|coral)
  DD=/work/cneale/exe/gromacs-3.3.1/exec/share/gromacs/template
  ;;
ccb)
  DD=/hpf/data/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/share/gromacs/template
  ;;
ccb_aix)
  DD=/hpf/data/pomes/cneale/exe/gromacs-3.3.1_aix/exec/share/gromacs/template
  ;;
greatwhite)
  DD=/home/cneale/exe/gromacs-3.3.1/exec/share/gromacs/template
  ;;
*)
  DD=/home/cneale/gromacs/template
  ;;
esac

# The NEW variable allows use of a special writing location (TMPDIR on ccb, could be /scratch/`whoami` on sharcnet)
case "$CLUSTER" in
ccb*)
  NEW="${TMPDIR}"
  ;;
*)
  NEW="."
  ;;
esac

case "$CLUSTER" in
ccb)
  CPP="\/home\/cneale\/exe\/cpp"
  ;;
ccb_aix)
  CPP="xlc -E"
  ;;
*)
  CPP="cpp"
  ;;
esac

PATH=$PATH:.

cd ${MD}

TINY_SLEEP=1
SHORT_SLEEP=10
LONG_SLEEP=60
EXTENDED_SLEEP=300

###############################################
# Startup tests:

if [ -e DO_NOT_RUN ]; then
  echo "ERROR error: file DO_NOT_RUN exists... exiting"
  exit
fi

for((v=0;v<2;v++)); do
  num=0;
  if [ -e finished_grompp ]; then
    let "num=$num+1"
  fi
  if [ -e finished_mdrun ]; then
    let "num=$num+1"
  fi
  if [ -e finished_desort ]; then
    let "num=$num+1"
  fi
  if [ -e finished_test ]; then
    let "num=$num+1"
  fi
  
  if((num==0)); then
    echo "Unsure how to start the run. Check this out"
    echo "$ls -l finished_grompp finished_mdrun finished_desort finished_test"
    ls -l finished_grompp finished_mdrun finished_desort finished_test
    echo "Perhaps you forgot to set a finished_XXX file upon starting your run?"
    echo "  - Otherwise there seems to be an error in the script."
  else 
    if((num!=1)); then
      echo "Unsure how to start the run. Check this out"
      echo "$ls -l finished_grompp finished_mdrun finished_desort finished_test"
      ls -l finished_grompp finished_mdrun finished_desort finished_test
      echo "Only one of these files should have existed."
    fi
  fi
  if((num==1)); then
    break; 
  fi
  echo "Will sleep then try one more time"
  sleep {$LONG_SLEEP}
done

if((num!=1)); then
  echo "Could not resolve multiple x for finished_x problem. Exiting"
  touch ./DO_NOT_RUN
  exit
fi

###############################################
# Initializations:

gromppProblemsInARow=0
reverted=0
MAX_CONSECUTIVE_GROMPP_ERRORS=2
MAX_CONSECUTIVE_MDRUN_ERRORS=2

function waitForExistNotEmpty {
  # First arg controls usage: 
  #   0  for existance test 
  #   1  for not empty test 
  #  -1  for must not exist test
  #   12 for not empty test plus require single value to equal third arg
  # Second arg is name of file/directory
  # Third arg is the expected single value in the file if First arg is = 12 
  #
  # Note: This overly complicated procedure is required for proper usage of the 
  #       CCB cluster where NFS delay can be significant and simple -s tests
  #       routinely fail to detect the fact that the file is empty as far as 
  #       val=`cat file` is concerned
  notEmpty="$1"
  case "$notEmpty" in
  1*)
    eneFlag="-s $2"
    ;;
  -1)
    eneFlag="! -e $2"
    ;;
  0)
    eneFlag="-e $2"
    ;;
  *)
    echo "ERROR error: incorrect argument to waitForExistNotEmpty = $notEmpty"
    exit
    ;;
  esac

  for((length=1;length<100;length++)); do
    if [ ${eneFlag} ]; then 
      break
    fi
    sleep ${length}
  done
  if((length==100)); then
    echo "ERROR error: Have slept for 30 minutes while waiting for the file $2 to meet conditions [ ${eneFlag} ] ... is there a problem in your script or is the NFS delay very very large?"
  fi

  # for all uses other than First Arg = 12 this function is over
  if((notEmpty!=12)); then
    return
  fi
  expectedVal="$3"
  for((length=1;length<100;length++)); do
    currentVal=`cat $2`
    if [ -n "$currentVal" ]; then
      # make sure variable is non-empty before making the comparison
      case "$currentVal" in
      $expectedVal)
        echo "NOTE: breaking from loop since currentVal($currentVal)=expectedVal($expectedVal)"
        break
        ;;
      esac
    fi
    sleep ${length}
  done
  if((length==100)); then
    echo "ERROR error: Have slept for 30 minutes while waiting for the file $2 to meet conditions `cat $2 = $expectedVal` ... is there a problem in your script or is the NFS delay very very large?"
  fi
}

###############################################
# The main loop:

for ((njump=0;njump<NJUMPS;njump++)); do

  # GROMPP (WITH SORTING)
  if [ -e finished_test ]; then
    NIN=`cat finished_test`
    TINIT=`cat finished_next_start_time`
    let "NOUT=$NIN+1"
    DIR=md${NOUT}_running
    PREV=md${NIN}_success
    if [ ! -e ${PREV} ]; then
      echo "There was some problem. Expected ${PREV} to exist, but it does not"
      touch ./DO_NOT_RUN
      exit
    fi
    if [ ! -e ${NEW}/${DIR} ]; then
      mkdir ${NEW}/${DIR}
    fi

    #nsteps=`echo "$TJUMP/0.002" | bc -l | awk -F '.' '{print $1}'`
    nsteps=`echo "$TJUMP/0.002" | bc` 

    #G: having NJUMPS_POSRE > 1 is probably not neccessary
    if((NOUT<=NJUMPS_POSRE)); then
      sed "s/TINIT/${TINIT}/" ${MYMOL}_posre.mdp | sed "s/NSTEPS/${nsteps}/" | sed "s/SAVE_FREQUENCY/${SAVE_FREQUENCY}/" | sed "s/CPP/${CPP}/" > ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp
      posreFlag="-r md0_success/${MYMOL}_md0_deshuffleddesorted.gro"
    else
      sed "s/TINIT/${TINIT}/" ${MYMOL}.mdp | sed "s/NSTEPS/${nsteps}/" | sed "s/SAVE_FREQUENCY/${SAVE_FREQUENCY}/" | sed "s/CPP/${CPP}/" > ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp
      posreFlag=""
    fi
    waitForExistNotEmpty 1 ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp
    if [ ! `tail -1 ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp | awk '{print $1}'` = ";EOF" ]; then sleep ${TINY_SLEEP}; fi
    if [ ! `tail -1 ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp | awk '{print $1}'` = ";EOF" ]; then sleep ${SHORT_SLEEP}; fi
    if [ ! `tail -1 ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp | awk '{print $1}'` = ";EOF" ]; then sleep ${LONG_SLEEP}; fi

    if((NIN!=0)); then
      if((NOUT<=NJUMPS_POSRE || NEVER_USE_SORT_SHUFFLE==1 || MYNP==1)); then
        # Can not do shuffle/sort with posre
        ${ED}/grompp -np ${MYNP} -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp ${posreFlag} -c ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -t ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.trr -p ${MYMOL}.top -n ${MYMOL}.ndx -e ${PREV}/${MYMOL}_md${NIN}.edr -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr
        rm mdout.mdp &
      else
        # Shuffle the .trr input file correctly. Assume that it is not currently shuffled
        ${ED}/grompp -np ${MYNP} -shuffle -sort -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp ${posreFlag} -c ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -p ${MYMOL}.top -n ${MYMOL}.ndx -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_a.tpr -deshuf ${NEW}/${DIR}/deshuffle_md${NOUT}_a.ndx 
        rm -f ${NEW}/${DIR}/deshuffle_md${NOUT}_a.ndx mdout.mdp &
        echo System | ${ED}/editconf -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}_a.tpr -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit_a.gro
        # g_desort -f original shuffled will unshuffle, therefore g_desort -f shuffled original will REshuffle
        ${DD}/g_desort -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit_a.gro ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -o ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}_a.ndx -n 6
        ${ED}/trjconv -f ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.trr -o ${PREV}/${MYMOL}_md${NIN}_reshuffleresort.trr -n ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}_a.ndx
        # Create the run input file
        ${ED}/grompp -np ${MYNP} -shuffle -sort -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp ${posreFlag} -c ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -t ${PREV}/${MYMOL}_md${NIN}_reshuffleresort.trr -p ${MYMOL}.top -n ${MYMOL}.ndx -deshuf ${NEW}/${DIR}/deshuffle_md${NOUT}.ndx -e ${PREV}/${MYMOL}_md${NIN}.edr -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr
        rm -f ${NEW}/${DIR}/deshuffle_md${NOUT}.ndx mdout.mdp &

        # In the future: implement this check. Note that this will require rethinking the _a postfixes
        #                since the files without the _a postfixes are the ones that I actually use
        #if((DOUBLE_CHECK_DESORT!=0)); then

          # If a new reshuffle.ndx file differs then the run is invalid.
          echo System | ${ED}/editconf -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro
          ${DD}/g_desort -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -o ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}.ndx -n 6
          case "$CLUSTER" in
          ccb_aix)
            # does not recognize the -q option
            look=`diff ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}_a.ndx ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}.ndx | head -1`
            ;;
          *)
            look=`diff -q ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}_a.ndx ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}.ndx`
            ;;
          esac
          if [ -n "$look" ]; then
            echo There was a big problem. ${NEW}/${DIR}/reshuffleresort_md${NOUT}_a.ndx and ${NEW}/${DIR}/reshuffleresort_md${NOUT}.ndx are different.
            mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr ${NEW}/${DIR}/${MYMOL}_md${NOUT}_notValid.tpr
          fi

        # End of the new test
        #fi

        #Create the deshuffle file to properly handle the next run
        ${DD}/g_desort -f ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro -o ${NEW}/${DIR}/deshuffledesort${MYMOL}_md${NOUT}.ndx -n 6
        rm -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}_a.tpr ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit_a.gro ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}_a.ndx ${PREV}/${MYMOL}_md${NIN}_reshuffleresort.trr ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro ${NEW}/${DIR}/reshuffleresort${MYMOL}_md${NOUT}.ndx mdout.mdp &
      fi  
    else
      if((NOUT<=NJUMPS_POSRE || NEVER_USE_SORT_SHUFFLE==1 || MYNP==1)); then
        # Can not do shuffle/sort with posre
        ${ED}/grompp -np ${MYNP} -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp ${posreFlag} -c ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -p ${MYMOL}.top -n ${MYMOL}.ndx -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr
        rm mdout.mdp &
      else
        ${ED}/grompp -np ${MYNP} -shuffle -sort -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.mdp ${posreFlag} -c ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro -p ${MYMOL}.top -n ${MYMOL}.ndx -deshuf ${NEW}/${DIR}/deshuffle_md${NOUT}.ndx -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr
        rm -f ${NEW}/${DIR}/deshuffle_md${NOUT}.ndx mdout.mdp &
        echo System | ${ED}/editconf -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro
        ${DD}/g_desort -f ${PREV}/${MYMOL}_md${NIN}_deshuffleddesorted.gro ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro -o ${NEW}/${DIR}/deshuffledesort${MYMOL}_md${NOUT}.ndx -n 6
        rm -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}_shuffledsortedInit.gro &
      fi
    fi

    if [ -e ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr ]; then
      echo ${NOUT} > finished_grompp
      waitForExistNotEmpty 12 finished_grompp ${NOUT}
      rm -f finished_test
      gromppProblemsInARow=0
    else
      # can not revert
      let "gromppProblemsInARow=$gromppProblemsInARow+1"
      if((gromppProblemsInARow>MAX_CONSECUTIVE_GROMPP_ERRORS)); then
        touch ./DO_NOT_RUN
        exit
      fi
    fi

  fi

  # MDRUN
  if [ -e finished_grompp ]; then
    NOUT=`cat finished_grompp`
    TINIT=`cat finished_next_start_time`
    DIR=md${NOUT}_running

    # Reversion is important in cases where a crash or time overrun leads to loss of data in ${TMPDIR}
    if [ ! -e ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr ]; then
      sleep ${SHORT_SLEEP}
      if [ ! -e ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr ]; then
        rm -f finished_grompp
        echo "ERROR error: finished_grompp existed for NOUT=$NOUT but ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr did not exist"
        if((reverted==0)); then
          echo "       reverting to finished_test"
          reverted=1
          let "NIN=$NOUT-1"
          echo "$NIN" > finished_test
          if [ ! -s finished_test ]; then sleep ${SHORT_SLEEP}; fi
          continue
        else
          touch ./DO_NOT_RUN
          exit
        fi
      else
        reverted=0;
      fi
    fi

    if((MYNP==1)); then
      returnValue=`${ED}/mdrun -deffnm ${NEW}/${DIR}/${MYMOL}_md${NOUT}`
    else
      case "$CLUSTER" in
      ccb)
        returnValue=`${mpirunProg} ${PED}/${mdrun_mpiProg} -np ${MYNP} -deffnm ${NEW}/${DIR}/${MYMOL}_md${NOUT}`
        ;;
      ccb_aix)
        export MP_INFOLEVEL=6
        returnValue=`${mpirunProg} ${PED}/${mdrun_mpiProg} -np ${MYNP} -deffnm ${NEW}/${DIR}/${MYMOL}_md${NOUT} -procs ${MYNP} -hfile /home/cneale/host.list`
        ;;
      greatwhite)
        returnValue=`prun -n ${MYNP} ${PED}/mdrun_mpi -np ${MYNP} -deffnm ${NEW}/${DIR}/${MYMOL}_md${NOUT}`
        ;;
      *)
        # assume sharcnet 
        returnValue=`/opt/hpmpi/bin/mpirun -srun -n ${MYNP} ${PED}/mdrun_mpi -np ${MYNP} -deffnm ${NEW}/${DIR}/${MYMOL}_md${NOUT}`
        ;;
      esac
    fi
    if((returnValue!=0)); then
      echo "ERROR error: mpirun for mdrun_mpi returned non-zero (${returnValue}). Exiting"
      exit
    fi
    echo ${NOUT} > finished_mdrun
    waitForExistNotEmpty 12 finished_mdrun ${NOUT}
    rm -f finished_grompp
  fi

  # DESORT
  if [ -e finished_mdrun ]; then
    NOUT=`cat finished_mdrun`
    TINIT=`cat finished_next_start_time`
    DIR=md${NOUT}_running

    # Reversion is important in cases where a crash or time overrun leads to loss of data in ${TMPDIR}
    if [ ! -e ${NEW}/${DIR}/${MYMOL}_md${NOUT}.xtc ]; then
      rm -f finished_mdrun
      echo "ERROR error: finished_mdrun existed for NOUT=$NOUT but ${NEW}/${DIR}/${MYMOL}_md${NOUT}.xtc did not exist"
      if((reverted==0)); then
        echo "       reverting to finished_test"
        reverted=1
        let "NIN=$NOUT-1"
        echo "$NIN" > finished_test
        waitForExistNotEmpty 12 finished_test ${NIN}
        continue
      else
        touch ./DO_NOT_RUN
        exit
      fi
    else
      reverted=0;
    fi

    if((NOUT<=NJUMPS_POSRE || NEVER_USE_SORT_SHUFFLE==1 || MYNP==1)); then
      # Can not do shuffle/sort with posre
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}.xtc ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}.trr ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.trr
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}.gro ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.gro
    else
      echo System | ${ED}/trjconv -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.xtc -s ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr -n ${NEW}/${DIR}/deshuffledesort${MYMOL}_md${NOUT}.ndx -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc &
      echo System | ${ED}/trjconv -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.trr -s ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr -n ${NEW}/${DIR}/deshuffledesort${MYMOL}_md${NOUT}.ndx -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.trr &
      echo System | ${ED}/trjconv -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}.gro -s ${NEW}/${DIR}/${MYMOL}_md${NOUT}.tpr -n ${NEW}/${DIR}/deshuffledesort${MYMOL}_md${NOUT}.ndx -o ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.gro &
      # Wait for all 3 desorts to finish
      wait
    fi

    echo ${NOUT} > finished_desort
    waitForExistNotEmpty 12 finished_desort ${NOUT}
    rm -f finished_mdrun
  fi

  # TEST
  if [ -e finished_desort ]; then
    TINIT=`cat finished_next_start_time`
    runHasNoErrors=1
    NOUT=`cat finished_desort`
    DIR=md${NOUT}_running

    # Reversion is important in cases where a crash or time overrun leads to loss of data in ${TMPDIR}
    if [ ! -e ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc ]; then
      rm -f finished_desort
      echo "ERROR error: finished_desort existed for NOUT=$NOUT but ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc did not exist"
      if((reverted==0)); then
        echo "       reverting to finished_test"
        reverted=1
        let "NIN=$NOUT-1"
        echo "$NIN" > finished_test
        waitForExistNotEmpty 12 finished_test ${NIN}
        continue
      else
        touch ./DO_NOT_RUN
        exit
      fi
    else
      reverted=0;
    fi

    ${ED}/gmxcheck -f ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc 2> ${NEW}/${DIR}/checkXTC
    # Ensure no magic number error
    magicNumberError=`grep Error ${NEW}/${DIR}/checkXTC | wc -l`
    if((magicNumberError==1));then
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted_magicNumberError.xtc
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.trr ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted_magicNumberError.trr
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.gro ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted_magicNumberError.gro
      runHasNoErrors=0
    fi
    ## Ensure expected number of frames is reached -- this is system and mdp file specific
    numFramesXTC=`grep '^Step' ${NEW}/${DIR}/checkXTC | awk '{print $2}'`
    if((numFramesXTC!=NFRAMES_IN_XTC)); then
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.xtc ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted_incompleteFrames.xtc
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.trr ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted_incompleteFrames.trr
      mv ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted.gro ${NEW}/${DIR}/${MYMOL}_md${NOUT}_deshuffleddesorted_incompleteFrames.gro
      runHasNoErrors=0
    fi
    if((runHasNoErrors)); then
      mv ${NEW}/${DIR} md${NOUT}_success
      waitForExistNotEmpty -1 ${NEW}/${DIR}
      while [ -e ${NEW}/${DIR} ]; do
        mv ${NEW}/${DIR} md${NOUT}_success
        sleep ${LONG_SLEEP}
      done
      
      echo ${NOUT} > finished_test
      nextTime=`echo "${TINIT}+${TJUMP}" | bc -l`
      echo ${nextTime} > finished_next_start_time
      waitForExistNotEmpty 12 finished_test ${NOUT}
      waitForExistNotEmpty 12 finished_next_start_time ${nextTime} 

      # Now do some clean up operations
      let "NCLEAN=$NOUT-2"
      if [ -e md${NCLEAN}_success ]; then
        if [ ! -e DATA ]; then
          mkdir DATA
        fi
        if((NCLEAN!=0)); then
          # Don't remove or modify the starting directory
          mkdir DATA/md${NCLEAN}_success
          mv md${NCLEAN}_success/${MYMOL}_md${NCLEAN}_deshuffleddesorted.xtc DATA/md${NCLEAN}_success
          waitForExistNotEmpty -1 md${NCLEAN}_success/${MYMOL}_md${NCLEAN}_deshuffleddesorted.xtc
          while [ -e md${NCLEAN}_success/${MYMOL}_md${NCLEAN}_deshuffleddesorted.xtc ]; do
            mv md${NCLEAN}_success/${MYMOL}_md${NCLEAN}_deshuffleddesorted.xtc DATA/md${NCLEAN}_success
            sleep ${LONG_SLEEP}
          done
          mv md${NCLEAN}_success/${MYMOL}_md${NCLEAN}.edr DATA/md${NCLEAN}_success
          waitForExistNotEmpty -1 md${NCLEAN}_success/${MYMOL}_md${NCLEAN}.edr
          while [ -e md${NCLEAN}_success/${MYMOL}_md${NCLEAN}.edr ]; do
            mv md${NCLEAN}_success/${MYMOL}_md${NCLEAN}.edr DATA/md${NCLEAN}_success
            sleep ${LONG_SLEEP}
          done
          rm -rf md${NCLEAN}_success &
        fi
      fi
    else
      # Send the run back to do the mdrun
      for((i=1;i<=MAX_CONSECUTIVE_MDRUN_ERRORS;i++)); do
        if [ ! -e md${NOUT}_failure${i} ]; then
          break
        fi
      done
      mv ${NEW}/${DIR} md${NOUT}_failure${i}
      if((i>MAX_CONSECUTIVE_MDRUN_ERRORS)); then
        echo "Too many failure for run $NOUT"
        touch ./DO_NOT_RUN
        exit
      fi
      # Send it back by setting as if grompp just finished
      mkdir ${NEW}/${DIR} 
      mv md${NOUT}_failure${i}/${MYMOL}_md${NOUT}.tpr md${NOUT}_failure${i}/deshuffledesort${MYMOL}_md${NOUT}.ndx ${NEW}/${DIR}
      echo ${NOUT} > finished_grompp
      waitForExistNotEmpty 12 finished_grompp ${NOUT}
    fi
    rm -f finished_desort
  fi

  wait
done 

wait


echo -n "TIMING TEST (end): "
date

