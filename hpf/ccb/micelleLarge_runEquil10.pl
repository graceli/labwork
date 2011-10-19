#!/bin/bash
LD_LIBRARY_PATH=/usr/lib64/tls:/hpf/tools/n1ge/lib/lx24-amd64:/tools/gcc/4.1.1/lib64:/tools/openmpi/1.2.1/lib:/tools/local/lib
PATH=/hpf/tools/n1ge/bin/lx24-amd64:/usr/local/bin:/usr/bin:/usr/X11R6/bin:/bin:/usr/games:/opt/gnome/bin:/tools/local/bin:/home/cneale/exe:/projects/pomes/cneale/exe/gromacs-3.3.1/exec/fftw-3.1.2/bin:.:/tools/openmpi/1.2.1/bin

### for lam-mpi LD_LIBRARY_PATH includes /tools/lam/lam-7.1.2/lib
### for openmpi_v1.2.1 LD_LIBRARY_PATH includes tools/openmpi/1.2.1/lib and /tools/local/lib
### for lam-mpi PATH includes /tools/lam/lam-7.1.2/bin
### for openmpi_v1.2.1 PATH includes /tools/openmpi/1.2.1/bin

## The hpf8p queue is no longer in use, it was used like this:
## qsub -pe hpf8p 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh

MD=/projects/pomes/cneale/micelle/2large
cd ${MD}

MYMOL=sys225
MYSYS=largeM

TINIT=153700
TJUMP=200
MYNP=4

### Upon regular running of this script, must set a file ./finished in order to start right away
while [ ! -e ./finished ]; do
  sleep 60
done
rm -f ./finished ./gromppfinished
while [ -e ./finished ]; do sleep 5; done
while [ -e ./gromppfinished ]; do sleep 5; done

### This is the major loop
for ((NIN=658; NIN<1000; ++NIN)); do
  NOUT=$(expr $NIN + 1)
  IN=md${NIN}
  OUT=md${NOUT}

  ### Prepare the MYNP variable
  if [ ! -e noeight ]; then
    npindicator=`qstat -g c | grep ompi8p2.q | awk '{print $3}'`
    if((npindicator==0)); then
      ## extra check to give the larger system priority on NP=8
      sleep 360
      npindicator=`qstat -g c | grep ompi8p2.q | awk '{print $3}'`
      if((npindicator==0)); then
        MYNP=8
      else
        MYNP=4
      fi
    else
      MYNP=4
    fi
  else
    MYNP=4
  fi

  ### Create the input files for the next job, make executable, and wait for completion
  sed "s/TINIT/${TINIT}/" ${MYMOL}.mdp > ${MYMOL}_${OUT}.mdp
  sed "s/_out/_${OUT}/" ${MYSYS}_template_grompp | sed "s/_in/_${IN}/" | sed "s/_mynp/${MYNP}/" > ${MYSYS}_${OUT}_grompp
  sed "s/_out/_${OUT}/" ${MYSYS}_template_mdrun.sh | sed "s/_in/_${IN}/" | sed "s/_mynp/${MYNP}/" > ${MYSYS}_${OUT}_mdrun.sh
  while [ ! -s ${MYSYS}_${OUT}_grompp ]; do sleep 5; done
  while [ ! `tail -1 ${MYSYS}_${OUT}_grompp | awk '{print $1}'` = "###EOF" ]; do sleep 5; done
  chmod +x ${MYSYS}_${OUT}_grompp
  while [ ! -s ${MYSYS}_${OUT}_mdrun.sh ]; do sleep 5; done
  while [ ! `tail -1 ${MYSYS}_${OUT}_mdrun.sh | awk '{print $1}'` = "###EOF" ]; do sleep 5; done
  chmod +x ${MYSYS}_${OUT}_mdrun.sh
  while [ ! -x ${MYSYS}_${OUT}_grompp ]; do sleep 5; done
  while [ ! -x ${MYSYS}_${OUT}_mdrun.sh ]; do sleep 5; done
  while [ ! -s ${MYMOL}_${OUT}.mdp ]; do sleep 5; done
  while [ ! `tail -1 ${MYMOL}_${OUT}.mdp | awk '{print $1}'` = ";EOF" ]; do sleep 5; done

  ### This appears to be the only way to be sure!
  sleep 30

  ### Ensure necessary files exist
  if [ ! -e ${MYMOL}_${IN}_deshuffleddesorted.gro ]; then
    echo "./${MYMOL}_${IN}_deshuffleddesorted.gro did not exist!" >> ./ERROR_MESSAGE
    exit
  fi
  if [ ! -e ${MYMOL}_${IN}_deshuffleddesorted.trr ]; then
    echo "./${MYMOL}_${IN}_deshuffleddesorted.trr did not exist!" >> ./ERROR_MESSAGE
    exit
  fi
  if [ ! -e ${MYMOL}_${IN}.edr ]; then
    echo "./${MYMOL}_${IN}.edr did not exist!" >> ./ERROR_MESSAGE
    exit
  fi
  if [ ! -e ${MYMOL}.top ]; then
    echo "./${MYMOL}.top did not exist!" >> ./ERROR_MESSAGE
    exit
  fi
  if [ ! -e ${MYMOL}.ndx ]; then
    echo "./${MYMOL}.ndx did not exist!" >> ./ERROR_MESSAGE
    exit
  fi

  ### Run the grompp preprocessor and ensure output files exist
  gottpr=0
  while [ "$gottpr" = 0 ]; do
    /home/cneale/scripts/conditionalAllowSub.sh
    mydate=`date`
    if [ ${MYNP} = 16 ]; then
      qsub -pe ompi16p-4-31 1 ${MD}/${MYSYS}_${OUT}_grompp
      echo "executing qsub -pe ompi16p-4-31 1 ${MD}/${MYSYS}_${OUT}_grompp at ${mydate}" >> ./qsubbed.info
    else
      if [ ${MYNP} = 8 ]; then
        qsub -pe ompi8p2 1 ${MD}/${MYSYS}_${OUT}_grompp
        echo "executing qsub -pe ompi8p2 1 ${MD}/${MYSYS}_${OUT}_grompp at ${mydate}" >> ./qsubbed.info
      else
        if [ ${MYNP} = 4 ]; then
          qsub -pe ompi4p-* 1 ${MD}/${MYSYS}_${OUT}_grompp
          echo "executing qsub -pe ompi4p-* 1 ${MD}/${MYSYS}_${OUT}_grompp at ${mydate}" >> ./qsubbed.info
        else
          if [ ${MYNP} = 1 ]; then
            qsub ${MD}/${MYSYS}_${OUT}_grompp
            echo "executing qsub ${MD}/${MYSYS}_${OUT}_grompp at ${mydate}" >> ./qsubbed.info
          else
            echo "The value of MYNP=${MYNP} is not valid" >> ./ERROR_MESSAGE
            exit
          fi
        fi
      fi
    fi
    while [ ! -e ./gromppfinished ]; do
      sleep 10
    done
    rm -f ./gromppfinished
    while [ -e ./gromppfinished ]; do sleep 5; done
    k=0  
    while [ ! -e ./${MYMOL}_${OUT}.tpr ]; do
      sleep 60
      k=$(expr $k + 1)
      if [ "$k" -gt 10 ]; then
        break
      fi
    done
    if [ -e ./${MYMOL}_${OUT}.tpr ]; then
      gottpr=1
    else
      echo "./${MYMOL}_${OUT}.tpr was not created!" >> ./ERROR_MESSAGE
    fi
  done
  if [ ! -e ./${MYMOL}_${OUT}.tpr ]; then
    echo "./${MYMOL}_${OUT}.tpr was not created!" >> ./ERROR_MESSAGE
    exit
  fi

  mkdir ${IN}_success
  while [ ! -e ${IN}_success ]; do sleep 5; done
  mv ${MYMOL}*${IN}* ${MYSYS}_*${IN}* output* deshuffledesort${MYMOL}_${IN}.ndx ${IN}_success
  # Remove the shuffled and sorted gro trr and xtc :: DONT'T DO THIS IF ITS NOT A SHUFFLED RUN
  rm -f ${IN}_success/${MYMOL}_${IN}.gro ${IN}_success/${MYMOL}_${IN}.xtc ${IN}_success/${MYMOL}_${IN}.trr ${IN}_success/${MYMOL}_${IN}.tpr

  # The log files are moved by ${MYMOL}*${IN}*
  #Don't copy the .ndx file, it is too large
  cp *.top ${IN}_success

  ### Ensure that the finished file does not exits prior to submission
  if [ -e ./finished ]; then rm -f ./finished; fi
  while [ -e ./finished ]; do sleep 5; done

  ### Submit the job and recycle
  /home/cneale/scripts/conditionalAllowSub.sh
  mydate=`date`
  if [ ${MYNP} = 16 ]; then
    qsub -pe ompi16p-4-31 16 ${MD}/${MYSYS}_${OUT}_mdrun.sh
    echo "executing qsub -pe ompi16p-4-31 16 ${MD}/${MYSYS}_${OUT}_mdrun.sh at ${mydate}" >> ./qsubbed.info
  else
    if [ ${MYNP} = 8 ]; then
      qsub -pe ompi8p2 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh
      echo "executing qsub -pe ompi8p2 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh at ${mydate}" >> ./qsubbed.info
    else
      if [ ${MYNP} = 4 ]; then
        qsub -pe ompi4p-* 4 ${MD}/${MYSYS}_${OUT}_mdrun.sh
        echo "executing qsub -pe ompi4p-* 4 ${MD}/${MYSYS}_${OUT}_mdrun.sh at ${mydate}" >> ./qsubbed.info
      else
        if [ ${MYNP} = 1 ]; then
          qsub ${MD}/${MYSYS}_${OUT}_mdrun.sh
          echo "executing qsub ${MD}/${MYSYS}_${OUT}_mdrun.sh at ${mydate}" >> ./qsubbed.info
        else
          echo "The value of MYNP=${MYNP} is not valid" >> ./ERROR_MESSAGE
          exit
        fi
      fi
    fi
  fi
  ### Wait until the last job has completed
  while [ ! -e ./finished ]; do
    sleep 60
  done
  rm -f ./finished
  while [ -e ./finished ]; do sleep 5; done

  ### Failure recovery 1
  if [ ! -e ${MYMOL}_${OUT}.gro ]; then
    echo "./${MYMOL}_${OUT}.gro did not exist!" >> ./ERROR_MESSAGE
    echo "There must have been an error in the run, resubmitting first recovery" >> ./ERROR_MESSAGE
    mkdir ${OUT}_failure1
    mv ${MYMOL}*${OUT}* *log output* ${OUT}_failure1
    cp ${MYSYS}_*${OUT}* ${OUT}_failure1
    cp *.ndx *.top ${OUT}_failure1
    cp ${OUT}_failure1/${MYMOL}_${OUT}.tpr ./
    /home/cneale/scripts/conditionalAllowSub.sh
    mydate=`date`
    if [ ${MYNP} = 16 ]; then
      qsub -pe ompi16p-4-31 16 ${MD}/${MYSYS}_${OUT}_mdrun.sh
      echo "executing qsub -pe ompi16p-4-31 16 ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 1) at ${mydate}" >> ./qsubbed.info
    else
      if [ ${MYNP} = 8 ]; then
        qsub -pe ompi8p2 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh
        echo "executing qsub -pe ompi8p2 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 1) at ${mydate}" >> ./qsubbed.info
      else
        if [ ${MYNP} = 4 ]; then
          qsub -pe ompi4p-* 4 ${MD}/${MYSYS}_${OUT}_mdrun.sh
          echo "executing qsub -pe ompi4p-* 4 ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 1) at ${mydate}" >> ./qsubbed.info
        else
          if [ ${MYNP} = 1 ]; then
            qsub ${MD}/${MYSYS}_${OUT}_mdrun.sh
            echo "executing qsub ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 1) at ${mydate}" >> ./qsubbed.info
          else
            echo "The value of MYNP=${MYNP} is not valid" >> ./ERROR_MESSAGE
            exit
          fi
        fi
      fi
    fi

    while [ ! -e ./finished ]; do
      sleep 60
    done
    rm -f ./finished
    while [ -e ./finished ]; do sleep 5; done
  fi

  ### Failure recovery 2
  if [ ! -e ${MYMOL}_${OUT}.gro ]; then
    echo "./${MYMOL}_${OUT}.gro did not exist!" >> ./ERROR_MESSAGE
    echo "There must have been an error in the run, resubmitting second recovery" >> ./ERROR_MESSAGE
    mkdir ${OUT}_failure2
    mv ${MYMOL}*${OUT}* *log output* ${OUT}_failure2
    cp ${MYSYS}_*${OUT}* ${OUT}_failure2
    cp *.ndx *.top ${OUT}_failure2
    cp ${OUT}_failure2/${MYMOL}_${OUT}.tpr ./
    /home/cneale/scripts/conditionalAllowSub.sh
    mydate=`date`
    if [ ${MYNP} = 16 ]; then
      qsub -pe ompi16p-4-31 16 ${MD}/${MYSYS}_${OUT}_mdrun.sh
      echo "executing qsub -pe ompi16p-4-31 16 ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 2) at ${mydate}" >> ./qsubbed.info
    else
      if [ ${MYNP} = 8 ]; then
        qsub -pe ompi8p2 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh
        echo "executing qsub -pe ompi8p2 8 ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 2) at ${mydate}" >> ./qsubbed.info
      else
        if [ ${MYNP} = 4 ]; then
          qsub -pe ompi4p-* 4 ${MD}/${MYSYS}_${OUT}_mdrun.sh
          echo "executing qsub -pe ompi4p-* 4 ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 2) at ${mydate}" >> ./qsubbed.info
        else
          if [ ${MYNP} = 1 ]; then
            qsub ${MD}/${MYSYS}_${OUT}_mdrun.sh
            echo "executing qsub ${MD}/${MYSYS}_${OUT}_mdrun.sh (Failure Recovery 2) at ${mydate}" >> ./qsubbed.info
          else
            echo "The value of MYNP=${MYNP} is not valid" >> ./ERROR_MESSAGE
            exit
          fi
        fi
      fi
    fi
    while [ ! -e ./finished ]; do
      sleep 60
    done
    rm -f ./finished
    while [ -e ./finished ]; do sleep 5; done
  fi

  ### Failure Quit
  if [ ! -e ${MYMOL}_${OUT}.gro ]; then
    echo "./${MYMOL}_${OUT}.gro did not exist!" >> ./ERROR_MESSAGE
    mkdir ${OUT}_failure3
    mv ${MYMOL}*${OUT}* *log output* ${OUT}_failure3
    cp ${MYSYS}_*${OUT}* ${OUT}_failure3
    cp *.ndx *.top ${OUT}_failure3
    cp ${OUT}_failure3/${MYMOL}_${OUT}.tpr ./
    echo "There must have been an error in the run, NOT resubmitting third recovery" >> ./ERROR_MESSAGE
    rm -f ./finished
    while [ -e ./finished ]; do sleep 5; done
    exit
  fi

  TINIT=`cnadd $TINIT $TJUMP | awk -F "." '{print $1}'`
  #TINIT=`cnadd $TINIT $TJUMP | awk '{print $1}'`

done 


