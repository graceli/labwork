#!/bin/sh

function check {
	if [ -e "common/nosol.tpr" ]; then
		echo "passed"
	else
		echo "common/nosol.tpr not found"
		return 1
	fi
	
	if [ -e "common/nonpolar.ndx" ]; then
		echo "passed"
	else
		echo "common/nonpolar.ndx not found"
		return 1
	fi
	
	if [ -e "common/pp_nonpolar.ndx" ]; then
		echo "passed"
	else
		echo "common/pp_nonpolar.ndx not found"
		return 1
	fi
	
	if [ -e "common/inositol_100mM_g_parse_index.ndx" ]; then
		echo "passed"
	else
		echo "common/inositol_100mM_g_parse_index.ndx found"
		return 1
	fi
	
}

function pp_nonpolar {
	DATA=.
	EXE=g_pp_nonpolar
	OUTPUT=.

	for file in `ls $1/*.xtc`; do
	       echo 14 | $EXE -f $file -s common/nosol.tpr -n common/pp_nonpolar.ndx -deffnm $file 2> /dev/null >&2 &

	       echo "running pp_nonpolar ... $n"

	       let n=n+1
	       if [ "$n" == 5 ]; then
	                wait
	       fi
	done

	mkdir pp_nonpolar
	mv *.dat pp_nonpolar
	tar cvfz pp_nonpolar.tgz pp_nonpolar
}

function residue_nonpolar {
	
	  echo `date` > residue_nonpolar.log
      EXE=g_inositol_residue_nonpolar_v2
      for file in `ls $1/*.xtc`; do
	      b=`basename $file`
	      echo 14 12 | $EXE -f $file -s ../common/nosol.tpr -n common/nonpolar.ndx -per_residue_contacts $b -per_inositol_contacts $b -dist 0.45 2>> residue_nonpolar.log >>&2 &

	      let n=$n+1
	      if [ "$n" == 5 ]; then
	               wait
	      fi
      done
       
      mkdir residue_nonpolar
      mv *.log *.xvg *.dat residue_nonpolar
      tar cvfz residue_nonpolar.tgz residue_nonpolar
      # rm -rf residue_nonpolar
}

function polar {
    EXE=/home/grace/bin/g_parse_index_oct21
    GRO=common/nosol.tpr
    NDX=common/inositol_100mM_g_parse_index.ndx
    CWD=$1

    for file in `ls $CWD/*_whole.xtc`; do
      b=`basename $file`
       seq 17 65 | $EXE -f $file -s $CWD/$GRO -n $CWD/$NDX -num_peptides 4 -num_inositol 45 -deffnm ${b}_ &


       let n=$n+1
       if [ "$n" == 5 ]; then
                wait
       fi
    done
        wait
 
    mkdir polar
    mv *.xvg *.dat polar
    tar cvfz polar.tgz polar
}

function do_dssp {
	DATA=.
	EXE=do_dssp
	OUTPUT=.

	echo `date`
	CWD=$PWD
	cd /dev/shm
	for file in `ls $CWD/*_whole.xtc`; do
		echo "do_dssp for $file"
		b=`basename $file`
		echo 1 | $EXE -f $file -s $CWD/nosol.tpr -sc /dev/shm/$b -o /dev/shm/$b -dt 2 2> /dev/null >&2 &


		let n=$n+1
		if [ "$n" == 5 ]; then
		         wait
		fi
	done

	wait

	cd $CWD
	mkdir dssp
	mv /dev/shm/*.xvg /dev/shm/*.xpm dssp
	tar cvfz dssp.tgz dssp

	trap "rm -rf /dev/shm/*; exit" TERM KILL INT
	
	
}

# assuming that all the needed files for analysis are there
# there are no unneeded xtcs in the directory
#variables are global

notpresent = check
if [ "$notpresent" == "1"]; then
	echo " Error "
	exit 0
fi

#residue_nonpolar $PWD
#wait

#polar $PWD
#wait