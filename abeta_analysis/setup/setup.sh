#!/bin/bash

function debug {
    echo "$1"
}

function add_binder {
    ERROR='add_binder.log'
    protein=$1
    binder=$2
    maxsol=$3
    repeat=$4
    
    if [ "$binder" == "water" ]; then
        #cmd="genbox -cp $protein -cs ~/gro/tip3.gro -o abeta42_${binder}_${repeat}.gro"
        cmd="mv $protein abeta42_${binder}_${repeat}.gro"
        return="abeta42_${binder}_${repeat}"
    else
        cmd="genbox -cp $protein -ci /home/grace/gro/${binder}_em.gro -nmol ${maxsol} -o abeta42_${binder}_${maxsol}_${repeat}.gro"
        return="abeta42_${binder}_${ratio}_${repeat}"
    fi

    $cmd 2>> $ERROR >&2
    echo $return
}

function add_ions {
    ERROR='add_ions.log'
    top=$1
	protein=$2

	# resolvate and add ions
	echo "genbox -cp $protein -cs /home/grace/gro/tip3.gro -p $top -o out1.gro"
	genbox -cp $protein -cs /home/grace/gro/tip3.gro -p $top -o out1.gro 2>> $ERROR >&2
	touch empty.mdp
	
	echo "grompp -f empty.mdp -c out1.gro -p $top -o out1.tpr"
	grompp -f empty.mdp -c out1.gro -p $top -o out1.tpr 2>> $ERROR >&2
	echo "genion -s out1.tpr -conc 0.15 -pname NA+ -nname CL- -o ions -p $top"
	echo "SOL" | genion -s out1.tpr -conc 0.15 -pname NA+ -nname CL- -o ions -p $top 2>> $ERROR >&2

	# add counter ions
	grompp -f empty.mdp -c ions.gro -p $top -o ions.tpr 2>> $ERROR >&2
    echo "genion -s ions.tpr -pname NA+ -nname CL- -np 10 -nn 0 -o ions_counter -p $top"
	echo "SOL" | genion -s ions.tpr -pname NA+ -nname CL- -np 10 -nn 0 -o ions_counter -p $top 2>> $ERROR >&2
}

structure=$1
if [ -z $structure ]; then
    echo "Error: please provide structure gro file name"
    exit;
fi

# parameters
MAX_REPEAT=20
base_dir=`pwd`
#/Users/grace/scratch/abeta_pydr/abeta42.gro

editconf -f $structure -o start.gro -box 8 8 8 
for ratio in 15 64; do 
    if [ ! -e $ratio ]; then
        mkdir $ratio
    fi
    cd $ratio
    for binder in scyllo chiro glycerol water; do
        if [ ! -e $binder ]; then
            mkdir $binder
        fi

        cd $binder
        for repeat in `seq 1 $MAX_REPEAT`; do 
            mkdir sys${repeat}
            
            if [ "$binder" == "water" ]; then
                cp $base_dir/abeta42_${binder}.top abeta42_${binder}_${repeat}.top
            else
                cp $base_dir/abeta42_${binder}_${ratio}.top abeta42_${binder}_${ratio}_${repeat}.top
            fi

			cp $base_dir/start.gro  start${repeat}.gro
		    cp $base_dir/vdwradii.dat .	
			echo "add_binder start${repeat}.gro $binder $ratio $repeat"
            top=`add_binder start${repeat}.gro $binder $ratio $repeat`

			echo "top file $top"
			echo "add_ions $top $top"
			add_ions $top $top 

            cp *.gro *.top *.mdp *.tpr sys${repeat}
            rm \#*
        done
        cd ../
    done
    cd ../
done

