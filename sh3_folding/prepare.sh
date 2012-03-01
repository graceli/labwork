#!/bin/sh

function make_systems {
    NAME=$1
    START=$2
    END=$3
	
	if [ ! -e "runs_$NAME" ]; then
		mkdir runs_${NAME}
	fi

    for i in `seq $START $END`; do
        mkdir runs_${NAME}/${i}
        cp initial/${NAME}_ions_final.gro runs_${NAME}/${i}/${NAME}_${i}.gro
        cp initial/${NAME}.top runs_${NAME}/${i}/${NAME}_${i}.top
        grompp -f params/em.mdp -c  runs_${NAME}/${i}/${NAME}_${i}.gro -o runs_${NAME}/${i}/em.tpr -p runs_${NAME}/${i}/${NAME}_${i}.top -po ${i}
    done
}

make_systems "max_mrsd" 1 100
make_systems "max_q" 1 100
