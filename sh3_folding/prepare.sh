#!/bin/sh

function make_check_dir() {
    mkdir runs_${NAME}
    NAME=$1
    START=$2
    END=$3
    for i in `seq $START $END`; do
        mkdir runs_${NAME}/${i}
        cp ${NAME}_ions_final.gro runs_${NAME}/${i}/${NAME}.gro
        cp ${NAME}.top ${NAME}/sys${i}/${i}/${NAME}.top
        grompp -f params/em.mdp -c  runs_${NAME}/${i}/${NAME}.gro -o runs_${NAME}/${i}/em.tpr -p runs_${NAME}/sys${i}/${i}/${NAME}.top
    done
}

make_systems "max_mrsd" 1 1
make_systems "max_q" 1 1