#!/bin/sh

function add_binder {
    protein=$1
    binder=$2
    maxsol=$3
    index=$4
    for maxsol in 15 64; do 
            for binder in scyllo chiro; do 
                    genbox -cp $protein -ci ~/gro/${binder}_em.gro -nmol ${maxsol} -o abeta42_${binder}_${maxsol}_${index}.gro 
            done 
    done
}

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
        for repeat in {1..1}; do 
            mkdir sys${repeat}
            add_binder /Users/grace/scratch/abeta_pydr/abeta42.gro $binder $ratio $repeat
            mv *.gro sys${repeat}
        done
        cd ../
    done
    cd ../
done

