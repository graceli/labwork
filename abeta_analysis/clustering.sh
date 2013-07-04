#!/bin/sh
#PBS -l nodes=1:ppn=8,walltime=08:00:00
set -x

#cluster size calculation with g_clustsize
function cluster {
    #TEST='-b 0 -e 100'
    GRP=135 # ligand or inositol isomer
    xtc=$1
    tpr=$2
    ndx=$3
    output_dir=$4
    mkdir -p $output_dir
    for file in `ls $xtc/*.xtc`; do
        base=`basename $file .xtc`
        echo $GRP | $HOME/bin/clustering -f $file -s $tpr -n $ndx -o $output_dir/${base}_o -ow $output_dir/${base}_ow -nc $output_dir/${base}_nclust.xvg -mc $output_dir/${base}_mc.xvg -ac $output_dir/${base}_ac.xvg -hc $output_dir/${base}_hc.xvg -clust-info $output_dir/${base}_clust_info.dat -cut $CUTOFF $TEST > /dev/null 2>&1 &
    done
    wait
}


cd $PBS_O_WORKDIR

CUTOFF=0.35
ratio=$RATIO
isomer=$ISOMER

xtc="$ratio/${isomer}_nonsolvent"
tpr="${isomer}_${ratio}_nosol.tpr"
ndx="g_hbond_${ratio}_${isomer}_by_residue.ndx"

mkdir ${isomer}_${ratio}_clust_${CUTOFF}
cluster $xtc $tpr $ndx ${isomer}_${ratio}_clust_${CUTOFF}

