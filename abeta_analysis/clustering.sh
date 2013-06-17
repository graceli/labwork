#cluster size calculation with g_clustsize
function cluster {
    GRP=1 #protein
    task=0
    xtc=$1
    tpr=$2
    ndx=$3
    output_dir=$4/cluster
    mkdir -p $output_dir
    for file in `ls $xtc/*.xtc`; do
        base=`basename $file .xtc`
        echo $GRP | g_clustsize -f $file -s $tpr -n $ndx -o $output_dir/${base}_o -ow $output_dir/${base}_ow -nc $output_dir/${base}_nclust.xvg -mc $output_dir/${base}_mc.xvg -ac $output_dir/${base}_ac.xvg -hc $output_dir/${base}_hc.xvg -cut $CUTOFF -noxvgr $TEST > $OUTPUT 2>&1 &

        let task=$task+1
        if [ "$task" == "16" ]; then
            wait
            task=0
        fi
    done
}


cd $PBS_O_WORKDIR

ratio=$1
sys=$2
isomer=$3

xtc="$ratio/${isomer}_nonsolvent"
tpr="${isomer}_${ratio}_nosol.tpr"
ndx="${isomer}_${ratio}_nosol.ndx"
output_base="/dev/shm/analysis"

cluster $xtc $tpr $ndx $output_base
