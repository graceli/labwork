#cluster size calculation with g_clustsize
function cluster {
    GRP=135 # ligand or inositol isomer
    xtc=$1
    tpr=$2
    ndx=$3
    output_dir=$4
    mkdir -p $output_dir
    for file in `ls $xtc/*.xtc`; do
        base=`basename $file .xtc`
        echo $GRP | clustering -f $file -s $tpr -n $ndx -o $output_dir/${base}_o -ow $output_dir/${base}_ow -nc $output_dir/${base}_nclust.xvg -mc $output_dir/${base}_mc.xvg -ac $output_dir/${base}_ac.xvg -hc $output_dir/${base}_hc.xvg -clust-info $output_dir/${base}_clust_info.dat -cut $CUTOFF $TEST > $OUTPUT 2>&1 &
    done
}


cd $PBS_O_WORKDIR

CUTOFF=0.35
ratio=$1
sys=$2
isomer=$3

xtc="$ratio/${isomer}_nonsolvent"
tpr="${isomer}_${ratio}_nosol.tpr"
ndx="g_hbond_${ratio}_${isomer}_by_residue.ndx"

mkdir ${isomer}_clust_${CUTOFF}
cluster $xtc $tpr $ndx ${isomer}_clust
