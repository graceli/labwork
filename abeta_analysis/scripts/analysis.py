
# These are the analysis configuration variables. For now leave them as global variables for convenience
DATA = data_dir
class Hbonds:
    def __init__():
        pass

    def analyze_hbonds_by_inositol(iso, ratio, output_dir):
       shell_command = "seq $chain1 $chain5 | parallel -j 5 "seq $insgrp_first $insgrp_second | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_${iso}_nonpolar_revised.ndx -per_residue_contacts $output_dir/${s}_chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/${s}_chain{}_inositol_np_contact.dat -per_residue_table $output_dir/${s}_chain{}_table.dat -num_inositols $num_inos $TEST"
        
        for index in range(START, FINISH):
            xtc = index + "_final" 
            if os.exists(xtc):
                # execute shell command
                os.comand(shell_command)
        # clean up analysis -- what is this?
 
 #   function nonpolar_revised {
    iso=$1
    ratio=$2
    output_dir=$3/nonpolar_revised
    mkdir -p $output_dir

    # We don't need make the index over and over again for the peptides, because they are all the same
     #echo -e "'SideChain'&aC*&!rACE\nsplitch16\nq" | make_ndx -f ${iso}_${ratio}_nosol.tpr -o ab_${ratio}_nonpolar.ndx

    for s in `seq 0 9`; do
        xtc="${s}_final"
        if [ -e "$DATA/${xtc}.xtc" ]; then
            seq $chain1 $chain5 | parallel -j 5 "seq $insgrp_first $insgrp_second | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_${iso}_nonpolar_revised.ndx -per_residue_contacts $output_dir/${s}_chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/${s}_chain{}_inositol_np_contact.dat -per_residue_table $output_dir/${s}_chain{}_table.dat -num_inositols $num_inos $TEST"
        fi
    done
    clean "${iso}_${ratio}_nonpolar"
}
