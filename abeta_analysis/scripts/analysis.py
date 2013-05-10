import os
import subprocess
import logging

logging.basicConfig(filename='analysis.log', level=logging.DEBUG)

# TODO: Implement these using configuration. For now leave them as python global variables for convenience
RATIO = 15
ISO = "scyllo"
DATA_BASE_DIR = "{0}/inositol/abeta42/current/{1}/{2}_nonsolvent".format(os.environ['SCRATCH'], RATIO, ISO)
OUTPUT_BASE_DIR = os.environ['pwd']
START = 0
FINISH = 9

class Analysis(object):
    def _clean(file_tag):
        shell_commands = "cd /dev/shm/grace;tar cvfz analysis_%(file_tag)s.tgz *;cp analysis_%(file_tag)s.tgz %(BASE_DIR)s;rm -rf /dev/shm/grace" % vars()

        for cmd in shell_commands.split(";"):
            logging.debug(cmd)

        try:
            # Copy files back to base_dir
            subprocess.check_call("cd /dev/shm/grace", shell=True)
            subprocess.check_call("tar cvfz analysis_%(file_tag)s.tgz *" % vars(), shell=True)
            subprocess.check_call("cp analysis_%(file_tag)s.tgz %(BASE_DIR)s" % vars(), shell=True)
            # Deletes the files in /dev/shm/grace
            subprocess.check_call("rm -rf /dev/shm/grace", shell=True)
        except subprocess.CalledProcessError as error:
            # If the subprocess call fails log the error
            logging.debug("Command {0} finished with retcode={1} and output {2}", 
                         (error.cmd, error.returncode, error.output))
        except OSError as error:
            print "Command finished with system error errno=", error.errno 

class NonpolarContactAnalysis(Analysis):
    def analyze_contact_by_inositol(self, iso, ratio, num_ligands, output_dir, testing=False):
        # Index of the first protein chain
        chain_first = 0
        # Index of the second protein chain
        chain_last = 4

        # Index of the first inositol group
        insgrp_first = 1
        # Index of the second inositol group
        insgrp_last = 2

        shell_command = "seq %(chain_first)s %(chain_last)s | parallel -j 16 \"seq %(insgrp_first)s %(insgrp_last)s | g_inositol_residue_nonpolar_v2 -f %(DATA)s/{0} -s %(iso)s_%(ratio)s_nosol.tpr -n ab_%(ratio)s_%(iso)s_nonpolar_revised.ndx -per_residue_contacts %(output_dir)s/{1}_chain{}_residue_np_contact.dat -per_inositol_contacts %(output_dir)s/{1}_chain{}_inositol_np_contact.dat -per_residue_table %(output_dir)s/{1}_chain{}_table.dat -num_inositols %(num_ligands)s\"" % vars()

        if testing:
            logging.debug("Running command as a test")
            shell_command += " -dt 0 100"

        try:
            for sys_index in range(START, FINISH):
                xtc = sys_index + "_final"
                if os.exists(xtc):
                    shell_command.format(xtc, sys_index)
                    logging.debug("Running: " + shell_command)
                    subprocess.check_call(shell_command, shell=True)
                else:
                    logging.debug("The file {0} does not exist", (xtc))

        except subprocess.CalledProcessError as error:
            # If the subprocess call fails log the error
            print "Command", error.cmd, "finished with retcode=", error.returncode, "and output", error.output
        except OSError as error:
            print "Command finished with system error errno=", error.errno

        _clean("hbonds_by_inositol")

#   function nonpolar_revised {
#     iso=$1
#     ratio=$2
#     output_dir=$3/nonpolar_revised
#     mkdir -p $output_dir

#     # We don't need make the index over and over again for the peptides, because they are all the same
#      #echo -e "'SideChain'&aC*&!rACE\nsplitch16\nq" | make_ndx -f ${iso}_${ratio}_nosol.tpr -o ab_${ratio}_nonpolar.ndx

#     for s in `seq 0 9`; do
#         xtc="${s}_final"
#         if [ -e "$DATA/${xtc}.xtc" ]; then
#             seq $chain1 $chain5 | parallel -j 5 "seq $insgrp_first $insgrp_second | g_inositol_residue_nonpolar_v2 -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_${iso}_nonpolar_revised.ndx -per_residue_contacts $output_dir/${s}_chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/${s}_chain{}_inositol_np_contact.dat -per_residue_table $output_dir/${s}_chain{}_table.dat -num_inositols $num_inos $TEST"
#         fi
#     done
#     clean "${iso}_${ratio}_nonpolar"
# }

def __main__():
    analysis = NonpolarContactAnalysis()
    analysis.analyze_contact_by_inositol(ISO, RATIO, RATIO, OUTPUT_BASE_DIR + "/test", 
        testing=True)
