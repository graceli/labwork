import os
import subprocess
import logging
import string

logging.basicConfig(filename='analysis.log', level=logging.DEBUG)

# TODO: Implement these using configuration. For now leave them as python global variables for convenience

class Analysis(object):
    def __init__(self):
        if not os.path.exists('/dev/shm/grace'):
            os.mkdir('/dev/shm/grace')

    RATIO = 15
    ISO = "scyllo"
    DATA_BASE_DIR = "{0}/inositol/abeta42/current/{1}/{2}_nonsolvent".format(os.environ['SCRATCH'], RATIO, ISO)
    OUTPUT_BASE_DIR = os.environ['PWD']
    START = 0
    FINISH = 9

    logging.debug("RATIO=%s", RATIO)
    logging.debug("ISO=%s", ISO)
    logging.debug("DATA_BASE_DIR=%s", DATA_BASE_DIR)
    logging.debug("OUTPUT_BASE_DIR=%s", OUTPUT_BASE_DIR)
    logging.debug("START=%s, FINISH=%s", START, FINISH)


    # This is the equivalent of performing teardown. But it might be better to perform this clean up action
    # when the analysis object is destroyed
    def _clean(self, file_tag):
        shell_command = "cd /dev/shm/grace; tar cvfz analysis_{0}.tgz * --remove-files; cp analysis_{0}.tgz {1}; rm /dev/shm/grace/analysis_{0}.tgz".format(file_tag, self.OUTPUT_BASE_DIR)
        logging.debug("Cleaning up ...")
        for cmd in shell_command.split(";"):
            logging.debug(cmd)

        try:
            # Copy files back to base_dir
            subprocess.check_call(shell_command, shell=True)
            # Deletes the files in /dev/shm/grace
            # subprocess.check_call("rm -rf /dev/shm/grace", shell=True)
        except subprocess.CalledProcessError as error:
            # If the subprocess call fails log the error
            logging.debug("Command %s finished with retcode=%s and output %s", 
                         error.cmd, error.returncode, error.output)
        except OSError as error:
            print "Command finished with system error errno=", error.errno 

class NonpolarContactAnalysis(Analysis):
    def _command(self, chain_first, chain_last, insgrp_first, insgrp_last, xtc_filename, system_index, output_dir, testing=False):
        insgroups = ' '.join([str(i) for i in range(insgrp_first, insgrp_last+1)])
        test_command = ''
        if testing:
            test_command = '-b 0 -e 1000'
 
        shell_command = string.Template("seq $chain_first $chain_last | parallel -j 16 \"echo {} $index_groups | g_inositol_residue_nonpolar_v2 -f $data_dir/$xtc_name -s ${iso}_${ratio}_nosol.tpr -n ab_${ratio}_${iso}_nonpolar_revised.ndx -per_residue_contacts $output_dir/ab_${iso}_${ratio}_${sys_index}_chain{}_residue_np_contact.dat -per_inositol_contacts $output_dir/ab_${iso}_${ratio}_${sys_index}_chain{}_inositol_np_contact.dat -per_residue_table $output_dir/ab_${iso}_${ratio}_${sys_index}_chain{}_table.dat -num_inositols $num_ligands $testing\"").substitute(chain_first=chain_first, chain_last=chain_last, insgrp_first=insgrp_first, insgrp_last=insgrp_last, data_dir=self.DATA_BASE_DIR, xtc_name=xtc_filename, iso=self.ISO, ratio=self.RATIO, output_dir=output_dir, sys_index=system_index, num_ligands=self.RATIO, index_groups=insgroups, testing=test_command)

        return shell_command
 
    def analyze_contact_by_ligand(self, iso, ratio, num_ligands, output_dir, testing=False):
        logging.debug("Analyzing nonpolar contact for each inositol")

        # Index of the first protein chain
        chain_first = 0
        # Index of the second protein chain
        chain_last = 4

        # Index of the first inositol group
        insgrp_first = 5
        # Index of the second inositol group
        insgrp_last = 19
        
        if testing:
            logging.debug("Running command as a test")

        try:
            for sys_index in range(self.START, self.FINISH):
                xtc = str(sys_index) + "_final.xtc"
                xtc_full_path = os.path.join(self.DATA_BASE_DIR, xtc)
                if os.path.exists(xtc_full_path):
                    shell_command = self._command(chain_first, chain_last, insgrp_first, insgrp_last, xtc, sys_index, output_dir, testing=testing)
                    logging.debug("Running: %s ", shell_command)
                    print "Running", shell_command
                    subprocess.check_call(shell_command, shell=True)
                else:
                    logging.debug("The file %s does not exist", xtc)
        except subprocess.CalledProcessError as error:
            # If the subprocess call fails log the error
            print "Command", error.cmd, "finished with retcode=", error.returncode, "and output", error.output
        except OSError as error:
            print "Command finished with system error errno=", error.errno

        self._clean("nonpolar_contacts_%(iso)s" % vars())


class HBondContactAnalysis(Analysis):
    def _command_by_ligand(self, xtc_filename, system_index, output_dir, testing=False):
        test_command = ''
        if testing:
            test_command = '-b 0 -e 100'
 
        shell_command = string.Template("seq 1 ${ratio} | parallel -j 16 \"echo {} Protein | g_hbond -f $data_dir/$xtc_name -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}_by_ligand.ndx -nonitacc -nomerge -num $output_dir/ab_${iso}_${ratio}_${sys_index}_ins{} $testing > /dev/null 2>&1\"").substitute(data_dir=self.DATA_BASE_DIR, xtc_name=xtc_filename, iso=self.ISO, ratio=self.RATIO, output_dir=output_dir, sys_index=system_index, testing=test_command)

        return shell_command

    def _command_by_residue(self, xtc_filename, system_index, output_dir, testing=False):
        res_start=0
        res_end=129
        ligand_grp=130

        test_command = ''
        if testing:
            test_command = '-b 0 -e 100'
 
        shell_command = string.Template("seq $res_start $res_end | parallel -j 8 \"echo {} $ligand_grp | g_hbond -f $data_dir/$xtc_name -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_by_residue.ndx -nonitacc -nomerge -num $output_dir/ab_${iso}_${ratio}_${sys_index}_residue{} $testing > /dev/null 2>&1\"").substitute(res_start=res_start, res_end=res_end, ligand_grp=ligand_grp, data_dir=self.DATA_BASE_DIR, xtc_name=xtc_filename, iso=self.ISO, ratio=self.RATIO, output_dir=output_dir, sys_index=system_index, testing=test_command)

        return shell_command

    # Calculates the number of hydrogen bonds made with each residue.
    # g_bond is ran for each ligand molecule
    def analyze_contact_by_ligand(self, iso, ratio, num_ligands, output_dir, testing=False):
        logging.debug("Analyzing nonpolar contact for each inositol")
        
        if testing:
            logging.debug("Running command as a test")

        try:
            for sys_index in range(self.START, self.FINISH):
                xtc = str(sys_index) + "_final.xtc"
                xtc_full_path = os.path.join(self.DATA_BASE_DIR, xtc)
                
                if os.path.exists(xtc_full_path):
                    shell_command = self._command_by_ligand(xtc, sys_index, output_dir, testing=testing)

                    logging.debug("Running: %s ", shell_command)
                    print "Running", shell_command
                    
                    subprocess.check_call(shell_command, shell=True)
                else:
                    logging.debug("The file %s does not exist", xtc)
        except subprocess.CalledProcessError as error:
            # If the subprocess call fails log the error
            print "Command", error.cmd, "finished with retcode=", error.returncode, "and output", error.output
        except OSError as error:
            print "Command finished with system error errno=", error.errno

        self._clean("hbonds_%(iso)s_by_ligand" % vars())

    # Calculates the number of hydrogen bonds made with each residue.
    # g_hbond is ran each time for each residue-ligands in the protein
    def analyze_contact_by_residue(self, iso, ratio, num_ligands, output_dir, testing=False):
        logging.debug("Analyzing nonpolar contact for each inositol")
        
        if testing:
            logging.debug("Running command as a test")

        try:
            for sys_index in range(self.START, self.FINISH):
                xtc = str(sys_index) + "_final.xtc"
                xtc_full_path = os.path.join(self.DATA_BASE_DIR, xtc)
                
                if os.path.exists(xtc_full_path):
                    shell_command = self._command_by_residue(xtc, sys_index, output_dir, testing=testing)

                    logging.debug("Running: %s ", shell_command)
                    print "Running", shell_command
                    
                    subprocess.check_call(shell_command, shell=True)
                else:
                    logging.debug("The file %s does not exist", xtc)
        except subprocess.CalledProcessError as error:
            # If the subprocess call fails log the error
            print "Command", error.cmd, "finished with retcode=", error.returncode, "and output", error.output
        except OSError as error:
            print "Command finished with system error errno=", error.errno

        self._clean("hbonds_%(iso)s_by_residue" % vars())


# # This is a bit of a hack to get g_hbond to work.
# res_start=0
# res_end=129
# INS_grp=130
# num=0
# function hbonds {
#     iso=$1
#     ratio=$2
#     output_dir=$3/hbonds
#     mkdir -p $output_dir
  
#     for s in `seq 0 9`; do
#       xtc="${s}_final"
#       if [ -e "$DATA/${xtc}.xtc" ]; then
#           mkdir -p $output_dir/$s
#           seq $res_start $res_end | parallel -j 8 "echo {} $INS_grp | g_hbond -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}.ndx -nonitacc -nomerge -num $output_dir/$s/{} $xvgr $TEST > /dev/null 2>&1"
#       fi
#       # python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
#     done

#     # Cleaning up
#     clean "${iso}_${ratio}_hbonds"
# }



# function hbonds_inositol {
#     iso=$1
#     ratio=$2
#     output_dir=$3/hbonds_inositol
#     mkdir -p $output_dir
  
#     for s in `seq 0 9`; do
#       xtc="${s}_final"
#       if [ -e "$DATA/${xtc}.xtc" ]; then
#           mkdir -p $output_dir/$s
#           seq 1 $ratio | parallel -j 8 "echo {} Protein | g_hbond -f $DATA/$xtc -s ${iso}_${ratio}_nosol.tpr -n g_hbond_inositol_${iso}_${ratio}.ndx -nonitacc -nomerge -num $output_dir/$s/{} $xvgr $TEST > /dev/null 2>&1"
#       fi
#       # python /home/grace/AnalysisScripts/abeta_analysis/abeta_analysis.py sys${s}.h5
#     done

#     # Cleaning up
#     clean "${iso}_${ratio}_hbonds_inositol"
# }




if __name__ == "__main__":
    hbonds = HBondContactAnalysis()
    print "starting nonpolar analysis"
    scratch = Analysis.OUTPUT_BASE_DIR + "/test"
    temp = '/dev/shm/grace'
    hbonds.analyze_contact_by_ligand(Analysis.ISO, Analysis.RATIO, Analysis.RATIO, temp, testing=True)
    hbonds.analyze_contact_by_residue(Analysis.ISO, Analysis.RATIO, Analysis.RATIO, temp, testing=True)

    nonpolar = NonpolarContactAnalysis()
    nonpolar.analyze_contact_by_ligand(Analysis.ISO, Analysis.RATIO, Analysis.RATIO, temp, testing=True)
    


