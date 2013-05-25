import os
import subprocess
import logging
import string

logging.basicConfig(filename=os.path.join(os.environ['PWD'], 'analysis_' + str(os.environ['PBS_JOBID']) + '.log'), level=logging.DEBUG)

# TODO: Implement these using configuration. For now leave them as python global variables for convenience

class Analysis(object):

    def __init__(self, iso, ratio):
        self.iso = iso
        self.ratio = ratio
        self.DATA_BASE_DIR = "{0}/inositol/abeta42/current/{1}/{2}_nonsolvent".format(os.environ['SCRATCH'], ratio, iso)
        self.START = 0
        self.FINISH = 10
        self.TEMP_DIR = '/dev/shm/grace'

        logging.debug("RATIO=%s", ratio)
        logging.debug("ISO=%s", iso)
        logging.debug("DATA_BASE_DIR=%s", self.DATA_BASE_DIR)
        logging.debug("START=%s, FINISH=%s", self.START, self.FINISH)
        logging.debug("Initializing temp dir %s", self.TEMP_DIR)

        try:
            os.mkdir(self.TEMP_DIR)
        except OSError as error:
            logging.debug("Directory %s already exists", self.TEMP_DIR)

    # This is the equivalent of performing teardown. But it might be better to perform this clean up action
    # when the analysis object is destroyed
    def _clean(self, file_tag, output_abs_path):
        shell_command = "cd {2}; tar cvfz analysis_{0}.tgz * --remove-files; cp analysis_{0}.tgz {1}; rm {2}/analysis_{0}.tgz".format(file_tag, output_abs_path, self.TEMP_DIR)

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
            logging.debug("Command %s finished with retcode=%s and output %s", error.cmd, error.returncode, error.output)
        except OSError as error:
            print "Command finished with system error errno=", error.errno 

class NonpolarContactAnalysis(Analysis):
    def __init__(self, iso, ratio, num_atoms, abs_path):
        super(NonpolarContactAnalysis, self).__init__(iso, ratio)
        self.output_abs_path = abs_path
        self.num_atoms = num_atoms
        logging.debug("Analysis results will be outputted to %s", self.output_abs_path) 
        if not os.path.exists(self.output_abs_path):
            logging.debug("%s does not exist. Creating ...", self.output_abs_path)
            os.makedirs(self.output_abs_path)

    def _command(self, insgrp_first, insgrp_last, xtc_filename, system_index, testing=False):
        insgroups = ' '.join([str(i) for i in range(insgrp_first, insgrp_last+1)])
        test_command = '-b 0 -e 250000'
        if testing:
            test_command = '-b 0 -e 1000'
 
        shell_command = string.Template("seq 0 ${insgrp_last} | g_inositol_residue_nonpolar_v2 -f $data_dir/$xtc_name -s ${iso}_${ratio}_nosol.tpr -n nonpolar_${iso}_${ratio}.ndx -per_residue_contacts $output_dir/ab_${iso}_${ratio}_${sys_index}_residue_np_contact.dat -per_inositol_contacts $output_dir/ab_${iso}_${ratio}_${sys_index}_inositol_np_contact.dat -per_residue_table $output_dir/ab_${iso}_${ratio}_${sys_index}_table.dat -num_inositols $num_ligands -num_atoms $num_atoms $testing").substitute(insgrp_last=insgrp_last, data_dir=self.DATA_BASE_DIR, xtc_name=xtc_filename, iso=self.iso, ratio=self.ratio, num_atoms=self.num_atoms, output_dir='/dev/shm/grace', sys_index=system_index, num_ligands=self.ratio, index_groups=insgroups, testing=test_command)

        return shell_command
 
    def analyze_contact_by_ligand(self, testing=False):
        logging.debug("Analyzing nonpolar contact for each inositol")

        # Index of the first protein chain
        chain_first = 0
        # Index of the second protein chain
        chain_last = 4

        # Index of the first inositol group
        insgrp_first = 1
        # Index of the second inositol group
        insgrp_last = 64
        
        if testing:
            logging.debug("Running command as a test")

        try:
            for sys_index in range(self.START, self.FINISH):
                xtc = str(sys_index) + "_final.xtc"
                xtc_full_path = os.path.join(self.DATA_BASE_DIR, xtc)
                logging.debug("Looking for file %s", xtc_full_path)
                if os.path.exists(xtc_full_path):
                    shell_command = self._command(insgrp_first, insgrp_last, xtc, sys_index, testing=testing)
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

        self._clean("nonpolar_contacts_{0}_{1}".format(self.iso, self.ratio), self.output_abs_path)


class HBondContactAnalysis(Analysis):
    # TODO: This is duplicated code that should be put into the super class ... if I understood how to call super in python.
    # Read essential python reference and writing idiomatic python ..
    def __init__(self, iso, ratio, abs_path):
        super(HBondContactAnalysis, self).__init__(iso, ratio)
        self.iso = iso
        self.ratio = ratio
        self.output_abs_path = abs_path
        logging.debug("Analysis results will be outputted to %s", self.output_abs_path) 
        if not os.path.exists(self.output_abs_path):
            logging.debug("%s does not exist. Creating ...", self.output_abs_path)
            os.makedirs(self.output_abs_path)

    def _command_by_ligand(self, xtc_filename, system_index, testing=False):
        test_command = '-b 0 -e 250000'
        if testing:
            test_command = '-b 0 -e 100'
 
        shell_command = string.Template("seq 1 ${ratio} | parallel -j 16 \"echo {} Protein | g_hbond -f $data_dir/$xtc_name -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}_by_ligand.ndx -nonitacc -nomerge -num $output_dir/ab_${iso}_${ratio}_${sys_index}_ins{} $testing > /dev/null 2>&1\"").substitute(data_dir=self.DATA_BASE_DIR, xtc_name=xtc_filename, iso=self.iso, ratio=self.ratio, output_dir='/dev/shm/grace', sys_index=system_index, testing=test_command)

        return shell_command

    def _command_by_residue(self, xtc_filename, system_index, testing=False):
        res_start=0
        res_end=134
        ligand_grp=135

        test_command = '-b 0 -e 250000'
        if testing:
            test_command = '-b 0 -e 100'
 
        shell_command = string.Template("seq $res_start $res_end | parallel -j 8 \"echo {} $ligand_grp | g_hbond -f $data_dir/$xtc_name -s ${iso}_${ratio}_nosol.tpr -n g_hbond_${ratio}_${iso}_by_residue.ndx -nonitacc -nomerge -num $output_dir/ab_${iso}_${ratio}_${sys_index}_residue{} $testing > /dev/null 2>&1\"").substitute(res_start=res_start, res_end=res_end, ligand_grp=ligand_grp, data_dir=self.DATA_BASE_DIR, xtc_name=xtc_filename, iso=self.iso, ratio=self.ratio, output_dir='/dev/shm/grace', sys_index=system_index, testing=test_command)

        return shell_command

    # Calculates the number of hydrogen bonds made with each residue.
    # g_bond is ran for each ligand molecule
    def analyze_contact_by_ligand(self, testing=False):
        logging.debug("Analyzing nonpolar contact for each inositol")
        
        if testing:
            logging.debug("Running command as a test")

        try:
            for sys_index in range(self.START, self.FINISH):
                xtc = str(sys_index) + "_final.xtc"
                xtc_full_path = os.path.join(self.DATA_BASE_DIR, xtc)
                
                if os.path.exists(xtc_full_path):
                    shell_command = self._command_by_ligand(xtc, sys_index, testing=testing)

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

        self._clean("hbonds_{0}_{1}_by_ligand".format(self.iso, self.ratio), self.output_abs_path)

    # Calculates the number of hydrogen bonds made with each residue.
    # g_hbond is ran each time for each residue-ligands in the protein
    def analyze_contact_by_residue(self, testing=False):
        logging.debug("Analyzing nonpolar contact for each inositol")
        
        if testing:
            logging.debug("Running command as a test")

        try:
            for sys_index in range(self.START, self.FINISH):
                xtc = str(sys_index) + "_final.xtc"
                xtc_full_path = os.path.join(self.DATA_BASE_DIR, xtc)
                
                if os.path.exists(xtc_full_path):
                    shell_command = self._command_by_residue(xtc, sys_index, testing=testing)

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

        self._clean("hbonds_{0}_{1}_by_residue".format(self.iso, self.ratio), self.output_abs_path)


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
    OUTPUT_BASE_DIR = os.environ['PWD']
    hbonds = HBondContactAnalysis(os.path.join(OUTPUT_BASE_DIR, "test/hbonds"))

    hbonds.analyze_contact_by_ligand(Analysis.ISO, Analysis.RATIO, Analysis.RATIO, testing=True)
    hbonds.analyze_contact_by_residue(Analysis.ISO, Analysis.RATIO, Analysis.RATIO, testing=True)

    nonpolar = NonpolarContactAnalysis(os.path.join(OUTPUT_BASE_DIR, "test/nonpolar"))
    nonpolar.analyze_contact_by_ligand(Analysis.ISO, Analysis.RATIO, Analysis.RATIO, testing=True)
   
