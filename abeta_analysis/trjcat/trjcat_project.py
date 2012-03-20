import shlex
import subprocess
import logging
import os
import sys
import glob

# TODO this is very similar in nature to Oliver Beckstein's GromacsWrapper
class GromacsCommand:
    # how do pass in multiple kwargs?
    def __init__(self, executable, **kwargs):
        self.gromacs_exe = executable
        self.gromacs_args = kwargs  

    def run(self):
        """ Executes the command as a subprocess """
        
        command = self._build_command()
        # args = shlex.split(command)
        # p = subprocess.Popen(args, stdout=open(os.devnull)).communicate(index_group)
        print "executed", command
        logging.info("executed %s", command)
        
        # if p.returncode == None or p.returncode != 0:
        #     logging.warn("FAILED with", p.returncode)
        #     raise CallProcessError    
    
    def _build_command(self):
        """ Returns the Gromacs command as a string """
        
        files_str = self.gromacs_args["xtc"]
        tpr = self.gromacs_args["tpr"]
        temp_name = self.gromacs_args["output"]
        index_file = self.gromacs_args["index"]
            
        command = "%s -f %s -s %s -o %s -n %s" % (self.gromacs_exe, files_str, tpr, temp_name, index_file)
        
        if "custom" in self.gromacs_args:
            command = command + " " + self.gromacs_args["custom"]
            
        return command
   

class Trajectory:
    def __init__(self, name, project_path, traj_path, trajs):
        # names of the trajectory pieces
        self.name = name 
        self.project_path = project_path
        self.traj_path = traj_path

        # location of the trajectories
        self.files_to_cat = trajs

    def build(self, tpr, index_group):
        files_str = " ".join([ os.path.join(self.traj_path, t) for t in self.files_to_cat])

        # input is index_group                                                                    
        index_file = os.path.join(self.project_path, "index.ndx")
        temp_name = self.name + "_temp"
        
        trjcat = GromacsCommand('trjcat', xtc=files_str, tpr=tpr, output=temp_name, index=index_file)
        trjcat.run()

        final_name = self.name
        trjconv = GromacsCommand('trjconv', xtc=temp_name, tpr=tpr, output=final_name, index=index_file, custom='-pbc mol')
        trjconv.run()
        
    def check(self):
        command = "gmxcheck -f %s" % (self.name)
        print command  

# Represents a simulation project generated using gromacs
# What do you have to have at the minimum to consider having some data?
class Project:
    def __init__(self, name, tpr, traj_set):
        self.project_name = name
        self.tpr = tpr
        self.trajectories = traj_set

    def data_info(self):
        # log a list of trajectories produced and their file sizes
        # Project: ProjectName
        # N directories trjcatted
        # Group selected: "Non-Solvent"
        # traj    size
        # 1       100 M
        # ...
        pass

def main():
    FORMAT = '%(asctime)s %(levelname)s %(message)s'
    logging.basicConfig(filename='trjcat_project.log', format=FORMAT, level=logging.DEBUG)
    logger = logging.getLogger('trjcat')

    # Usage: trjcat_project <project_name>

    # TODO error checking
    # TODO refactor to configuration file
    project_name = "TestProject"
    
    if not os.path.exists(project_name):
        print "project {0} does not exist".format(project_name)
        sys.exit(1)

    logging.info("Initializing trjcatting for project %s", project_name)
    N = 10     

    # Read project directory and build a list of files to trjcat
    trajectories = []
    for dir_idx in range(N):
        trj_dir_path = os.path.join(project_name, str(dir_idx))
        results = []
        if os.path.exists(trj_dir_path):
            try:
                results = glob.glob("%(trj_dir_path)s/*.xtc" % vars())
                print results
            except OSError as e:
                results = []
                logging.log("Exception:", e)    

            trjs = Trajectory(str(dir_idx) + "_final", project_name, trj_dir_path, results)
            trajectories.append(trjs)
        else:
            logging.info("Directory %s does not exist", dir_idx)                               

    p = Project(project_name, "TestProject.tpr", trajectories)
    print p
    
    for t in trajectories:
            t.build(p.tpr, "Protein")    
         

if __name__ == '__main__':
    main()


