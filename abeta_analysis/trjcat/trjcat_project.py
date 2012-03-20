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
        self._files_to_cat = trajs

    def build(self, tpr, index_group, project_output):
        files_str = " ".join([ os.path.join(self.traj_path, t) for t in self._files_to_cat])
        
        if files_str == "":
            print "Nothing to concatenate for", self.name
            return
            
        # input is index_group                                                                    
        index_file = os.path.join(self.project_path, "index.ndx")
        temp_name = self.name + "_temp"
        
        trjcat = GromacsCommand('trjcat', xtc=files_str, tpr=tpr, output=os.path.join(project_output, temp_name), index=index_file)
        trjcat.run()

        outfile = os.path.join(project_output, self.name)
        trjconv = GromacsCommand('trjconv', xtc=temp_name, tpr=tpr, output=outfile, index=index_file, custom='-pbc mol')
        trjconv.run()               
        
        # place holder for file creation
        os.system("touch {0}".format(outfile))
        
    def check(self):
        command = "gmxcheck -f %s" % (self.name)
        print command  

# Represents a simulation project generated using gromacs
# What do you have to have at the minimum to consider having some data?
class Project:
    def __init__(self, name, tpr, output, trajectories=[]):
        self.project_name = name
        self.tpr = tpr
        self.trajectories = trajectories
        self.subdirectories = []
        self.project_output = output
        
    def build_trajectories(self, index_group):
        self._prepare_for_build()
        
        for i in range(len(self.trajectories)):
            self.trajectories[i].build(self.tpr, index_group, self.project_output)

    def data_info(self):
        # log a list of trajectories produced and their file sizes
        # Project: ProjectName
        # N directories trjcatted
        # Group selected: "Non-Solvent"
        # traj    size
        # 1       100 M
        # ...
        pass
        
    def _prepare_for_build(self):
        logging.info("Preparing for trajectory building")
        # create the output directory
        try:
            os.mkdir(self.project_output)
        except OSError as e:
            print e
            logging.info("Exception: %s", str(e))
            
        
        # return True
        
    def add_trajectory(self, traj):
        self.trajectories.append(traj)

    def add_directory(self, subdir):
        self.subdirectories.append(subdir)
        
def list_xtcs(directory):
    results = []
    if os.path.exists(directory):
        results = glob.glob("%(directory)s/*.xtc" % vars())
        if results == []:
            logging.info("No trajectories were found")
    else:
        logging.info("Directory %s does not exist", directory)                               
    
    return results
    
def main():
    FORMAT = '%(asctime)s %(levelname)s %(message)s'
    logging.basicConfig(filename='trjcat_project.log', format=FORMAT, level=logging.DEBUG)
    logger = logging.getLogger('trjcat')

    # Usage: trjcat_project <project_name>

    # TODO error checking
    # TODO refactor to configuration file
    # component such as protein, solvent etc defined in the index file
    system_component = "Protein"
    project_name = "TestProject"
    project_output = project_name + "_" + system_component
    
    if not os.path.exists(project_name):
        print "project {0} does not exist".format(project_name)
        sys.exit(1)

    logging.info("Initializing trjcatting for project %s", project_name)
    N = 10     

    p = Project(project_name, "TestProject.tpr", project_output)
        
    # Read project directory and build a list of files to trjcat
    for dir_idx in range(N):
        # trj_dir_path = os.path.join(project_name, str(dir_idx))
        project_subdir = str(dir_idx)
        p.add_directory(project_subdir)
        
        traj_path = os.path.join(project_name, project_subdir)
        
        results = list_xtcs(traj_path)
        print results
        
        p.add_trajectory(Trajectory(str(dir_idx) + "_final", project_name, traj_path, results))
 
    p.build_trajectories(system_component)
    

if __name__ == '__main__':
    main()


