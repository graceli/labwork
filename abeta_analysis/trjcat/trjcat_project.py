import shlex
import subprocess
import logging
import os
import sys
import glob
import optparse

# TODO this is very similar in nature to Oliver Beckstein's GromacsWrapper
class GromacsCommand:
    def __init__(self, executable, **kwargs):
        self.gromacs_exe = executable
        self.gromacs_args = kwargs  

    def run(self):
        """ Executes the command as a subprocess """        
        command = self._build_command()
        args = shlex.split(command)
        process = subprocess.Popen(args, stdout=open(os.devnull)).communicate(self.gromacs_args["pipe"])
        logging.info("Executed %s", command)
        logging.info("Finished with %s", process)
    
    def _build_command(self):
        """ Returns the Gromacs command as a string """        
        files_str = self.gromacs_args["xtc"]
        temp_name = self.gromacs_args["output"]
        index_file = self.gromacs_args["index"]
            
        command = "%s -f %s -o %s -n %s" % (self.gromacs_exe, files_str, temp_name, index_file)
        
        if "custom" in self.gromacs_args:
            command = command + " " + self.gromacs_args["custom"]
        
        # TODO fix this. This is a hack. Handle differing arguments properly
        if "tpr" in self.gromacs_args:
            command = command + " " + self.gromacs_args["tpr"]
                        
        return command
   

class Trajectory:
    def __init__(self, name, project_path, traj_path, trajs):
        self.name = name 
        self.project_path = project_path
        self.traj_path = traj_path
        # location of the trajectories
        self._files_to_cat = trajs

    def build(self, tpr, index_file, index_group, project_output):
        files_str = " ".join([ os.path.join(self.traj_path, t) for t in self._files_to_cat])
        
        if files_str == "":
            print "Nothing to concatenate for", self.name
            return
            
        index_file = os.path.join(self.project_path, index_file)
        temp_outfile = os.path.join(project_output, self.name + "_temp")
        trjcat = GromacsCommand('trjcat', xtc=files_str, output=temp_outfile, index=index_file, pipe=index_group)
        trjcat.run()

        final_output = os.path.join(project_output, self.name)
        trjconv = GromacsCommand('trjconv', xtc=temp_outfile, tpr="-s " + os.path.join(project_output, tpr), output=final_output, index=index_file, custom='-pbc mol', pipe=index_group)
        trjconv.run()               
                
    def check(self):
        command = "gmxcheck -f %s" % (self.name)
        print command  

# Represents a simulation project generated using gromacs
# What do you have to have at the minimum to consider having some data?
class Project:
    def __init__(self, name, tpr, output, index="index.ndx", trajectories=[]):
        self.project_name = name
        self.tpr = tpr
        
        # TODO raise an error if the index file was not provided
        self.index_file = index
        self.trajectories = trajectories
        self.subdirectories = []
        self.project_output = output
        
    def build_trajectories(self, index_group):
        self._prepare_for_build()
        
        for i in range(len(self.trajectories)):
            self.trajectories[i].build(self.tpr, self.index_file, index_group, self.project_output)

    def data_info(self):
        # log a list of trajectories produced and their file sizes
        # Project: ProjectName
        # N directories trjcatted
        # Group selected: "Non-Solvent"
        # traj    size
        # 1       100 M
        # ...
        pass
        
    def add_trajectory(self, traj):
        self.trajectories.append(traj)

    def add_directory(self, subdir):
        self.subdirectories.append(subdir)

    def _prepare_for_build(self):
        logging.info("Preparing for trajectory building")
        try:
            os.mkdir(self.project_output)
        except OSError as e:
            print e
            logging.info("Exception: %s", str(e))
        
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

    # TODO error checking
    # TODO refactor to configuration file
    # component such as protein, solvent etc defined in the index file
    
    # system_component = "Protein"
    # project_name = "TestProject"
    # project_output = project_name + "_" + system_component
    # subdir_prefix = "sys"
    # N = 10     
    
    usage = "usage: %prog [options] project_name"
    parser = optparse.OptionParser(usage, description='Trjcat some trajectories')                                          
    
    # parser.add_option("-p", "--project_name", dest="project_name", help='Existing project name')
    parser.add_option("-o", "--project_output", dest="project_output", 
        help='New project directory', default="Test")
    parser.add_option("-f", "--subdir_prefix", dest="subdir_prefix", 
        help='Optional prefix for the project subdirectory', default="")
    parser.add_option("-N", "--num_replicas", type=int, dest="N", 
        help="Number of subdirectories (for testing purposes only)")
    parser.add_option("-n", "--component", dest="system_component",
        help="The component of the system to get out (protein, etc)", default="System")
        
    (options, args) = parser.parse_args()

    # TODO error handling for add_option -- look into how to properly do this 
    # http://docs.python.org/library/optparse.html
    if len(args) != 1:
          parser.error("Incorrect number of arguments")

    project_name = args[0]
    if not os.path.exists(project_name):
        print "project {0} does not exist".format(project_name)
        sys.exit(1)

    logging.info("Initializing trjcatting for project %s", project_name)

    p = Project(project_name, "TestProject.tpr", options.project_output, index="index.ndx")
        
    # Read the project directory and build a list of files to trjcat
    for dir_idx in range(options.N):
        project_subdir = options.subdir_prefix + str(dir_idx)
        p.add_directory(project_subdir)
        
        traj_path = os.path.join(project_name, project_subdir)
        
        results = list_xtcs(traj_path)
        print results
        
        p.add_trajectory(Trajectory(str(dir_idx) + "_final", project_name, traj_path, results))
 
    p.build_trajectories(options.system_component)
    

if __name__ == '__main__':
    main()


