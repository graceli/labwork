import shlex
import subprocess
import logging
import os
import sys
import glob
import optparse
import datetime

# TODO this is very similar in nature to Oliver Beckstein's GromacsWrapper
class GromacsCommand:
    def __init__(self, executable, **kwargs):
        self.gromacs_exe = executable
        self.gromacs_args = kwargs  

    def run(self):
        """ Executes the command as a subprocess """        
        command = self._build_command()
        args = shlex.split(command)
        logging.info("Executing %s with input %s", command, self.gromacs_args["pipe"])
        process = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=open(os.devnull)).communicate(self.gromacs_args["pipe"])
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

    def build(self, tpr, index_file, index_group, project_output, temp_output="/dev/shm/grace", center=False):
        # check if the trajectory exists or not before building
        wildcard = os.path.join(project_output, "{0}*".format(self.name))

        if len(glob.glob(wildcard)) != 0:
            print self.name, ".xtc already exists on disk ... skipping trjcat ..."
            logging.info("%s.xtc already exists on disk ... skipping trjcat ...", self.name)
            return

        print self._files_to_cat
        
        num_trajs = len(self._files_to_cat)

    	files_str = ""
        print num_trajs

    	if num_trajs == 0 or num_trajs == 1:
    	    print "Nothing to concatenate for", self.name
    	    return
    	else:
    	    files_str = " ".join(self._files_to_cat)
            logging.debug("%s to be trjcatted", files_str) 
            index_file = os.path.join(self.project_path, index_file)
            temp_outfile = os.path.join(temp_output, self.name + "_temp")

            trjcat = GromacsCommand('trjcat', xtc=files_str, output=temp_outfile, index=index_file, pipe=index_group)
            trjcat.run()

            final_output = os.path.join(project_output, self.name)
            custom_command = "-pbc whole"
            pipe_command = index_group
            if center:
                custom_command = "-pbc res -center"
                pipe_command = "{0} {1}".format("center_group", index_group)
            
            trjconv = GromacsCommand('trjconv', xtc=temp_outfile, tpr="-s " + os.path.join(self.project_path, tpr), output=final_output, index=index_file, custom=custom_command, pipe=pipe_command)
            trjconv.run()
        
            # Remove temp files to avoid overflow if writing to /dev/shm
            # Bit of a hack fix
            # Never got this to work ... Silly
            # os.system("rm -f %(temp_outfile)s" % vars())        

    def check(self):
        command = "gmxcheck -f %s" % (self.name)
        print command  

    def __repr__(self):
        representation = " ".join(["Trajectory:", "name:", self.name, "path:", self.traj_path ])
        return representation

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
        
    def build_trajectories(self, index_group, temp="/dev/shm/grace", center=False):
        self._prepare_for_build()
        
        for i in range(len(self.trajectories)):
            self.trajectories[i].build(self.tpr, self.index_file, index_group, self.project_output, temp_output=temp, center=center)

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
        results = glob.glob("%(directory)s/*prod*.xtc" % vars())
        if results is None:
            logging.info("No trajectories were found")
    else:
        logging.info("Directory %s does not exist", directory)                               
    
    return results
    
def main():
    # TODO refactor to configuration file    
    usage = "usage: %prog [options] name start_idx end_idx"
    parser = optparse.OptionParser(usage, description='Trjcat some trajectories')                                          

    parser.add_option("-o", "--project_output", dest="project_output", 
        help='New project directory', default="Test")
    parser.add_option("-e", "--temp_dir", dest="temp_dir",
        help="Temp directory (default=/dev/shm/grace)", default="/dev/shm/grace")
    parser.add_option("-f", "--subdir_prefix", dest="subdir_prefix",
        help='Optional prefix for the project subdirectory', default="")
        
    # Note that if this option is set and a center_group index group is 
    # not defined in index file, then trjconv will fail
    parser.add_option("-c", "--center_system", dest="center_system", 
        help="Use a centering group and output by -pbc res mol", default=False)
    parser.add_option("-n", "--component", dest="system_component",
        help="The component of the system (in Gromacs index group language) to extract",        default="System")

    (options, args) = parser.parse_args()
    
    # TODO error handling for add_option -- look into how to properly do this 
    # http://docs.python.org/library/optparse.html
    if len(args) != 3:
          parser.error("Incorrect number of arguments")


    project_name = args[0]
    start_idx = int(args[1])
    end_idx = int(args[2])

    if not os.path.exists(project_name):
        print "project {0} does not exist".format(project_name)
        sys.exit(1)

    # Setup logging
    FORMAT = '%(asctime)s %(levelname)s %(message)s'
    now = datetime.datetime.now()
    logging.basicConfig(filename='trjcat_project_' + project_name + '_' + now.strftime("%Y-%m-%d-%H-%M") + '.log', format=FORMAT, level=logging.DEBUG)

    logging.info("Initializing trjcatting for project")
    logging.info("Ran with options=%s and args=%s", options, args)

    p = Project(project_name, project_name + ".tpr", options.project_output, index=project_name + ".ndx")

    # Read the project directory and build a list of files to trjcat
    all_dirs_failed = True
    
    # loop through project directories in the interval [start_idx, end_idx]
    for dir_idx in range(start_idx, end_idx + 1):
        project_subdir = options.subdir_prefix + str(dir_idx)
        p.add_directory(project_subdir)
        
        traj_path = os.path.join(project_name, project_subdir)
        
        results = list_xtcs(traj_path)
        if results is not None:
            all_dirs_failed = False
            
        traj = Trajectory(str(dir_idx) + "_final", project_name, traj_path, results)
        logging.info("Added %s", traj)
        p.add_trajectory(traj)
 
    if all_dirs_failed:
        print "No numerical subdirectories were found ... perhaps you are missing a prefix?"
        sys.exit(1)
    
    p.build_trajectories(options.system_component, center=options.center_system)
    

if __name__ == '__main__':
    main()


