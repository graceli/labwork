import logging
import os
import sys
import glob

class GromacsCommand:
    # how do pass in multiple kwargs?
    def __init__(self, executable, args):
        self.gromacs_exe = executable
        self.gromacs_args = args
        
    def run():
        # run command
        # wrapper around 
        pass
    
    # does python have some sort of default to_string method for a class that I can override ?
    # like in ruby or Java?
    def to_string():
        # returns the command as a string
        # string version of what would be ran
        pass
        
class Trajectory:
    def __init__(self, name, dir_name, trajs):
        # name of the single trajectory
        self.name = name

        # names of the trajectory pieces
        self.files_to_cat = trajs

        # location of the trajectories
        self.dir_name = dir_name
    
    def build_traj(self, index_group):
        files_str = " ".join(self.files_to_cat)
        # input is index_group                                                                    
        index_file = os.path.join(self.dir_name, "index.ndx")
        temp_name = self.dir_name
        command = "trjcat -f %s -o %s -n %s" % (files_str, temp_name, index_file)
        args = shlex.split(command)
        p = subprocess.Popen(args, stdout=os.devnull).communicate(index_group)
        logging.info("%s executed", command)
        
        if p.returncode == None || p.returncode != 0:
            raise CallProcessError
        
        final_name = temp_name + "_processed"
        command = "trjconv -f %s -o %s -pbc mol" % (self.dir_name, final_name)
        
        logging.info("%s executed", command)
        
    def check_traj(self):
        command = "gmxcheck -f %s" % (self.name)
        print command

class Project:
    def __init__(self, trjs):
        self.trajectories = trjs 

    def data_info(self):
        # log a list of trajectories produced and their file sizes
        # Project: ProjectName
        # N directories trjcatted
        # Group selected: "Non-Solvent"
        # traj    size
        # 1       100 M
        # ...
        
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
    files_for_trjcat = []
    for dir_idx in range(N):
        trj_dir_path = os.path.join(project_name, str(dir_idx))
        results = []
        if os.path.exists(trj_dir_path):
            try:
                results = glob.glob("%(trj_dir_path)s/*.xtc" % vars())
                print results
            except OSError as e:
                results = []
                
            trjs = Trajectory(trj_dir_path, results, "Non-Solvent")
            files_for_trjcat.append(trjs)
        else:
            logging.info("Directory %s does not exist", dir_idx)
    

    for t in files_for_trjcat:
        t.build_traj()
        t.traj_info()


if __name__ == '__main__':
    main()


