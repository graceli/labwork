import shlex
import subprocess
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
    def __init__(self, name, project_path, traj_path, trajs):
        # names of the trajectory pieces
        self.name = name 
        self.project_path = project_path
        self.traj_path = traj_path

        # location of the trajectories
        self.files_to_cat = trajs

    def build(self, index_group):
        files_str = " ".join([ os.path.join(self.traj_path, t) for t in self.files_to_cat])
        # input is index_group                                                                    
        index_file = os.path.join(self.project_path, "index.ndx")
        temp_name = self.name + "_temp"
        command = "trjcat -f %s -o %s -n %s" % (files_str, temp_name, index_file)
        args = shlex.split(command)
        p = subprocess.Popen(args, stdout=open(os.devnull)).communicate(index_group)
        logging.info("%s executed", command)

        if p.returncode == None or p.returncode != 0:
            raise CallProcessError

        final_name = self.name

        command = "trjconv -f %s -o %s -pbc mol" % (self.dir_name, final_name)
        args = shlex.split(command)
        p = subprocess.Popen(args, stdout=os.devnull).communicate(index_group)
        logging.info("%s executed", command)

        if p.returncode == None or p.returncode != 0:
            raise CallProcessError  

    def check(self):
        command = "gmxcheck -f %s" % (self.name)
        print command  

class Project:
    def __init__(self, name, trjs, group):
        self.project_name = name
        self.trajectories = trjs
        self.contents = group

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

    # project = Project(project_name, trajectories, "Non-Solvent")
    for t in trajectories:
        t.build("Protein")    


if __name__ == '__main__':
    main()


