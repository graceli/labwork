import logging
import os
import sys
import glob

class Trajectory: 
    # for each directory i
        # files = read_xtc_files(directory)
        # trjcat -f " ".join(files) -n project_directory/index_file -o tempfile_i
        # trjconv -f tempfile_i -pbc mol -o final_i
    
    def __init__(self, dir_idx, files, group):
        self.dir_name = dir_idx
        self.files_to_cat = files
        self.system_group = group
    
    def trjcat(self):
        files_str = " ".join(self.files_to_cat)
        # input is group                                                                    
        idx_file_path = os.path.join(self.dir_name, "index.ndx")
        command = "trjcat -f %s -o %s -n %s" % (files_str, self.dir_name, idx_file_path)
        # subprocess.call
        print command

def main():
    FORMAT = '%(asctime)s %(levelname)s %(message)s'
    logging.basicConfig(filename='trjcat_project.log', format=FORMAT, level=logging.DEBUG)
    logger = logging.getLogger('trjcat')
    
    # Usage: trjcat_project <project_name>

    # TODO error checking
    # temporarily use a fixed project name
    project_name = "TestProject"
    
    if not os.path.exists(project_name):
        print "project {0} does not exist".format(project_name)
        sys.exit(1)
    
    logging.info("Initializing trjcatting for project %s", project_name)
    N = 10

    # read project and build a list of files to trjcat
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
        t.trjcat()
    
    # log a list of trajectories produced and their file sizes
    # Project: ProjectName
    # N directories trjcatted
    # Group selected: "Non-Solvent"
    # traj    size
    # 1       100 M
    # ...

if __name__ == '__main__':
    main()


