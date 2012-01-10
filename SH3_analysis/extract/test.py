import tarfile
import os
import glob

# Return the xtc files in the tarfile as a set
def xtcfiles(tarfilename):
    tar = tarfile.open(tarfilename)
    members = tar.getmembers()
    files = set([])
    for info in members:
        if os.path.splitext(info.name)[1] == '.xtc':
            files.add(info.name)
    return files

def main():
    # list and store a list of the data directories
    data_dirs = glob.glob(os.path.join(DATA, "PRIOR*")):
        
    print data_dirs
    
    # for d in data
    # LOGS="sh3_300_logs"
    #     logfiles = glob.glob(LOGS + "/*.csv")
    #     print logfiles 
    # for tarfile in run_dir:
    #     # build an in memory dictionary of tarfiles        
    #     files_set = xtcfiles(tarfile)
    #     
    #     # load corresponding log file
    #     for each file at T=300:
    #         if file is in tarfile
    #             extract file from tarfile
                
if __name__ == '__main__':
    main()
    # files = xtcfiles('lgw124.90188-lgw35.91189.tmp.455')
    # print len(files)
    # name = 'lgw124.910660000000000000_small.xtc'
    # if name in files:
    #     print name, "is in tar file"