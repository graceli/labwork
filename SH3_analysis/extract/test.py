import tarfile
import os

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
    # build an in memory dictionary of tarfiles
    tarfiles_dict = {}
    for each tar file:
        tarfiles_dict[tarfile] = xtcfiles(tarfile)
        
        for each file at T=300:
            if file is in tarfile
                extract file from tarfile
                
if __name__ == '__main__':
    files = xtcfiles('lgw124.90188-lgw35.91189.tmp.455')
    print len(files)
    name = 'lgw124.910660000000000000_small.xtc'
    if name in files:
        print name, "is in tar file"