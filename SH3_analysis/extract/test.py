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

def construct_fname(row):
	replica = row[0]
	sequence = row[1]
	return "lgw" + replica + "." + sequence + "_small.xtc"

def construct_logname(directory):
	pattern = re.compile(':')
	new = pattern.sub('_', directory)
	return new + '_database.dat.noforce.300K_clean.csv'

def parse_names_from_log(logname):
	cr = csv.reader(open(logname, 'rb'), delimiter=',')	
	fileslist = []
	for row in cr:
		filename = construct_fname(row)	
		fileslist.append(filename)
	return fileslist

def main():
    # list and store a list of the data directories
	DATA = '/scratch/p/pomes/grace/drSH3'
	data_dirs = glob.glob(os.path.join(DATA, "PRIOR*"))

	for d in data_dirs:
		# make a listing of the tar files in the run directory
		run_dir_tarfiles = glob.glob(os.join.path(d, 'lgw*tmp*'))
		# load the corresponding log file
		logname = os.path.join(LOGS, construct_logname(d))
		files_at_300 = parse_names_from_log(logname)

		# extract from each tar file in the run dir
		# the xtc file at T=300K (beta=1.67)
		for tarfile in run_dir_tarfiles:
			all_xtc_in_tar = xtcfiles(tarfile)
			for file in files_at_300:
				if file in all_xtc_in_tar: 
					# extract file from tar file	
					tar_obj.extract(file, path="/dev/shm/grace")

                
if __name__ == '__main__':
    main()
