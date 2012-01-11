import tarfile
import os
import glob
import re
import csv

# Return the xtc files in the tarfile as a set
def sanitize_name(name):
	pattern = re.compile(':')
	new = pattern.sub('_', name)
	
	return new
		
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
	basename = "lgw" + replica + "." + sequence
	return basename + "_small.xtc", basename + '.edr'

def construct_logname(directory):
	pattern = re.compile(':')
	new = pattern.sub('_', directory)
	base = os.path.basename(new)
	return base + '_database.dat.noforce.300K_clean.csv'

def parse_names_from_log(logname):
	cr = csv.reader(open(logname, 'rb'), delimiter=',')	
	fileslist = []
	for row in cr:
		xtcname, edrname = construct_fname(row)	
		fileslist.append(xtcname)
		fileslist.append(edrname)
	return fileslist

def main():
    # list and store a list of the data directories
	DATA = '/scratch/p/pomes/grace/drSH3'
	LOGS = '/scratch/p/pomes/grace/hpss/sh3/SH3_analysis/sh3_300_logs'
	CWD =  '/scratch/p/pomes/grace/sha/2012'

	data_dirs = glob.glob(os.path.join(DATA, "PRIOR*"))

	#data_dirs =[ os.path.join(DATA, 'PRIOR_TO_RESTART_Tue_Nov_2_17:35:28_EDT_2010') ]

	for d in data_dirs:
		# make a listing of the tar files in the run directory
		print "processing:", d
		data_path = os.path.join(d, "output/data")

		#run_dir_tarfiles = glob.glob(os.path.join(data_path, 'lgw*tmp*'))
		run_dir_tarfiles = [ os.path.join(data_path, 'lgw99.52166-lgw99.52749.tmp.458') ] #os.path.join(data_path, 'lgw96.11028-lgw101.11709.tmp.742') ]

		# load the corresponding log file
		logname = os.path.join(LOGS, construct_logname(d))
		
		print "scanning in log file:", logname

		files_at_300 = parse_names_from_log(logname)

		# extract from each tar file in the run dir
		# the xtc file at T=300K (beta=1.67)
		num_extracted = 0
		for tf in run_dir_tarfiles:
			print "searching in tarfile", tf
			all_xtc_in_tar = xtcfiles(tf)
			for file in files_at_300:
				if file in all_xtc_in_tar: 
					# extract file from tar file	
					tar_obj = tarfile.open(tf)
					tar_obj.extract(file, path="/dev/shm/grace")
					print "extracted", file
					num_extracted += 1
		print num_extracted, "files extracted" 
	
		new_tar = os.path.join('/dev/shm', sanitize_name(d) + '.tar.gz')
		ret_code = subprocess.Popen(['tar', 'xvfz', new_tar, '/dev/shm/grace/*.xtc'])
		if ret_code == 0:
			subprocess.Popen(['mv', new_tar,
			

if __name__ == '__main__':
    main()

