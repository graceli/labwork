import fileinput
import os
import time

for line in fileinput.input():
	new_path = os.path.join(os.environ['SCRATCH'], line[15:]).rstrip()
	if os.path.exists(new_path):
		file_atime = os.path.getatime(new_path)
		start_of_today = time.time()
		if file_atime < start_of_today:	
			# update time of file
			os.utime(new_path, None)	
			print "processed", new_path, "at new access time", time.ctime(file_atime)
		else:
			print new_path, "not updated", "last accessed at", time.ctime(file_atime)
	else:
		print "file", new_path, "not found"


