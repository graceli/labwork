import sys
import os
import re

def delete(file):
	file_nonewline = file.rstrip()
	base = os.path.basename(file_nonewline)
	filename,ext = os.path.splitext(base)

	if len(filename) > 0:
		firstChar = filename[0]
		#if firstChar == '#':
		#	return file

		if re.match(".e[0-9]+", ext) or re.match(".o[0-9]+",ext):
			return file.rstrip()

def main():
	"""docstring for main"""
	
	# grab a list of files from STDIN
	count = 0
	fileslist = sys.stdin.readlines()

	todelete = filter(delete, fileslist)	
	for file in todelete:
		print "deleting", repr(file.rstrip())
		os.system("rm %(file)s" % vars())
	
if __name__ == '__main__':
	main()

