import os
import subprocess
from collections import defaultdict
import sys

def printtabs(ntabs):
	for i in range(0,ntabs-1):
		print "\t",

def show(formFlag,showFlag,project):
	if showFlag == "all":
		for key in project.keys():
			ntabs=showSubproject(formFlag,key)
			if(formFlag == "expanded"):
				printtabs(ntabs)
			else:
				print len(project[key]), "number of files to be deleted"
	else:
		for key in project.keys():
			index=key.find(showFlag)
			if index>0:
				#print "found", showFlag, "in", key, index
				ntabs=showSubproject(formFlag,key)
				if(formFlag == "expanded"):
					printtabs(ntabs)
				else:
	                        	print len(project[key]), "files to be deleted"
		
	
def removeProjectFiles(project, file):
	code = 0
	if os.path.exists(file):
		print "deleting files in project", project, file
        	code = subprocess.check_call(["rm", file])

	return code

def touchFile(project,file):
	code=0
	print "touched file in project", project, file	
	code = subprocess.check_call(["touch", file])


def showSubproject(formFlag, path):
	parts=path.split(os.sep)
	if formFlag == "collapsed":
		print path,
	elif formFlag =="expanded":
		for i in range(0,len(parts)):
			printtabs(i)
			if i != 0:
				print "|----->",
			print parts[i]
	return len(parts)+1


if len(sys.argv)<2:
	print "usage: all or <project dir name> <collapsed|expanded>"
	sys.exit()

showFlag = sys.argv[1]
formFlag = sys.argv[2]

f=open("11086__grace_______pomes__________1TB____58438files")

project = defaultdict(list)
size=0
for line in f:
	columns = line.split()
	fileTodelete = columns[12]
	
	if showFlag != "all":
		if line.find(showFlag)>0:
			size += int(columns[10])

	path, filename = os.path.split(fileTodelete)
	#parts = path.split(os.sep)
	project[path].append(filename)

show(formFlag,showFlag,project)
print "total removal size=", float(size)/(1024*1024*1024), "gigs"
