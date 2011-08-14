#!/usr/bin/env python
import glob
import os


fileslist = glob.glob("*")

for file in fileslist:
	if os.path.getsize(file) == 0:
		print file, "has size of zero. Deleting file."
		os.remove(file)


