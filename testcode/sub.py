from subprocess import Popen, PIPE

(stdout, stderr) = Popen(["cat","foo.txt"], stdout=PIPE).communicate()
print stdout
print stderr
