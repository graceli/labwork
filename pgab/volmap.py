import glob
import sys
import os

traj_list = []

output_temp = "/dev/shm/grace"
if not os.path.exists(output_temp):
	os.mkdir(output_temp)

output_disk = os.path.join(os.environ['PWD'], "volmap")
if not os.path.exists(output_disk):
	os.mkdir(output_disk)

#for dir_idx in range(37, 51):
#	traj_list.extend(glob.glob("pgab/%(dir_idx)s/*prod*.xtc" % vars()))

#traj_str = " ".join(traj_list)
#num_cs = len(traj_list)
#command = "`for i in `seq 1 %s`; do echo c; done` | trjcat -f %s -o %s/pgab_37-50_all_dt10.xtc -n pgab/pgab_volmap.ndx -dt 10 -settime -cat" % (str(num_cs), traj_str, output_temp)
#os.system(command)	
command = "cp pgab_37-50_all_dt10.xtc /dev/shm/grace"
os.system(command)	

command = "echo nonsolvent | trjconv -f %s/pgab_37-50_all_dt10.xtc -s pgab/pgab.tpr -o %s/pgab_37-50_all_dt10_mol.xtc -n pgab/pgab_volmap_nonsolvent.ndx -pbc mol" % (output_temp, output_disk)
os.system(command)

command = "echo Protein nonsolvent | trjconv -f %s/pgab_37-50_all_dt10_mol.xtc -s pgab/pgab.tpr -o %s/pgab_37-50_all_dt10_fit.xtc -n pgab/pgab_volmap_nonsolvent.ndx -fit rot+trans" % (output_disk, output_disk)
os.system(command)

