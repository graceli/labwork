# prompts
# try with synchronization and using 8 cores to run stress tests

#nohup bonnie++ -d /3ware/c0u0 -s 32768 -n 1 -m "t2:/3ware/c0u0" -r 16000 -yp -u root >> u0_large.csv &
#nohup bonnie++ -d /3ware/c0u1 -s 32768 -n 1 -m "t2:/3ware/c0u1" -r 16000 -yp -u root >> u1_large.csv &

#nohup bonnie++ -d /3ware/c0u0 -m "t2:/3ware/c0u0" -yp -u root -q >> unit0.csv &
#nohup bonnie++ -d /3ware/c0u1 -m "t2:/3ware/c0u1" -yp -u root -q >> unit1.csv &

NFILES=1024

bonnie -p4 -u root
for i in `seq 1 2`; do
	nohup bonnie++ -d /3ware/c0u0 -s 4000 -n $NFILES -m "t2:/3ware/c0u0-${i}" -r 2000 -ys -u root >> u0_sync_${i}_${NFILES}.csv &
	nohup bonnie++ -d /3ware/c0u1 -s 4000 -n $NFILES -m "t2:/3ware/c0u1-${i}" -r 2000 -ys -u root >> u1_sync_${i}_${NFILES}.csv &
done
