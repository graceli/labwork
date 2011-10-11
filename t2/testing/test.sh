dd if=/dev/zero of=/3ware/c0u0/dd_sdb1.out bs=1048576 count=1024 2> dd_sdb1.log
echo "completed /dev/sdb1 count=1024"
dd if=/dev/zero of=/3ware/c0u0/dd_sdb1.out bs=1048576 count=4096 2>> dd_sdb1.log
echo "completed /dev/sdb1 count=4096"
dd if=/dev/zero of=/3ware/c0u0/dd_sdb1.out bs=1048576 count=40960 2>> dd_sdb1.log
echo "completed /dev/sdb1 count=40960"

dd if=/dev/zero of=/3ware/c0u1/dd_sdc1.out bs=1048576 count=1024 2> dd_sdc1.log 
echo "completed /dev/sdc1 count=1024"
dd if=/dev/zero of=/3ware/c0u1/dd_sdc1.out bs=1048576 count=4096 2>> dd_sdc1.log
echo "completed /dev/sdc1 count=4096"
dd if=/dev/zero of=/3ware/c0u1/dd_sdc1.out bs=1048576 count=40960 2>> dd_sdc1.log
echo "completed /dev/sdc1 count=40960"
