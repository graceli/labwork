#!/bin/sh

for log in `ls */*.300K`; do
    # replace the pesky colons with an underscore
    string=$log
    new_log_name=${string//:/_}
    base=`basename $new_log_name`
    echo "cleaned $log"
    
    # transform log records from the form
    # record: replica#: 0   sequence#: 91067   w: 0.998446   w_nominal: 36
    # into a csv file like so:
    # #replica_num,sequence#,w,w_nominal
    # 0,91067,0.998446,36
    cat $log | awk 'BEGIN{print "#replica_num,sequence_num,w,w_nominal" }{print $3","$5","$7","$9 }' > ${base}_clean.csv
    
    # rename the original log file with the colons from the filename removed
    mv $log ${new_log_name}
done    
