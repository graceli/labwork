#print lines that are bidentate; columns 2-13 are bidentate states
awk '{sum=$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13; if(sum > 0) print FNR,"\t",$0}'
