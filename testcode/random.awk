#!/usr/bin/awk -f 

function randint(n) {
     return int(n * rand())
}



{
	if(NR%2){
		printf("%s %d %d %d\n", $1, $2, $3-randint(100000), $4+randint(100000));
	}else{
		printf("%s %d %d %d\n", $1, $2, $3+randint(100000), $4-randint(100000));
	}
}

