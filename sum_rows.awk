#!/usr/bin/awk -f

{
	for(i=1;i<=NF;i++){
		total_row[i]+=$i;
	}

}
END{
	for(i=1;i<=NF;i++){
		printf "%d ", total_row[i];
	}
	printf "\n";
}
