#!/usr/bin/awk -f
BEGIN{
	# print C1,C2;
	for(i=0;i<=24;i++){
		for(j=0;j<=24;j++){
			#printf "%d %d %.4f\n",i,j,matrix[i,j]*100/NR;
			matrix[i,j]=0;
		}
	}
}
{
	matrix[$1,$2]++;	
}
END{
	for(i=0;i<=5;i++){
		for(j=0;j<=5;j++){
			if(matrix[i,j]) {
				printf "%d %d %.4f\n",i,j,matrix[i,j];
			}else{
				printf "%d %d %.4f\n",i,j,0;
			}
		}
		printf "\n";
	}
}
