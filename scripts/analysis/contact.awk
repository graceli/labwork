#!/usr/bin/awk -f

{
	if($1<=100000){
		contact_inter_inos[$3,$4]++;
		contact_intra_inos[$2,$4]++;
		contact_inter_intra[$3,$2]++;
	}
}
END{
		
#	for(j=0;j<=12; j++){
#		printf "%12d", j;
#	}
#	printf "\n";
	for(i=0; i<=24; i++){
#		printf "%d", i;
		for(j=0;j<=12; j++){
			printf "%d %d %.4f\n",i,j,contact_inter_inos[i,j]*100/100000 > "contact_inter_inos.map";
			printf "%d %d %.4f\n",i,j,contact_inter_intra[i,j]*100/100000 > "contact_inter_intra.map";
		}
		printf "\n" > "contact_inter_inos.map";
		printf "\n" > "contact_inter_intra.map";
	}

	for(i=0;i<=12; i++){
		for(j=0;j<=12;j++){
			printf "%d %d %.4f\n",i,j,contact_intra_inos[i,j]*100/100000 > "contact_intra_inos.map";
		}
		printf "\n" > "contact_intra_inos.map";
	}
}
