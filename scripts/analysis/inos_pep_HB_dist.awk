#computes the distribution of inositol-peptide HBs and outputs to a file
{
	dist[ $4 ] ++;	
}
END{
	for(i=0; i<=12; i++  ){
		if(i in dist){
			printf "%d %f\n",i,dist[i]/NR > "inos_pep_dist_pd.hist";
			printf "%d %d\n",i,dist[i] > "inos_pep_dist.hist";
		}else{
			printf "%d %d\n",i,0 > "inos_pep_dist_pd.hist";
			printf "%d %d\n",i,0 > "inos_pep_dist.hist";
		}	
	}
}
