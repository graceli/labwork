#!/usr/bin/awk -f 

#inositol histogram is a 2D matrix where each column index corresponds to a inositol molecule
BEGIN{
	for(i=2; i<=NINS+1; i++){
		inositol_histogram[i]=0;
	}
}
{
	for(i=2;i<=NINS+1; i++){
		if($i > 0){
			#if there is one contact, count as bound
			inositol_histogram[i]++;	
		}
	}
}
END{
	print "#fraction bound unbound";
	for(i=2;i<=NINS+1;i++){
		fract_bound = inositol_histogram[i]/NR;
		fract_unbound = 1-fract_bound;
		printf "%d %f %f ", i-1,fract_bound, fract_unbound;
		if(fract_bound > 0.60) {
			print "> 0.6 bound";
		}else{
			printf "\n";
		}
			
	}
}
