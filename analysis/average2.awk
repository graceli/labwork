#!/usr/bin/awk -f

# columns of the data file
# time intra inter total_inos_pep inos1_inos-pep inos1_num_chains inos2_inos-pep inos2_num_chains ...
# works for arbitrary number of inositols
# must be called with -v nins=X, where X is the number of inositols

BEGIN{print ARGV[1];intra=0;inter=0;inos=0}
{
	if($1<=100000){
		intra+=$2;
		inter+=$3;
		inos+=$4;
		# take into consideration of more than 2 inositols
		for(i=5; i<=NF; i+=2){
			hbs[ $i,$(i+1) ]++;
		}
	}
}
END{
        user=ENVIRON["USER"];
        print "file analyzed :",user"@:"ENVIRON["PWD"]"/"FILENAME"\n";
	printf "%20s%12d\n", "# intra-peptide HBs = ",intra;
	printf "%20s%12d\n", "# inter-peptide HBs = ",inter;
	printf "%20s%12d\n", "# inositol-peptide HBs = ",inos;

	avg_intra=intra/100000;
	avg_inter=inter/100000;
	avg_inos=inos/100000;

	print "\n";

	print "intra/frame=",avg_intra;
	print "inter/frame=",avg_inter;
	print "inos/frame=",avg_inos;
	print "(intra+inter+inos)/frame=", avg_intra+avg_inter+avg_inos;

	print "\n";
	######################## prints number of inositols #########################
	print "Number of inositols";
	printf "%16d%12d%12d%12d%12d%c","0","1","2","3","4","\n";
	for(i=0;i<=12;i++){
		printf "%d ",i;
		for(j=0;j<=4;j++){
			if(hbs[i,j]){
				printf "%d ",hbs[i,j];
				hbs_sum[i]+=hbs[i,j];
				chains_sum[j]+=hbs[i,j];
			}else{
				printf "%s ","0";
			}
		}
		printf "%d\n",hbs_sum[i];
	}
	printf "%s ","0";
	for(j in chains_sum){
		printf "%d ",chains_sum[j];
	}
	print "\n";
        #############################################################################	
	#prints each entry as a percentage of total number of inositols (2*num_frames)
	print "% total number inositols in simulation";
        printf "%16d%12d%12d%12d%12d%c","0","1","2","3","4","\n";
        for(i=0;i<=12;i++){
                printf "%d ",i;
                for(j=0;j<=4;j++){
                        if(hbs[i,j]){
                                printf "%.4f ",hbs[i,j]*100/(nins*100000);
                        }else{
                                printf "%s ","0";
                        }
                }
                printf "\n";
        }
	print "\n";
        #############################################################################	
	#prints each entry as a percentage of total n number of inos-pep HB formed
	print "% of n number of inositol-peptide HBs formed";
        printf "%16d%12d%12d%12d%12d%c","0","1","2","3","4","\n";
        for(i=0;i<=12;i++){
                printf "%d ",i;
                for(j=0;j<=4;j++){
                        if(hbs[i,j]){
                                printf "%.4f ",hbs[i,j]*100/hbs_sum[i];
                        }else{
                                printf "%s ","0";
                        }
                }
                printf "\n";

        }
	print "\n";
	#############################################################################
        print "% of number of inositol chains HBs formed for m chains bound";
        printf "%16d%12d%12d%12d%12d%c","0","1","2","3","4","\n";
        for(i=0;i<=12;i++){
                printf "%d ",i;
                for(j=0;j<=4;j++){
                        if(hbs[i,j]){
                                printf "%.4f ",hbs[i,j]*100/chains_sum[j];
                        }else{
                                printf "%s ","0";
                        }
                }
		printf "\n";
        }
	print "\n";	
	#############################################################################
	print "Summary:\n";
	print "#inositol vs #HB" > "hbs.stats";	
	for(k=0; k<=12; k++){
		#PRINT
		#printf "%4d%12d%12.4f\n", k,hbs_sum[k],hbs_sum[k]*100/(nins*100000) > "hbs.stats";
		print k,hbs_sum[k], hbs_sum[k]*100/(nins*100000) > "hbs.stats";
		sum_hb_percent += hbs_sum[k]*100/(nins*100000);
	}	
	printf "%4s%24.4f\n\n","sum=",sum_hb_percent;

	print "#inositol vs #chains bound\n" > "chains.stats";
	for(l in chains_sum){
		#printf "%4d%12d%12.4f\n", l,chains_sum[l],chains_sum[l]*100/(nins*100000) > "chains.stats";
		print l,chains_sum[l],chains_sum[l]*100/(nins*100000) > "chains.stats";
		sum_chain_percent += chains_sum[l]*100/(nins*100000);
	}
	printf "%4s%24.4f\n","sum=",sum_chain_percent;
}
