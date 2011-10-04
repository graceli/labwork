#!/bin/sh                                                                      
function submit_inositol {                                                      
                                                                                
                                        
for iso in glycerol; do
        for r in 15 64; do
                qsub -v ANALYSIS="rmsd",ISO="$iso",RATIO="$r"
/home/grace/AnalysisScripts/abeta_analysis/analysis.sh                          
                                                           
                #qsub -v ANALYSIS="rmsf_calpha",ISO="$iso",RATIO="$r"
/home/grace/AnalysisScripts/abeta_analysis/analysis.sh
                #qsub -v ANALYSIS="hbonds",ISO="$iso",RATIO="$r"
/home/grace/AnalysisScripts/abeta_analysis/analysis.sh
                #qsub -v ANALYSIS="nonpolar",ISO="$iso",RATIO="$r"
/home/grace/AnalysisScripts/abeta_analysis/analysis.sh
                #qsub -v ANALYSIS="chain_hbonds",ISO="$iso",RATIO="$r"
/home/grace/AnalysisScripts/abeta_analysis/analysis.sh
        done;
done
}

# function submit_water {
# for iso in water; do 
#         for r in 64; do 
#                 qsub -v ANALYSIS="rmsd",ISO="$iso",RATIO="$r"
# /home/grace/AnalysisScripts/abeta_analysis/analysis.sh -N analysis_w
#                 #qsub -v ANALYSIS="rmsf_calpha",ISO="$iso",RATIO="$r"
# /home/grace/AnalysisScripts/abeta_analysis/analysis.sh -N analysis_w
#                 #qsub -v ANALYSIS="chain_hbonds",ISO="$iso",RATIO="$r"
# /home/grace/AnalysisScripts/abeta_analysis/analysis.sh -N analysis_w
#                 #qsub -v ANALYSIS="hbonds",ISO="$iso",RATIO="$r"
# /home/grace/AnalysisScripts/abeta_analysis/analysis.sh
#                 #qsub -v ANALYSIS="nonpolar",ISO="$iso",RATIO="$r"
# /home/grace/AnalysisScripts/abeta_analysis/analysis.sh
#         done;
# done
# }

submit_inositol
#submit_water
