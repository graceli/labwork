#!/usr/bin/perl

#this script generates a  with "n" number of strands with gromacs genconf
#and a topology file for the sheet
#Currently only works for gagagaga
#to implement: automatic resizing of boxes
#length of x-dimension=distance between strand ends images + length of strand
#length of z-dimension= distance between sheet layers) + thickness of sheet

if(@ARGV < 1){
	print STDERR "usage: perl generate_sheet.pl <num-strands> [-apsheet] [-ap]\n";
	print STDERR "note: specify distance in nanometers, -ap to generate antiparallel stacking\n";
	print STDERR "default: parallel sheet, parallel stacking\n";
	exit;
}

print "You have specified:", "nstrands=$ARGV[0] ", "sheet-dist = $ARGV[1] ", "strand-dist = $ARGV[2] ","\n";

$nstrands=$ARGV[0];
$sheet_opt=$ARGV[1];
$stack_opt=$ARGV[2];

$temp="GEN_SHEET_TEMP";

$sheet="parasheet";
#grow fibril
if($sheet_opt eq "-psheet"){
	system("genconf -f gagagaga_remodelled_final.gro -nbox 1 $nstrands 1 -dist 0.1 -o ${temp}_gagagaga_${sheet}_${nstrands}.gro > /dev/null 2>&1");
	#removes the velocities outputted by genconf 
	system("editconf -f ${temp}_gagagaga_${sheet}_${nstrands}.gro -o ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro -c");
}elsif($sheet_opt eq "-apsheet"){
	$sheet="apsheet";
	$nstrands_apsheet=$nstrands/2;
	system("genconf -f gagagaga_anti_onefaced.gro -nbox 1 $nstrands_apsheet 1 -dist 0.05 -o ${temp}_gagagaga_${sheet}_${nstrands}.gro > /dev/null 2>&1");
	system("editconf -f ${temp}_gagagaga_${sheet}_${nstrands}.gro -o ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro -c");
}
$new_nstrands = 2*$nstrands;


#stack
if($stack_opt eq "-pstack"){
	#stack parallel-ly (if -dist is 0.0 then packing b/n sheets is too gappy; specifying -0.1 nm allows the boxes of the images to stack to overlap a bit
	# which results in near perfect packing (visually verified)
	system("genconf -f ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro -nbox 1 1 2 -dist -0.1 -o gagagaga_${sheet}_parastack_${new_nstrands}.gro > /dev/null 2>&1");

}elsif($stack_opt eq "-apstack"){
	#stack antiparallelly by rotating about z by |theta| = 180
	#get number atoms (head -2 gets first 2 lines and tail -1 gets the last line)
	$num_atoms=`head -2 ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro | tail -1`;
	#print STDERR "$num_atoms\n";
	#get box size
	$box=`tail -1 ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro`;
	#print $box;
	@box_dims=split(/\s+/,$box);

	#length of y dimension or the "width of the sheet (fibril axis) is fixed (if this is changed, the sheet will not be infinite)
	$y_fixed = $box_dims[2];
	$x_fixed = $box_dims[1];
	$z = 2*$box_dims[3]; # double the z-dimension to accomodate the new sheet
	#rotate and center in box
	system("editconf -f ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro -rotate 0 0 180 -o ${temp}_gagagaga_${sheet}_${nstrands}_rz180.gro -c > /dev/null 2>&1"); 
	#print "gagagaga_sheet_${nstrands}_rz180.gro originally $x_fixed $y_fixed $box_dims[3]; resize to $x_fixed $y_fixed $z\n";
	#resize box to new box size
	system("editconf -f ${temp}_gagagaga_${sheet}_${nstrands}_rz180.gro -o ${temp}_gagagaga_${sheet}_${nstrands}_rz180_resized.gro -box $x_fixed $y_fixed $z -c > /dev/null 2>&1"); 
	# XY plane of the sheet being stacked; needs to be shifted to fit top sheet into surface grooves defined by the bottom sheet if sheet parallel
	# if sheet is antiparallel, do not shift by 0.2 nm
	$top_sheet_z = $box_dims[3]/2-0.1;
	if($sheet_opt eq "-apsheet"){
		system("editconf -f ${temp}_gagagaga_${sheet}_${nstrands}_rz180_resized.gro -translate 0 0 $top_sheet_z -o ${temp}_gagagaga_${sheet}_${nstrands}_rz180_tz.gro > /dev/null 2>&1"); 
	}else{
		system("editconf -f ${temp}_gagagaga_${sheet}_${nstrands}_rz180_resized.gro -translate 0.2 0 $top_sheet_z -o ${temp}_gagagaga_${sheet}_${nstrands}_rz180_tz.gro > /dev/null 2>&1"); 
	}
	#combine the two gro files -- is there a better way of doing this?
	#double number of atoms
	$new_num_atoms = 2*$num_atoms;
	chomp($date = `date`);
	system("echo \"combined stacked ${sheet} generated on $date \"> ${temp}_cutout1.gro");
	system("echo \'  $new_num_atoms\' >> ${temp}_cutout1.gro");
	#cut out the coordinate part of the bottom ${sheet} and paste into cutout1.gro
	system("tail -n +3 ${temp}_gagagaga_${sheet}_${nstrands}_novel.gro | head -n -1 >> ${temp}_cutout1.gro");
	#remove the first two lines of the top ${sheet} and paste into cutout2.gro
	system("tail -n +3 ${temp}_gagagaga_${sheet}_${nstrands}_rz180_tz.gro > ${temp}_cutout2.gro");
	#combine into one file
	system("cat ${temp}_cutout1.gro  ${temp}_cutout2.gro > ${temp}_gagagaga_ap_stacked_${new_nstrands}.gro");
	system("editconf -f ${temp}_gagagaga_ap_stacked_${new_nstrands}.gro -o gagagaga_${sheet}_apstack_${new_nstrands}.gro -c");
}

system("rm $temp\*");




