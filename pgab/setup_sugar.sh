
function make_box {
    editconf -f $GRO -box 3 3 3 -o nag_box.gro
    genbox -cp nag_box.gro -cs ~/gro/tip3.gro -p nag_glycam.top
}

make_box glycam_3_12_2012_752_GMX.gro