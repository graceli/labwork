# March 25 2012
# This is yet another dumb script for a data migration from colosse to scinet
# When I migrated Abeta project to colosse to extend my runs, I renamed "0" to "10" as the array on colosse did not have "0"

# move the files from directory 10 to 0

set -x

function move_15_directories {
    for iso in scyllo chiro glycerol water; do
        cd $iso
        mv 10/sys10_prod.part0019.xtc 0/sys0_prod.part0019.xtc
        mv 10/sys10_prod.part0019.log 0/sys0_prod.part0019.log
        mv 10/sys10_prod.part0019.edr 0/sys0_prod.part0019.edr
        mv 10/sys10_prod.cpt 0/sys0_prod.cpt
        mv 10/sys10_prod_prev.cpt 0/sys0_prod_prev.cpt
        cd ../
    done
}

function move_64_directories {
    for iso in scyllo chiro glycerol water; do
        cd $iso
        mv 10/sys10_prod.part0015.xtc 0/sys0_prod.part0015.xtc
        mv 10/sys10_prod.part0015.log 0/sys0_prod.part0015.log
        mv 10/sys10_prod.part0015.edr 0/sys0_prod.part0015.edr
        
        if [ -e 10/sys10_prod.part0016.xtc ]; then
            mv 10/sys10_prod.part0016.xtc 0/sys0_prod.part0016.xtc            
            mv 10/sys10_prod.part0016.log 0/sys0_prod.part0016.log
            mv 10/sys10_prod.part0016.edr 0/sys0_prod.part0016.edr
        fi
        
        mv 10/sys10_prod.cpt 0/sys0_prod.cpt
        mv 10/sys10_prod_prev.cpt 0/sys0_prod_prev.cpt
        cd ../
    done    
}

move_64_directories