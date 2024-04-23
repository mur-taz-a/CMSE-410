# Usage vmd -e loadtrajectory.tcl

# no acetylation - old model
mol new ../../CellWallHydration-1/Simulations/30/system.psf
set dcd_list [lsort -increasing [glob ../../CellWallHydration-1/Simulations/30/run*dcd]]
foreach dcd $dcd_list {
	mol addfile $dcd step 1000 waitfor all 
}
mol rename top noacetylation

# my acetylation set 
set acetylationindex [list 5 10 15 18]
foreach dirname $acetylationindex {
	mol new ${dirname}percent/acetylatedsystem.psf
	set dcd_list [lsort -increasing [glob ${dirname}percent/run*dcd]]
	foreach dcd $dcd_list {
		mol addfile $dcd step 1000 waitfor all 
	}
	mol rename top ${dirname}percent
}

