# Usage: vmd -dispdev text -e applypatchpsfgen.tcl
package require psfgen
resetpsf
# Load molecule
mol load psf commons/system.psf pdb commons/system.pdb
# # Load topology files
topology top_all36_prot.rtf
topology top_all36_lipid.rtf
topology top_all36_carb.rtf
topology top_lignin.top
topology acetylationpatches.top

readpsf commons/system.psf
coordpdb commons/system.pdb
# puts "I'm here after reading topologies.."
# exit
# Proc to apply the specific patch based on the match between the polymer typer and hydroxyl hydrogen atom name.
# Removed selidx as an argument from the proc as we can load the entire pdb, instead of writng a fragment sel of the pdb: tmpfrag.pdb
proc applypatch {polyclass atomname psegname presid} {
	if { [string equal $polyclass XY] } {
		set patchname "H$atomname"
		# set segmentname "$patchname"
		puts "PATCHNAME: $patchname"
		# puts "SEGNAME: $segmentname"
	} elseif { [string equal $polyclass XX] } {
		set patchname "L$atomname"
		# set segmentname "$patchname"
		puts "PATCHNAME: $patchname"
		# puts "SEGNAME: $segmentname"
	} else {
		puts "Note, we are only doing this for lignin and hemicellulose. Check!"
	}
	
	patch $patchname $psegname:$presid 
	
}
# Make individual lists of each lignin and hemicellulose hydroxyl hydrogen atoms with solvent exposed hydroxyl groups
set lighosel [atomselect top "segname \"XX.*\" and name \"HO.*\""]
set lighoidxlist [$lighosel get index]
llength $lighoidxlist
set hemhosel [atomselect top "segname \"XY.*\" and name \"HO.*\""]
set hemhoidxlist [$hemhosel get index]
llength $hemhoidxlist
# exit
# After making the list we get the atom name for the given polymer type, as this will be used to check and apply the appropriate patch from the topology file.
# First for lignin
foreach ligselidx $lighoidxlist {
	set tmpsel [atomselect top "index $ligselidx"]
	set hatomname [$tmpsel get name]
	set hsegname [$tmpsel get segname]
	set hresid [$tmpsel get resid]
	set polyclassi [string range $hsegname 0 1]
	puts "$polyclassi $hatomname $hsegname $hresid"
	applypatch $polyclassi $hatomname $hsegname $hresid  
	# writepsf test.psf
	# writepdb test.pdb
	# exit
}
# And then for hemicellulose, but really order doesn't matter.
foreach hemselidx $hemhoidxlist {
	set tmpsel [atomselect top "index $hemselidx"]
	set hatomname [$tmpsel get name]
	set hsegname [$tmpsel get segname]
	set hresid [$tmpsel get resid]
	set polyclassi [string range $hsegname 0 1]
	puts "$polyclassi $hatomname $hsegname $hresid"
	applypatch $polyclassi $hatomname $hsegname $hresid 
}
regenerate angles dihedrals
guesscoord
writepsf acetylatedsystemfull.psf
writepdb acetylatedsystemfull.pdb
exit
