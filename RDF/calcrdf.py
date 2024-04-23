from vmd import Molecule, atomsel
from vmd import vmdnumpy as vnp
import numpy as np
import glob
from vmd import evaltcl
from scipy import stats
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
evaltcl("package require pbctools")
def loadtraj(directory):
    mol = Molecule.Molecule()
    mol.load("%s/acetylatedsystem.psf" % directory)
    dcdlist = sorted(glob.glob("%s/run*dcd" % directory))
    # Now, doing this for last 10 trajectory files - since system has equilibrated
    for dcd in dcdlist[-10:]: 
    	mol.load(dcd, step=10)
    return mol
def unwrap():
	evaltcl("fastpbc unwrap [atomselect top \"all\"]")
def getrdf(mol, seltext):
	# print('Hey')
	Rref = vnp.timestep(int(mol))
	# print('hey')
	evaltcl("set celsel [atomselect %d \"%s\"]" % (int(mol), seltext))
	# evaltcl("set wsel [atomselect %d \"water\"]" % (int(mol)))
	evaltcl("set isel [atomselect %d \"resname SOD\"]" % (int(mol)))
	rdfoutput = evaltcl("measure gofr $celsel $wsel rmax 10 usepbc 1 first 0 last %d" % (mol.numFrames()-1))
	# print(rdfoutput)
	outputlist = rdfoutput[1:-1].split("} {")
	# print(outputlist)
	# exit()
	numpyoutputlists = []
	for output in outputlist:
		numpyoutputlists.append(np.array(output.split(), dtype=float))
	output = np.vstack(numpyoutputlists[:-1]).T
	return output
	
acetylationlist = ['5percent','10percent','15percent','18percent']
for degree in acetylationlist:
	mol = loadtraj("../%s" %(degree))
	unwrap()
	# for name, seltext in zip(['lignin','hemicellulose','cellulose','lignin-acetyl','hemicellulose-acetyl','lignin-hydroxyl','hemicellulose-hydroxyl', 'carbonyl-oxygen','bridge-oxygen','methyl-carbon','carbonyl-carbon'], ['segname \\\"XX.*\\\"','segname \\\"XY.*\\\"','segname \\\"CEL.*\\\"and (x*x + y*y) > 18*18','segname \\\"XX.*\\\" and withinbonds 1 of name \\\".AC.*\\\" \\\".AY.*\\\"','segname \\\"XY.*\\\" and withinbonds 1 of name \\\".AC.*\\\" \\\".AY.*\\\"','segname \\\"XX.*\\\" and withinbonds 1 of name \\\"HO.*\\\"', 'segname \\\"XY.*\\\" and withinbonds 1 of name \\\"HO.*\\\"','name \\\"OAC.*\\\"','name \\\"O.*\\\" and (withinbonds 1 of name \\\".AC.*\\\" \\\".AY.*\\\")','name \\\"CAY.*\\\"', 'name \\\"CAC.*\\\"']):
	for name, seltext in zip(['carbonyl-oxygen','bridge-oxygen','methyl-carbon','carbonyl-carbon'], ['name \\\"OAC.*\\\"','type OSL and (withinbonds 1 of name \\\".AC.*\\\" \\\".AY.*\\\")','name \\\"CAY.*\\\"', 'name \\\"CAC.*\\\"']):	
		result = getrdf(mol, seltext)
		np.save("npz_rdf/%s-cationrdf-%s.npy"% (name, degree), result)
	mol.delete()
exit()







