from vmd import Molecule, atomsel, molecule
from vmd import vmdnumpy as vnp
import numpy as np
import glob
from vmd import evaltcl
from scipy import stats
from scipy.spatial import KDTree
def shortestvector(Ri,pbc):
	return Ri - (np.round(Ri / pbc) * pbc)
evaltcl("package require pbctools")
def loadtraj(directory):
	mol = Molecule.Molecule()
	mol.load("%s/acetylatedsystem.psf" % directory)
	dcdlist = sorted(glob.glob("%s/run*dcd" % directory))
	for dcd in dcdlist:
		mol.load(dcd, step=100)			# x value weight depends on stride length; 1000 is just for checking if code runs
	return mol
def wrap(mol):
	evaltcl("fastpbc wrap [atomselect top \"all\"] center [list 0 0 0]")

def waterfarawaycounter(mol):
	sel = atomsel("noh and water and not pbwithin 5 of (noh and not (water or ion))")
	outarray = np.empty(mol.numFrames(), dtype=int)
	x = (0.5 * np.arange(mol.numFrames()))
	for f in range(mol.numFrames()):
		sel.frame = f
		sel.update()
		# print(f, len(sel))
		outarray[f] = len(sel)
	return (x, outarray)
acetylationlist = ['0percent','5percent', '10percent', '15percent', '18percent']
for degree in acetylationlist:
	#print("../"+str(degree))
	mol = loadtraj("../"+str(degree))
	x, outdata = waterfarawaycounter(mol)
	np.savez("npz_cont/contacts_waterfar%s.npz" % (degree), x, outdata)
	mol.delete()

exit()

