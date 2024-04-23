
from vmd import Molecule, atomsel
from vmd import vmdnumpy as vnp
import numpy as np
import glob
from vmd import evaltcl
from scipy import stats
evaltcl("package require pbctools")
def loadtraj(directory):
	mol = Molecule.Molecule()
	mol.load("%s/acetylatedsystem.psf" % directory)
	dcdlist = sorted(glob.glob("%s/run*dcd" % directory))
	#dcdlist = [dcdlist[0]] # this is just for testing, change it later to include all dcds.
	for dcd in dcdlist:
		print(dcd)
		# exit()
		mol.load(dcd, step=100)			# x value weight depends on stride length
	return mol
def unwrap(mol):
	evaltcl("fastpbc unwrap [atomselect top \"all\"]")
	ref = atomsel("segname \"CEL.*\"", frame=0)
	asel = atomsel("all")
	sel = atomsel("segname \"CEL.*\"")
	for f in range(mol.numFrames()):
		asel.frame=f
		sel.frame=f
		asel.move(sel.fit(ref))
def getwaterdiffusion(mol, seltext):
	Rref = vnp.timestep(int(mol), 0)
	idxs = np.nonzero(vnp.atomselect("noh and " + seltext,int(mol), 0))[0]
	msds = np.empty(mol.numFrames(), dtype=float)
	for f in range(mol.numFrames()):
		R = vnp.timestep(int(mol), f)
		displacement = np.sum(np.square(R[idxs] - Rref[idxs]), axis=1)
		msds[f] = np.mean(displacement)
	#8 frames equivalent to 4 ns, 2 frames per ns, each frame is 500ps.
	x = 0.5 * np.arange(mol.numFrames())
	res = stats.linregress(x, msds)
	#print(res.intercept, res.slope)
	#print(1e-7 * res.slope/6) #Units would be in A^2/ns... cm^2/s needs a multiplication by 10^9/10^16=1e-7
	return (x, msds, res)
def getdiffusion(mol, seltext):
	offset=100
	idxs = np.nonzero(vnp.atomselect("noh and " + seltext,int(mol), 0))[0]
	msds = np.empty(mol.numFrames(), dtype=float)
	for f in range(offset,mol.numFrames()):
		Rref = vnp.timestep(int(mol), f-offset)
		R = vnp.timestep(int(mol), f)
		displacement = np.sum(np.square(R[idxs] - Rref[idxs]), axis=1)
		msds[f] = np.mean(displacement)
	#8 frames equivalent to 4 ns, 2 frames per ns, each frame is 500ps.
	x = 0.5 * np.arange(mol.numFrames())
	res = stats.linregress(x, msds)
	#print(res.intercept, res.slope)
	#print(1e-7 * res.slope/6) #Units would be in A^2/ns... cm^2/s needs a multiplication by 10^9/10^16=1e-7
	return (x, msds, res)
def altgetdiffusion(mol, seltext):
	offset = 100
	idxs = np.nonzero(vnp.atomselect("noh and " + seltext,int(mol), 0))[0]
	meanslopes = []
	for f in range(offset,mol.numFrames(),offset):
		Rref = vnp.timestep(int(mol), f-offset)
		displacements = np.empty(offset, dtype=float)
		for i in range(offset):
			R = vnp.timestep(int(mol), f-offset+i)
			displacements[i] = np.mean(np.sum(np.square(R[idxs] - Rref[idxs]), axis=1))
		x = 0.5 * np.arange(offset)
		res = stats.linregress(x[int(offset*0.2):], displacements[int(offset*0.2):])
		meanslopes.append(res.slope/6)
	array = np.array(meanslopes)
	return (np.mean(array), np.std(array) / np.sqrt(len(array)))
acetylationlist = ['0percent','5percent', '10percent', '15percent', '18percent']
for degree in acetylationlist:
	mol = loadtraj("../%s" %(degree))
	unwrap(mol)
	for name, seltext in zip(['water', 'cellulose', 'lignin', 'hemicellulose', 'ion'],['water', 'segname \"CEL.*\"', 'segname \"XX.*\"', 'segname \"XY.*\"', 'resname SOD',]): 
		result2 = getdiffusion(mol, seltext)
		altdiffusion = altgetdiffusion(mol, seltext)
		result = getwaterdiffusion(mol, seltext) 
		if name == "water":
			np.savez("npz_diff/msds_NaCl_30percent_%s_%s.npz" % (name, degree), x=result[0], msds=result[1], altdiffusion=altdiffusion)
		else:
			np.savez("npz_diff/msds_NaCl_30percent_%s_%s.npz" % (name, degree), x2=result2[0], msds2=result2[1], altdiffusion=altdiffusion)
	mol.delete()
exit()
