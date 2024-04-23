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
def contactsfromballtree(balltree1, balltree2, pbcvec):
	query = balltree1.query_ball_tree(balltree2, 5)
	contact = 0
	for i, resultlist in enumerate(query):
		for j in resultlist:
			vec = shortestvector(balltree1.data[i]-balltree2.data[j], pbcvec)
			distance = np.linalg.norm(vec)
			contact += 1 / (1+ np.exp(5*(distance - 3)))
	return contact
def calccontacts(mol):
	cellsel = atomsel("noh and segname \"CEL.*\"")
	ligsel = atomsel('noh and segname \"XX.*\"')
	hemisel = atomsel('noh and segname \"XY.*\"')
	ionpsel = atomsel ('resname SOD')
	ligacetsel = atomsel ('noh and segname \"XX.*\" and (withinbonds 1 of name \".AC.*\" \".AY.*\")')
	hemiacetsel = atomsel ('noh and segname \"XY.*\" and (withinbonds 1 of name \".AC.*\" \".AY.*\")')
	lighydroxysel = atomsel('noh and segname \"XX.*\" and (withinbonds 1 of name \"HO.*\")')
	hemihydroxysel =atomsel('noh and segname \"XY.*\" and (withinbonds 1 of name \"HO.*\")')
	oacetylsel =atomsel('noh and segname \"XY.*\" \"XX.*\" and (name \"OAC.*\" \"OAY.*\")')
	cacetylsel =atomsel('noh and segname \"XY.*\" \"XX.*\" and (name \"CAC.*\")')
	cacetylmsel =atomsel('noh and segname \"XY.*\" \"XX.*\" and (name \"CAY.*\")')
	ohacetylsel = atomsel('noh and name \"O.\" and (withinbonds 1 of name \".AC.*\" \".AY.*\") ')
	ohydroxysel = atomsel('noh and segname \"XY.*\" \"XX.*\" and (name \"O.*\" and withinbonds 1 of name \"HO.*\")')
	hacetylsel = atomsel('noh and segname \"XY.*\" \"XX.*\" and (withinbonds 1 of type \"HA.*\")')
	# watersel = atomsel('water')
	cellidx = np.array(cellsel.index).astype(int)
	ligidx = np.array(ligsel.index).astype(int)
	hemiidx = np.array(hemisel.index).astype(int)
	ionpidx = np.array(ionpsel.index).astype(int)
	aligidx = np.array(ligacetsel.index).astype(int)
	ahemidx = np.array(hemiacetsel.index).astype(int)
	ohligidx = np.array(lighydroxysel.index).astype(int)
	ohhemiidx = np.array(hemihydroxysel.index).astype(int)
	oacetylidx = np.array(oacetylsel.index).astype(int)
	cacetylidx = np.array(cacetylsel.index).astype(int)
	cacetylmidx = np.array(cacetylmsel.index).astype(int)
	hacetylidx= np.array(hacetylsel.index).astype(int)
	ohydroxyidx = np.array(ohydroxysel.index).astype(int)
	ohacetylidx = np.array(ohacetylsel.index).astype(int)
	# wateridx = np.array(watersel.index).astype(int)
	clcontacts = np.empty(mol.numFrames(), dtype=float)
	chcontacts = np.empty(mol.numFrames(), dtype=float)
	lhcontacts = np.empty(mol.numFrames(), dtype=float)
	ipccontacts = np.empty(mol.numFrames(), dtype=float)
	iphcontacts = np.empty(mol.numFrames(), dtype=float)
	iplcontacts = np.empty(mol.numFrames(), dtype=float)
	ipalcontacts = np.empty(mol.numFrames(), dtype=float)
	ipahcontacts = np.empty(mol.numFrames(), dtype=float)
	ipolcontacts = np.empty(mol.numFrames(), dtype=float)
	ipohcontacts = np.empty(mol.numFrames(), dtype=float)
	iacocontacts = np.empty(mol.numFrames(), dtype=float)
	iacccontacts = np.empty(mol.numFrames(), dtype=float)
	iaccmcontacts = np.empty(mol.numFrames(), dtype=float)
	iachcontacts = np.empty(mol.numFrames(), dtype=float)
	iohycontacts = np.empty(mol.numFrames(), dtype=float)
	iohaccontacts = np.empty(mol.numFrames(),dtype=float)
	# wacocontacts = np.empty(mol.numFrames(), dtype=float)
	# wohycontacts = np.empty(mol.numFrames(), dtype=float)
	# wwcontacts = np.empty(mol.numFrames(),dtype=float)
	x = 0.5 * np.arange(mol.numFrames())
	for f in range(mol.numFrames()):
		pbc = molecule.get_periodic(molid=int(mol), frame=f)
		pbcvec = np.array([pbc['a'], pbc['b'], pbc['c']])
		#threshold = np.finfo(float32).eps
		threshold = 1e-3
		pbcvec = np.array([pbc['a']+threshold, pbc['b']+threshold, pbc['c']+threshold])
		R = vnp.timestep(int(mol), f)
		R += 0.5 * pbcvec
		R[R<0] = 0
		Rcell = R[cellidx]
		Rhemi = R[hemiidx]
		Rlig = R[ligidx]
		Rion = R[ionpidx]
		Ralig = R[aligidx]
		Rahem = R[ahemidx]
		Rolig = R[ohligidx]
		Rohem = R[ohhemiidx]
		Roace = R[oacetylidx]
		Rhace = R[hacetylidx]
		Rcace = R[cacetylidx]
		Rcacem = R[cacetylmidx]
		Rohyd = R[ohydroxyidx]
		Rohace = R[ohacetylidx]
		# Rwater = R[wateridx]
		celltree = KDTree(Rcell, boxsize=pbcvec)
		hemitree = KDTree(Rhemi, boxsize=pbcvec)
		ligtree = KDTree(Rlig, boxsize=pbcvec)
		ionptree = KDTree(Rion, boxsize=pbcvec)
		aligtree = KDTree(Ralig, boxsize=pbcvec)
		ahemtree = KDTree(Rahem, boxsize=pbcvec)
		oligtree = KDTree(Rolig, boxsize=pbcvec)
		ohemtree = KDTree(Rohem, boxsize=pbcvec)
		oacetyltree = KDTree(Roace, boxsize=pbcvec)
		cacetyltree = KDTree(Rcace, boxsize=pbcvec)
		cacetylmtree = KDTree(Rcacem, boxsize = pbcvec)
		hacetyltree = KDTree(Rhace, boxsize=pbcvec)
		ohydroxytree = KDTree(Rohyd, boxsize=pbcvec)
		ohacetyltree = KDTree(Rohace, boxsize=pbcvec)
		# watertree = KDTree(Rwater, boxsize=pbcvec) 
		clcontacts[f] = contactsfromballtree(celltree, ligtree, pbcvec)
		chcontacts[f] = contactsfromballtree(celltree, hemitree, pbcvec)
		lhcontacts[f] = contactsfromballtree(ligtree, hemitree, pbcvec)
		ipccontacts[f] = contactsfromballtree(ionptree, celltree, pbcvec)
		iphcontacts[f] = contactsfromballtree(ionptree, hemitree, pbcvec)
		iplcontacts[f] = contactsfromballtree(ionptree, ligtree, pbcvec)
		ipalcontacts[f] = contactsfromballtree(ionptree, aligtree, pbcvec)
		ipahcontacts[f] = contactsfromballtree(ionptree, ahemtree, pbcvec)
		ipolcontacts[f] = contactsfromballtree(ionptree, oligtree, pbcvec)
		ipohcontacts[f] = contactsfromballtree(ionptree, ohemtree, pbcvec)
		iacocontacts[f] = contactsfromballtree(ionptree, oacetyltree, pbcvec)
		iacccontacts[f] = contactsfromballtree(ionptree, cacetyltree, pbcvec)
		iaccmcontacts[f] = contactsfromballtree(ionptree, cacetylmtree, pbcvec)
		iachcontacts[f] = contactsfromballtree(ionptree, hacetyltree, pbcvec)
		iohycontacts[f] = contactsfromballtree(ionptree, ohydroxytree, pbcvec)
		iohaccontacts[f] = contactsfromballtree(ionptree, ohacetyltree, pbcvec)
		# wacocontacts[f] = contactsfromballtree(watertree, oacetyltree, pbcvec)
		# wohycontacts[f] = contactsfromballtree(watertree, ohydroxytree, pbcvec)
		# wwcontacts[f] = contactsfromballtree(watertree, watertree, pbcvec)
		# print(f, mol.numFrames())
	return (x, clcontacts, chcontacts, lhcontacts, ipccontacts, iphcontacts, iplcontacts, ipalcontacts, ipahcontacts, ipolcontacts, ipohcontacts, iacocontacts, iacccontacts, iaccmcontacts,iachcontacts, iohycontacts, iohaccontacts)
	#return (x, ipccontacts, iphcontacts, iplcontacts)
acetylationlist = ['0percent','5percent', '10percent', '15percent', '18percent']
for degree in acetylationlist:
	#print("../"+str(degree))
	mol = loadtraj("../"+str(degree))
	wrap(mol)
	x, clcontacts, chcontacts, lhcontacts, ipccontacts, iphcontacts, iplcontacts, ipalcontacts, ipahcontacts, ipolcontacts, ipohcontacts, iacocontacts, iacccontacts, iaccmcontacts, iachcontacts, iohycontacts, iohaccontacts = calccontacts(mol)
	#x, ipccontacts, iphcontacts, iplcontacts = calccontacts(mol)
	np.savez("npz_cont/contacts_3ang_cel_lig_%s.npz" % (degree), x, clcontacts)
	np.savez("npz_cont/contacts_3ang_cel_hem_%s.npz" % (degree), x, chcontacts)
	np.savez("npz_cont/contacts_3ang_lig_hem_%s.npz" % (degree), x, lhcontacts)
	np.savez("npz_cont/contacts_3ang_ion_cel_%s.npz" % (degree), x, ipccontacts)
	np.savez("npz_cont/contacts_3ang_ion_hem_%s.npz" % (degree), x, iphcontacts)
	np.savez("npz_cont/contacts_3ang_ion_lig_%s.npz" % (degree), x, iplcontacts)
	np.savez("npz_cont/contacts_3ang_ion_alig_%s.npz" % (degree), x, ipalcontacts)
	np.savez("npz_cont/contacts_3ang_ion_ahem_%s.npz" % (degree), x, ipahcontacts)
	np.savez("npz_cont/contacts_3ang_ion_olig_%s.npz" % (degree), x, ipolcontacts)
	np.savez("npz_cont/contacts_3ang_ion_ohem_%s.npz" % (degree), x, ipohcontacts)
	np.savez("npz_cont/contacts_3ang_ion_oacetyl_%s.npz" % (degree), x, iacocontacts)
	np.savez("npz_cont/contacts_3ang_ion_cacetyl_%s.npz"%(degree),x,iacccontacts)
	np.savez("npz_cont/contacts_3ang_ion_cacetylm_%s.npz"%(degree), x, iaccmcontacts)
	np.savez("npz_cont/contacts_3ang_ion_hacetyl_%s.npz" % (degree), x, iachcontacts)
	np.savez("npz_cont/contacts_3ang_ion_ohydroxy_%s.npz" % (degree), x, iohycontacts)
	np.savez("npz_cont/contacts_3ang_ion_ohacetyl_%s.npz" % (degree), x, iohaccontacts)
	# np.savez("npz_cont/contacts_3ang_water_oacetyl_%s.npz" % (degree), x, wacocontacts)
	# np.savez("npz_cont/contacts_3ang_water_ohydroxy_%s.npz" % (degree), x, wohycontacts)
	# np.savez("npz_cont/contacts_3ang_water_water_%s.npz" % (degree), x, wwcontacts)
	mol.delete()

exit()