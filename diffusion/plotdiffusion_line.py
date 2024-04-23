#/usr/bin/env python
import matplotlib
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}",
r"\usepackage{fixltx2e}",
	   r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
	   r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
	   r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
	   r'\sansmath' ,              # <- tricky! -- gotta actually tell tex to use!
	   r'\usepackage{textgreek}'
]

matplotlib.rcParams.update({
	"text.usetex": True,
	'font.size' : 18,
	"font.family": "sans-serif",
	"font.sans-serif": ["Helvetica"]})
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.stats import linregress
ticklist = ["o", "v", "<", "s", "p", "P", "*", "h"]
# def getname(i, j):
# 	if j == 0:
# 		if i == 0:
# 			return "K$^+$"
# 		else:
# 			return "Cu$^{2+}$"
# 	elif j == 1:
# 		return "Cl$^-$"
# 	else:
# 		print("This should never happen")
# 		exit()
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0, x>=x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
def plotalllog():
	tickcounter = 0
	fig, ax = plt.subplots(1,1)
	ionlist = ['sodium', 'fe3p']
	acetylationlist = [0,5,10,15,18]
	for i, name in enumerate(['water', 'lignin', 'hemicellulose', 'cellulose', 'ion']):
		plotarray = []
		plotpercents = []
		ploterror = []
		plotarray2 = []
		for degree in acetylationlist:
			data = np.load("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
			plotarray.append(data['altdiffusion'][0]*1e-7)
			ploterror.append(data['altdiffusion'][1]*1e-7)
			if name == 'water':
				print("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
				data = np.load("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
				x = data['x']
				msds = data['msds']
				#print(name, percent, len(x))
				result = stats.linregress(x[10:], msds[10:])
				diffusioncoeff = result.slope/6 * (1e-7)
				plotarray2.append(diffusioncoeff)
				#plotpercents.append(percent)
			else:
				print("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
				data = np.load("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
				msds = data['msds2']
				deltatime = 0.5 *100
				diffusioncoeff = np.mean(msds[100:])/(deltatime*6) * (1e-7)
				plotarray2.append(diffusioncoeff)
			plotpercents.append(degree)
		if name !='ion':
			print(plotpercents, plotarray2, ploterror)
			ax.errorbar(plotpercents, plotarray2, yerr=ploterror, label=name.capitalize(), marker=ticklist[tickcounter], color="C%d" % i)
			#ax.plot(plotpercents, plotarray2, color="C%d" % i)
		else:
			ax.errorbar(plotpercents, plotarray2, yerr=ploterror, label="Na$^{+}$", marker=ticklist[tickcounter], color="C%d" % tickcounter)
			#ax.plot(plotpercents, plotarray2, color="C%d" % i)
		tickcounter += 1
	# for i, salt in enumerate(['KCl', 'CuCl2']):
	# 	for j, name in enumerate(['ion+', 'ion-']):
	# 		if i == 1 and j == 1:
	# 			continue
	# 		plotarray = []
	# 		plotpercents = []
	# 		ploterror = []
	# 		for percent in percentlist[1:]:	
	# 			print("msds_%s_%dpercent.npz" % (name, percent))
	# 			data = np.load("../../Simulations_JZ/Analysis/msds/%s_%s_%dpercent.npz" % (salt, name, percent))
	# 			plotarray.append(data['altdiffusion'][0]*1e-7)
	# 			ploterror.append(data['altdiffusion'][1]*1e-7)
	# 			plotpercents.append(percent)
	# 		ax.errorbar(plotpercents, plotarray, yerr=ploterror, label=getname(i,j), marker=ticklist[tickcounter], color="C%d" % tickcounter)
	# 		tickcounter += 1
	ax.set_xticks(acetylationlist)
	ax.tick_params(axis='x', which='major', labelsize=15)
	# labels = [0,5,10,15,18]
	# ax.set_xticklabels(labels)
	ax.set_xlabel("\% Weight Gain")
	ax.set_ylabel("D (cm$^2$/s)")
	ax.set_yscale("log")
	plt.legend(loc='best', fontsize='xx-small', ncol=2)
	fig.tight_layout()
	fig.savefig("results/logdiffusion.png",dpi=300)
	fig.savefig("results/logdiffusion.pdf")
	# fig.savefig("results/diffusion.png", dpi=300)


# def plotcomponents(enumeratelist, filename):
# 	fig, ax = plt.subplots(1,1)
# 	acetylationlist = [0,5,10,15,18]
# 	for i, name in enumerate(enumeratelist):
# 		plotarray = []
# 		ploterror = []
# 		plotpercents = []
# 		for degree in acetylationlist:
# 			if name == 'water':	
# 				print("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
# 				data = np.load("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
# 				x = data['x']
# 				msds = data['msds']
# 				#print(name, percent, len(x))
# 				result = stats.linregress(x[100:], msds[100:])
# 				diffusioncoeff = result.slope/6 * (10**-7)
# 				deltatime = 0.5 *100
# 				#This print line compares all the different ways of estimating diffusion.
# 				print(diffusioncoeff, data['altdiffusion'][0]*1e-7, np.mean(data['msds'][100:])/(deltatime*6) * (10**-7))
# 				plotarray.append(data['altdiffusion'][0]*1e-7)
# 				ploterror.append(data['altdiffusion'][1]*1e-7)
# 				plotpercents.append(degree)
# 			else:
# 				print("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
# 				data = np.load("npz_diff/msds_NaCl_30percent_%s_%dpercent.npz" % (name, degree))
# 				print(data)
# 				x = data['x2']
# 				msds = data['msds2']
# 				#print(name, percent, len(x))
# 				# result = stats.linregress(x[100:], msds[100:])
# 				deltatime = 0.5 *100
# 				diffusioncoeff = np.mean(msds[100:])/(deltatime*6) * (10**-7)
# 				#This print line compares all the different ways of estimating diffusion.
# 				print(diffusioncoeff, data['altdiffusion'][0]*1e-7, np.mean(data['msds2'][100:])/(deltatime*6) * (10**-7))
# 				plotarray.append(data['altdiffusion'][0]*1e-7)
# 				ploterror.append(data['altdiffusion'][1]*1e-7)
# 				plotpercents.append(degree)
# 			# else:
# 			# 	print("msds_%s_%dpercent.npz" % (name, percent))
# 			# 	data = np.load("msds_%s_%dpercent.npz" % (name, percent))
# 			# 	msds = data['msds2']
# 			# 	#print(name, percent, len(x))
# 			# 	#result = stats.linregress(x[1:], msds[1:])
# 			# 	#8 loaded frames equivalent to 4 ns, 2 frames per ns, each frame is 500ps.
# 			# 	#calcdiffusion is using an offset of 100.
# 			# 	deltatime = 0.5 *100
# 			# 	diffusioncoeff = np.mean(msds2[100:])/(deltatime*6) * (10**-7)
# 			# 	plotarray.append(diffusioncoeff)
# 			# 	plotpercents.append(percent)
# 		print(plotarray, len(plotarray), len(acetylationlist))
# 		if name !='ion':
# 			#ax.plot(plotpercents, plotarray, label=name.capitalize(), marker='.')
# 			ax.errorbar(plotpercents, plotarray, yerr=ploterror, fmt='o', label=name.capitalize(), marker='.', zorder=2)
# 		else:
# 			#ax.plot(plotpercents, plotarray, label="Na$^{+}$", marker='.')
# 			ax.errorbar(plotpercents, plotarray, yerr=ploterror, fmt='o', label="Na$^{+}$", marker='.', zorder=2)
# 		if name != 'water':
# 			m, b, r, p, stderr = linregress(plotpercents[:2], plotarray[:2])
# 		else:
# 			# for water since there is no data for 0 weight%
# 			m, b, r, p, stderr = linregress(plotpercents[1:3], plotarray[1:3])
# 		#print("First part", i, j, r*r)
# 		ax.plot([plotpercents[0],plotpercents[3]], m*np.array([plotpercents[0],plotpercents[3]]) + b, color='C%d' %i, linestyle="--",linewidth=1, zorder = 1)
# 		if i == 0:
# 			if name =='water':
# 				ax.text(0.2, 0.15, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 			else:
# 				ax.text(0.2, 0.2, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 		elif i == 2:
# 			ax.text(0.2, 0.25, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 		else:
# 			ax.text(0.2, 0.3, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 		m, b, r, p, stderr = linregress(plotpercents[3:], plotarray[3:])
# 		ax.plot([plotpercents[3],plotpercents[-1]], m*np.array([plotpercents[3],plotpercents[-1]]) + b, color='C%d' %i, linestyle="--",linewidth=1, zorder = 1)
# 		if i == 0:
# 			if name == 'water':
# 				ax.text(0.85, 0.8, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 			else:
# 				ax.text(0.85, 0.1, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 		elif i == 2:
# 			ax.text(0.85, 0.72, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 		else:
# 			ax.text(0.85, 0.48, "$r^2=%.2f$" % (r*r), transform=ax.transAxes, ha='right', fontsize=9, color='C%d' %i)
# 	# ax.set_xticks(acetylationlist)
# 	ax.tick_params(axis='x', which='major', labelsize=15)
# 	# labels = [0,5,10,15,18]
# 	# ax.set_xticklabels(labels)
# 	ax.set_xlabel("Water Weight (\%)")
# 	ax.set_ylabel("Diffusion Constant, D (cm$^2$/s)")
# 	#ax.set_yscale("log")
# 	plt.legend(loc='best', fontsize='xx-small')
# 	fig.tight_layout()
# 	#fig.savefig("DiffusionIonNaCl.png", dpi=300)
# 	fig.savefig(filename, dpi=300)
plotalllog()
# plotcomponents(['lignin','cellulose', 'hemicellulose', 'water', 'ion'], "Diffusioncomponents.png")
# plotcomponents(['water'], "DiffusionWater.png")
# plotcomponents(['ion'], "DiffusionIonNaCl.png")