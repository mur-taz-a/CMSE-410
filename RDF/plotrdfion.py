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
    'font.size' : 16,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    'legend.fontsize': 12,
    'legend.handlelength': 2})
import matplotlib.pyplot as plt
import numpy as np

molecules = ['ion']
percentlist = ['0percent','5percent', '10percent', '15percent', '18percent']
fig, ax = plt.subplots(1, len(molecules), figsize=(10,6), sharex=True, sharey=True)
plt.subplots_adjust(left=0.12, right=0.88, wspace=0, hspace=0)
for i, name in enumerate(molecules):
	for j, percent in enumerate(percentlist):
		print("npz_rdf/%s-waterRDF-%s.npy"% (name, percent))
		data = np.load("npz_rdf/%s-waterRDF-%s.npy"% (name, percent))
		binvol = (4.0*np.pi*(data[1:,0]**3 - data[:-1,0]**3))
		binvol = np.hstack([[0],binvol])
		ax.plot(data[:,0], data[:,3]/binvol, color="C%s" % j, label="%s" % percent)
		ax.set_xlabel("r (\AA)")
		ax.set_title("%s" % name.capitalize(), fontsize = 15)
ax.set_xlim([0, 30])
ax.legend(ncol=2)
fig.text(0.01, 0.35, "g(r) (counts/\AA$^3$)", rotation='vertical')
fig.savefig("results/rdfion.png", dpi=300)
exit()