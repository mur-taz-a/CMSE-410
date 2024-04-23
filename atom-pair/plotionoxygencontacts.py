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
    "lines.linewidth": 2})

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
acetylationlist = [0,5,10,15,18]
hydroxyllist = [4098,2905,1931,738,0]
acetyllist = [0,1193,2197,3360,4098]
# clmean = dict()
# chmean = dict()
# lhmean = dict()
# ipcmean = dict()
# iphmean = dict()
# iplmean = dict()
# ipalmean = dict()
# ipahmean = dict()
iaceomean = dict()
ihydomean = dict()
feacemean = dict()
fehydmean = dict()
# iacehmean = dict()
# ipolmean = dict()
# ipohmean = dict()

# clmean = list()
# chmean = list()
# lhmean = list()
# ipcmean = list()
# iphmean = list()
# iplmean = list()
# ipalmean = list()
# ipahmean = list()
iaceomean = list()
ihydomean = list()
feacemean = list()
fehydmean = list()
# iacehmean = list()
# ipolmean = list()
# ipohmean = list()
    
fig, ax = plt.subplots(len(acetylationlist), 1, figsize=(8,6), sharex=True)
plt.subplots_adjust(left=0.12, right=0.84, wspace=0, hspace=0)
for i, degree in enumerate(acetylationlist):
    #datacl = np.load("cel_lig_%dpercent.npz" % (percent))
    #datach = np.load("cel_hem_%dpercent.npz" % (percent))
    #datalh = np.load("lig_hem_%dpercent.npz" % (percent))
    # dataipc = np.load("npz_cont/contacts_ion_cel_%dpercent.npz" % (degree))
    # dataiph = np.load("npz_cont/contacts_ion_hem_%dpercent.npz" % (degree))
    # dataipl = np.load("npz_cont/contacts_ion_lig_%dpercent.npz" % (degree))
    # dataipal = np.load("npz_cont/contacts_3ang_ion_alig_%dpercent.npz" % (degree))
    # dataipah = np.load("npz_cont/contacts_3ang_ion_ahem_%dpercent.npz" % (degree))
    dataiaco = np.load("npz_cont/contacts_3ang_ion_oacetyl_%dpercent.npz" % (degree))
    dataihyo = np.load("npz_cont/contacts_3ang_ion_ohydroxy_%dpercent.npz" % (degree))
    datafeaco = np.load("../fe3p/simulation/analysis/calc/contacts_3ang_ion_oacetyl_%d.npz" % (degree))
    datafehyo = np.load("../fe3p/simulation/analysis/calc/contacts_3ang_ion_ohydroxy_%d.npz" % (degree))
    # dataipol = np.load("npz_cont/contacts_ion_olig_%dpercent.npz" % (degree))
    # dataipoh = np.load("npz_cont/contacts_ion_ohem_%dpercent.npz" % (degree))
    #x = datacl['arr_0']
    #clcontacts = datacl['arr_1']
    #chcontacts = datach['arr_1']
    #lhcontacts = datalh['arr_1']
    # x = dataipc['arr_0']
    # ipccontacts = dataipc['arr_1']
    # iphcontacts = dataiph['arr_1']
    # iplcontacts = dataipl['arr_1']
    # x = dataipal['arr_0']
    # ipalcontacts = dataipal['arr_1']
    # ipahcontacts = dataipah['arr_1']
    # x = dataipol['arr_0']
    # ipolcontacts = dataipol['arr_1']
    # ipohcontacts = dataipoh['arr_1']
    iaceocontacts = dataiaco['arr_1']
    ihydocontacts = dataihyo['arr_1']
    feaceocontacts = datafeaco['arr_1']
    fehydocontacts = datafehyo['arr_1']
    #clmean.append(np.mean(clcontacts[200:]))
    #chmean.append(np.mean(chcontacts[200:]))
    #lhmean.append(np.mean(lhcontacts[200:]))
    # ipcmean.append(np.mean(ipccontacts[200:]))
    # iphmean.append(np.mean(iphcontacts[200:]))
    # iplmean.append(np.mean(iplcontacts[200:]))
    # ipalmean.append(np.mean(ipalcontacts[200:]))
    # ipahmean.append(np.mean(ipahcontacts[200:]))
    iaceomean.append(np.mean(iaceocontacts[200:]))
    ihydomean.append(np.mean(ihydocontacts[200:]))
    feacemean.append(np.mean(feaceocontacts[200:]))
    fehydmean.append(np.mean(fehydocontacts[200:]))
    # iacehmean.append(np.mean(iacehcontacts[200:]))
    # ipolmean.append(np.mean(ipolcontacts[200:]))
    # ipohmean.append(np.mean(ipohcontacts[200:]))
    # ax[i].hist(ipccontacts[200:], bins=50, color="blue", label='Ion-Cellulose')
    # ax[i].hist(iphcontacts[200:], bins=50, color="green", label='Ion-Hemicellulose')
    # ax[i].hist(iplcontacts[200:], bins=50, color="orange", label='Ion-Lignin')
    # ax[i].set_xlabel("Number of contacts")
# ax[i].legend(loc='best', fontsize='xx-small')
# fig.text(0.03, 0.4, "Bin counts", rotation='vertical')
# fig.text(0.85, 0.80, "3 Water Wt. (\%)" ,fontsize=10)
# fig.text(0.85, 0.70, "5 Water Wt. (\%)" ,fontsize=10)
# fig.text(0.85, 0.60, "10 Water Wt. (\%)" ,fontsize=10)
# fig.text(0.85, 0.50, "15 Water Wt. (\%)" ,fontsize=10)
# fig.text(0.85, 0.35, "20 Water Wt. (\%)" ,fontsize=10)
# fig.text(0.85, 0.25, "25 Water Wt. (\%)" ,fontsize=10)
# fig.text(0.85, 0.15, "30 Water Wt. (\%)" ,fontsize=10)
# fig.savefig("contactsalltracenew_{}.png".format('NaCl'), dpi=300)

fig, ax = plt.subplots(1, 1)
# plt.plot(acetylationlist,ipcmean, '-o', label='Ion-Cellulose')
# plt.plot(acetylationlist,iphmean, '-o' ,label='Ion-Hemicellulose')
# plt.plot(acetylationlist,iplmean, '-o', label='Ion-Lignin')
# plt.plot(acetylationlist, ipalmean,'-o', label='Ion-AcetylLignin',color = 'red')
# plt.plot(acetylationlist, ipahmean,'-v', label='Ion-AcetylHemi',color ='red')
plt.bar(np.array(acetylationlist)-0.5, iaceomean, label='Na$^{+}$ - Acetyl',color='dodgerblue')
plt.bar(np.array(acetylationlist)-0.5, ihydomean, bottom = iaceomean, label = 'Na$^{+}$ - Hydroxyl', color='orangered')
plt.bar(np.array(acetylationlist)+0.5, feacemean, label='Fe$^{3+}$ - Acetyl', color='darkblue')
plt.bar(np.array(acetylationlist)+0.5, fehydmean, bottom=feacemean, label='Fe$^{3+}$ - Hydroxyl', color='darkred')
# plt.plot(acetylationlist, iacehmean, '-2', label='IOn-Hacetyl', color='purple')
# plt.plot(acetylationlist, ipolmean, '-o',label='Ion-hydroxylLignin', color ='green')
# plt.plot(acetylationlist, ipohmean,'-v', label='Ion-hydroxylHemi',  color = 'green')
# labels = [0,5,10,15,18]
ax.set_xticks(acetylationlist)
plt.xlabel('\% Weight Gain ')
plt.ylabel('Avg. number of contacts')
# plt.title('Oxygen Ion contacts in different groups')
plt.legend(loc='best')
fig.tight_layout()
plt.savefig("results/oxygen_ion.png",dpi=300)
plt.savefig("results/oxygen_ion.pdf")

# iaceprop= np.array(iaceomean)/np.array(acetyllist)
# ihydprop=np.array(ihydomean)/np.array(hydroxyllist)
# print(iaceprop, ihydprop)
# fig,ax = plt.subplots(1,1)
# plt.plot(acetylationlist,iaceprop,'-o', label='Acetyl',color='brown')
# plt.plot(acetylationlist,ihydprop,'-x', label='Hydroxyl',color='purple')
# plt.plot(acetylationlist,np.add(iaceomean,ihydomean)/4098,'-1',label='Cumulative',color='grey')
# ax.set_xticks(acetylationlist)
# plt.xlabel('Degree of acetylation (\%)')
# plt.ylabel('\%\ of groups interacting')
# plt.title('Oxygen Ion contacts in different groups')
# plt.legend(loc='best')
# fig.tight_layout()
# plt.savefig("results/oxygen_ion_prop.png",dpi=300)
# plt.savefig("results/oxygen_ion_prop.pdf")