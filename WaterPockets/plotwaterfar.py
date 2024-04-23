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
# from scipy import stats
# acetylationlist = [0,5,10,15,18]
# # frame = dict()
# # water = dict()
# frame = []
# water = []
# fig, ax = plt.subplots(len(acetylationlist), 1, figsize=(8,6), sharex=True)
# plt.subplots_adjust(left=0.12, right=0.84, wspace=0, hspace=0)
# for i, degree in enumerate(acetylationlist):
#     datawater = np.load("npz_cont/contacts_waterfar_%dpercent.npz" % (degree))
#     frame.append(datawater['arr_0'])
#     water.append(datawater['arr_1'])
#     ax.plot(frame[i], water[i],label=f'{degree}% Acetylation')
# ax.set_ylabel('Avg. number of contacts')
# ax.legend()
# ax.set_xlabel("Frame count")
# plt.savefig("results/water_pockets_percent_combined.png", dpi=300)
    # print(frame, water)
# print(frame, water)
# exit()
# for i in range(len(water)):
#     plt.plot(frame[i],water[i])
#     # plt.legend(loc='best')
#     # plt.ylabel('Avg. number of contacts')
#     # plt.xlabel("Frame count")
#     # plt.tight_layout()
# plt.savefig("results/water_pockets_percent.png", dpi=300)
acetylationlist = [0, 5, 10, 15, 18]
frame = []
water = []

fig, ax = plt.subplots(figsize=(8, 6))
plt.subplots_adjust(left=0.12, right=0.84, wspace=0, hspace=0)

for i, degree in enumerate(acetylationlist):
    datawater = np.load("npz_cont/contacts_waterfar%dpercent.npz" % (degree))
    frame_data=(datawater['arr_0'])
    water_data=(datawater['arr_1'])
    frame.append(frame_data[100:])
    water.append(water_data[100:])
    ax.plot(frame[i], water[i], label=f'{degree}% Acetylation')

ax.set_xlabel("Time (ns)")
# ax.set_title("Impact of acetylation on Water pockets")
ax.set_ylabel('Number of water molecules in water pockets')
ax.set_xticks([0,200,400,600,800,1000])
#plt.title('Water Pockets (5 Å away from cell wall components)')
ax.legend()
plt.savefig("results/water_pockets_5ang.png", dpi=100)
plt.savefig("results/water_pockets_5ang.pdf")


