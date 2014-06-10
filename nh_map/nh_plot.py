from math import *
import os
import time
import datetime
import numpy
import sys

import numpy as np
import matplotlib as plt

from pylab import *
from lxml import html


fig = figure(figsize=(10,6))
ax1 = fig.add_axes([0.1,0.13,0.92,0.8])


xs, ys, zs = [], [], []

data = open("huge_nh_data.txt", "r")
for line in data:
	s = line.split()
	xs.append(float(s[0]))
	ys.append(float(s[1]))
	zs.append(float(s[2]))


grbs, grb_ra, grb_dec, grb_nh = [], [], [], []

grb_data_file = open("grb_data2010.txt", "r")
for line in grb_data_file:
		s = line.split()
		if "2009" not in line:
			grbs.append(s[3])
			grb_ra.append(float(s[9]))
			grb_dec.append(float(s[10]))
			grb_nh.append(float(s[11]))
grb_data_file.close()



ax1.set_xlim([0, 360])
ax1.set_ylim([-90, 90])
ax1.set_xlabel(r'$\rm{Right\, Ascension\, (J2000)}$', fontsize=18)
ax1.set_ylabel(r'$\rm{Declination\, (J2000)}$', fontsize=18)
ax1.xaxis.set_ticks(np.arange(0, 361, 60))
ax1.yaxis.set_ticks(np.arange(-80, 81, 20))

sc = ax1.scatter(xs, ys, c=zs, s = 20, cmap="Greys", marker = "o", edgecolor="none", vmax = 0.30e22)
cb = plt.colorbar(sc, pad=0.01, ticks=[0.5e21, 1.0e21, 1.5e21, 2.0e21, 2.5e21, 3.0e21])
cb.ax.set_yticklabels(["0.5", "1.0", "1.5", "2.0", "2.5", ">3"])
cb.set_label(r'$\rm{N_H\,10^{21}\,(cm^{-2})}$', fontsize=18)

grbs_plot = ax1.scatter(grb_ra, grb_dec, facecolors='none', edgecolors='red', s = 40)


fig.savefig("nh_map.pdf", format = "pdf")
fig.savefig("nh_map.jpeg", format = "jpeg")
#plt.show()
