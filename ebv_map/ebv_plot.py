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

data = open("ebv_data.txt", "r")
for line in data:
	s = line.split()
	xs.append(float(s[0]))
	ys.append(float(s[1]))
	zs.append(float(s[2]))
data.close()

grbs, grb_ra, grb_dec, grb_ebv, grb_col = [], [], [], [], []

grb_data_file = open("grb_data2010.txt", "r")
for line in grb_data_file:
		s = line.split()
		grbs.append(s[3])
		grb_ra.append(float(s[9]))
		grb_dec.append(float(s[10]))
		grb_ebv.append(float(s[14]))
		if s[12] == "green":
			grb_col.append("#74c476")
		elif s[12] == "red":
			grb_col.append("#d7301f")
		elif s[12] == "yellow":
			grb_col.append("#df65b0")
		elif s[12] == "blue":
			grb_col.append("#08306b")

grb_data_file.close()



ax1.set_xlim([0, 360])
ax1.set_ylim([-90, 90])
ax1.set_xlabel(r'$\rm{Right\, Ascension\, (J2000)}$', fontsize=18)
ax1.set_ylabel(r'$\rm{Declination\, (J2000)}$', fontsize=18)
ax1.xaxis.set_ticks(np.arange(0, 361, 60))
ax1.yaxis.set_ticks(np.arange(-80, 81, 20))

sc = ax1.scatter(xs, ys, c=zs, s = 20, cmap="PuBu", marker = "o", edgecolor="none", vmax = 0.5, vmin= 0.0)
cb = plt.colorbar(sc, pad=0.01, ticks=[0.08, 0.16, 0.24, 0.32, 0.4, 0.48])
cb.ax.set_yticklabels(["0.08", "0.16", "0.24", "0.32", "0.40", ">0.5"])
#grbs_plot = ax1.scatter(grb_ra, grb_dec, edgecolor=grb_col, facecolor=grb_ebv, s = 50, vmax = 1.0, vmin= 0.0, cmap="PuBu")
grbs_plot = ax1.scatter(grb_ra, grb_dec, c=grb_ebv, s = 60, cmap="PuBu", edgecolor=grb_col, vmax = 0.5, vmin = 0.0, marker = "o")



cb.set_label(r'$\rm{E_{B-V}\, (mag)}$', fontsize=18)


#grbs_plot = ax1.scatter(grb_ra, grb_dec, marker = "o", color="#31a354")


fig.savefig("ebv_map.pdf", format = "pdf")
fig.savefig("ebv_map.jpeg", format = "jpeg")
#plt.show()
