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
ax1 = fig.add_axes([0.1,0.13,0.95,0.8])


xs = []
ys = []
zs = []

#data = open("huge_nh_data.txt", "r")
data = open("ebv_data.txt", "r")
for line in data:
	s = line.split()
	if float(s[0]) < 361:
		xs.append(float(s[0]))
		ys.append(float(s[1]))
		zs.append(float(s[2]))

ax1.set_xlim([0, 360])
ax1.set_ylim([-90, 90])
ax1.set_xlabel(r'$\rm{Right\, Ascension\, (J2000)}$', fontsize=18)
ax1.set_ylabel(r'$\rm{Declination\, (J2000)}$', fontsize=18)
ax1.xaxis.set_ticks(np.arange(0, 361, 60))
ax1.yaxis.set_ticks(np.arange(-80, 81, 20))

#sc = ax1.scatter(xs, ys, c=zs, s = 20, cmap="PuBu", marker = "o", edgecolor="none")
sc = ax1.scatter(xs, ys, c=zs, s = 20, cmap="PuBu", marker = "o", edgecolor="none", vmax = 2.0)
cb = plt.colorbar(sc, pad=0.01)
#cb.set_label(r'$\rm{N_H\, (cm^{-2})}$', fontsize=18)
cb.set_label(r'$\rm{E_{B-V}\, (mag)}$', fontsize=18)

grbs = ["100423A",
"100425A",
"100615A",
"100621A",
"100724A",
"100728B",
"100814A",
"101023A",
"110223B",
"110305A",
"111107A",
"111008A",
"111229A",
"120119A",
"120311A",
"120701A",
"120728A",
"120815A",
"120909A",
"120922A",
"121024A",
"130408A",
"130514A",
"130615A",
"130816A"]

grb_ra = [136.45919,	
299.19659,	
177.205583,	
315.30466,	
194.543333,	
44.056254,	
22.47309,	
317.963492,	
150.233250,	
260.880417,	
129.47771,	
60.45104,	
76.28686,	
120.02895,	
273.09235,	
80.348000,	
137.094250,	
273.95762,	
275.73633,	
234.74850,	
70.47208,	
134.405417,	
296.28306,	
274.82962,	
197.14079]

grb_dec =[
+21.49833,
-26.43069,
-19.481111,
-51.10624,
-11.102583,
+0.28105,
-17.99541,
-65.387669,
-68.301667,
-15.802444,
-66.52011,
-32.70916,
-84.71073,
-9.08162,
+14.29619,
-58.549528,
-54.437361,
-52.13126,
-59.44836,
-20.18172,
-12.29069,
-32.360806,
-7.976220,
-68.16125,
-58.94497]

grbs_plot = ax1.scatter(grb_ra, grb_dec, marker = "o", color="#31a354")



plt.show()
