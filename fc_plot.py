#!/usr/bin/env python

"""
Create finding charts
\tUsage: fc.py #options#

Example:fc_class.py -l MJD 56013.32594 -stars stars.txt -c 246.864417 -29.415000 -t 246.864417 -29.415000 0.0007 "GRB 120327A" -cmap PuBu arcsinh -s 0.02 -b GROND r

Options:
\t-f \t "fits file"\t \t \t \tDefault: "r.fits"
\t-t \t "ra" "dec" "circle radius" "label"\tDefault: "ra_cen" "dec_cen" "0.0007" "GRB"
\t-c \t "ra_cen" "dec_cen" \t \t \tDefault: "ra central pixel" "dec central pixel"
\t-s \t "FC size"\t \t \t \tDefault: "0.033"
\t-b \t "band"\t \t \t \t \tDefault: " " / Top-Right
\t-cmap \t "cmap style" "stretch factor"\t \tDefault: "Greys" "arsinh"
\t-format  "format" \t \t \t \tDefault: "pdf"
\t-l \t "Additional Label" \t \t \tDefault: " " / Bottom-Left
\t-stars \t "Starlist" \t \t \t \tlist containing the stars for calibration or error circles: ra dec radius label

\t \t \tExample I - calibration Stars:

\t \t \t75.054527 -3.347666 0.0007 I
\t \t \t75.064527 -3.357666 0.0007 II
\t \t \t75.074527 -3.367666 0.0007 III
\t \t \t75.084527 -3.377666 0.0007 IV

\t \t \tExample II - error Circles:

\t \t \t75.054527 -3.347666 0.005 XRT
\t \t \t75.064527 -3.357666 0.003 GROND
\t \t \t75.074527 -3.367666 0.02  BAT


for colormaps check: http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps/
stretch factors: arsinh, power, linear, sqrt
"""

import matplotlib
matplotlib.use('Agg')
import aplpy
import sys
import math
import scipy
import pylab
import matplotlib.pyplot as pyplot
import numpy as np
from argparse import ArgumentParser


print __doc__

parser = ArgumentParser()
parser.add_argument("-t", "--target", dest="target", nargs="*", default=[900., 900., 0.0007, "GRB"])
parser.add_argument("-f", "--file", dest="file", default="GROND_r_OB_ana.fits")
parser.add_argument("-l", "--label", dest="label", nargs="*", default=" ")
parser.add_argument("-c", "--center", dest="center", nargs=2, default=[900., 900], type=float)
parser.add_argument("-s", "--size", dest="size", default= 4./120., type=float)
parser.add_argument("-b", "--band", dest="band", nargs="*", default=" ")
parser.add_argument("-cmap", "--colormap", dest="cmap", nargs="*", default= ["Greys", "arcsinh"])
parser.add_argument("-format", "--output_format", dest="format", default= "pdf")
parser.add_argument("-stars", "--starlist", dest="stars", default=0)
args = parser.parse_args()


fits_file = args.file
target_parameters = args.target
target_radius = float(target_parameters[2])
target_id = target_parameters[3]
additional_label = args.label
add_label = ""
for i in range(0, len(additional_label), 1):
	add_label += (str(additional_label[i])+ " ")

central_coordinates = args.center
fc_size = args.size
scaley = fc_size*0.08
band = args.band
band_name = ""
for n in range(0, len(band), 1):
	band_name += (str(band[n])+ " ")
colormap = args.cmap
cmap_style = colormap[0] 
cmap_stretch = colormap[1]
output_format = args.format


# --------------------------------------------------------------------- #
# 			Get star-parameters from the given starist					#
starlist = args.stars

stars_ra = []
stars_dec = []
stars_rad = []
stars_label = []
if starlist != 0:
	stars = open(starlist, "r")
	for line in stars:
		s = line.split()
		stars_ra.append(float(s[0]))
		stars_dec.append(float(s[1]))
		stars_rad.append(float(s[2]))
		stars_label.append(s[3])
else:
	stars_ra.append(0)
	



fc = aplpy.FITSFigure(str(fits_file))

def coords(x, y):
	x_cen = 0
	y_cen = 0
	if x == 900.:
		x_conv, y_conv = fc.pixel2world(x, y)
		x_cen += x_conv
		y_cen += y_conv
	else:
		x_cen += x
		y_cen += y
	return x_cen, y_cen


ra_target, dec_target = coords(float(target_parameters[0]), float(target_parameters[1]))
ra_cen, dec_cen = coords(central_coordinates[0], central_coordinates[1])


def initCompass(self, parent):
    self._ax1 = parent._ax1
    self._wcs = parent._wcs
    self.world2pixel = parent.world2pixel
    self.pixel2world = parent.pixel2world
    self._initialize_compass()


aplpy.overlays.Compass.__init__ = initCompass

fc.set_system_latex(True)
fc.set_theme('publication')
fc.tick_labels.set_xformat(r'$hh:mm:ss$')
fc.tick_labels.set_yformat(r'$dd:mm$')
fc.axis_labels.set_xtext(r'$\rm{Right\, Ascension\, (J2000)}$')
fc.axis_labels.set_ytext(r'$\rm{Declination\, (J2000)}$')
fc.axis_labels.set_font(size=18)
fc.ticks.set_minor_frequency(10)
fc.ticks.set_length(15)
fc.add_label(0.88, 0.95, band_name, relative=True, fontsize=35, color="#225ea8")
fc.add_label(0.15, 0.05, add_label, relative=True, fontsize=35, color="#225ea8")
fc.show_circles(ra_target, dec_target, target_radius, color="#225ea8", linewidth=1.5)
fc.add_label(ra_target, dec_target-scaley, target_id, color="#225ea8", fontsize=22)
fc.show_colorscale(cmap=cmap_style, stretch=cmap_stretch)
fc.add_scalebar(0.01)
fc.frame.set_linewidth(3)
fc.scalebar.show(30/3600., label= r'$\rm{30"}$', corner = 'top left',  linewidth=3, color="#225ea8")
fc.scalebar.set_font(size=16)
fc.recenter(ra_cen, dec_cen, fc_size)
fc.compass = aplpy.overlays.Compass(fc)
fc.compass.show_compass(color='#225ea8', corner = 4, length=0.1) 
fc.compass._compass[0].set_arrowstyle('-')
fc.add_label(0.95, 0.17, "N", relative=True, fontsize=35, color="#225ea8")

if stars_ra[0] != 0:
    for i in range(0, len(stars_ra), 1):
        fc.show_circles(stars_ra[i], stars_dec[i], stars_rad[i], color="black", linewidth=1.5)
        fc.add_label(stars_ra[i], stars_dec[i]+scaley, stars_label[i], color="black", fontsize=22)
else:
    print "no Starlist given"

fc.save(target_id+"."+output_format, format=output_format)


