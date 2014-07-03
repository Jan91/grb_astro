#!/usr/bin/env python

import matplotlib
import matplotlib.ticker as ticker
import aplpy
import sys
import math
import scipy
from pylab import *
import matplotlib.pyplot as pyplot
import numpy as np
from argparse import ArgumentParser

def xconv(x):
    return x

def yconv(y):
    return (2.5*(29 - math.log(y, 10)) - 48.6)

def READ_LC_OUT(lc_file):

	x_val, x_err, y_val, y_err, res_val = [], [], [], [], [] 

	for line in lc_file:
		s = line.split()
		flux = 10**(29 - ((float(s[2])+48.6)/2.5))
		x_val.append(float(s[0]))
		x_err.append(float(s[1]))
		y_val.append(flux)
		ferr = (10**(29 - ((float(s[3])+float(s[2])+48.6))/2.5)) - flux
		#res_val.append((flux - MODEL/-ferr) # Include MODEL here
 		y_err.append(ferr)

 	return x_val, x_err, y_val, y_err, res_val


def GET_DATA():

	grond_lcs = ["GRB_g_relativeLC.txt",
	"GRB_r_relativeLC.txt",
	"GRB_i_relativeLC.txt",
	"GRB_z_relativeLC.txt",
	"GRB_J_relativeLC.txt",
	"GRB_H_relativeLC.txt",
	"GRB_K_relativeLC.txt"]

	g_data = []
	r_data = []
	i_data = []
	z_data = []
	j_data = []
	h_data = []
	k_data = []

	for filt in grond_lcs:
		if "_g_" in filt:
			f = open(filt, "r")
			g, gerr, gt, gterr, gres = READ_LC_OUT(f)
			g_data.append([g, gerr, gt, gterr, gres])
		elif "_r_" in filt:
			f = open(filt, "r")
			r, rerr, rt, rterr, rres = READ_LC_OUT(f)
			r_data.append([r, rerr, rt, rterr, rres])
		elif "_i_" in filt:
			f = open(filt, "r")
			i, ierr, it, iterr, ires = READ_LC_OUT(f)
			i_data.append([i, ierr, it, iterr, ires])		
		elif "_z_" in filt:
			f = open(filt, "r")
			z, zerr, zt, zterr, zres = READ_LC_OUT(f)
			z_data.append([z, zerr, zt, zterr, zres])
		elif "_J_" in filt:
			f = open(filt, "r")
			J, Jerr, Jt, Jterr, Jres = READ_LC_OUT(f)
			j_data.append([J, Jerr, Jt, Jterr, Jres])
		elif "_H_" in filt:
			f = open(filt, "r")
			H, Herr, Ht, Hterr, Hres = READ_LC_OUT(f)
			h_data.append([H, Herr, Ht, Hterr, Hres])
		elif "_K_" in filt:
			f = open(filt, "r")
			K, Kerr, Kt, Kterr, Kres = READ_LC_OUT(f)
			k_data.append([K, Kerr, Kt, Kterr, Kres])


	return g_data, r_data, i_data, z_data, j_data, h_data, k_data


def PLOT_LC():

	g, r, i, z, j, h, k = GET_DATA()

	fig = figure(figsize=(8,10))


	xmin = r[0][0][0] - (0.2 * r[0][0][0])
	xmax = r[0][0][len(r[0][0])-1] + (0.2 * r[0][0][len(r[0][0])-1])

	# RESIDUALS
	ax1 = fig.add_axes([0.12,0.1,0.78,0.1])

	xlabel(r"$\rm Time\, since\, BAT\, Trigger\, [s]$", fontsize=18)
	ylabel(r"$\rm Residuals\, [\sigma]$", fontsize=18)
	xlim([xmin, xmax]) 
	ax1.set_xscale('log')
	ax1.xaxis.tick_bottom()


	# GROND Light-curve

	ax2 = fig.add_axes([0.12,0.2,0.78,0.378])


	ylabel(r"$\rm Flux\, [\mu Jy]$", fontsize=18)
	ax2.yaxis.tick_left()
	ax2.set_xscale('log')
	ax2.set_yscale('log')
	xlim([xmin, xmax])


	gpl = plt.errorbar(g[0][0], g[0][2], xerr=g[0][1], yerr=g[0][3], color='#0c2c84', fmt='o', capsize=0, markersize=5)
	rpl = plt.errorbar(r[0][0], r[0][2], xerr=r[0][1], yerr=r[0][3], color='#225ea8', fmt='o', capsize=0, markersize=5)
	ipl = plt.errorbar(i[0][0], i[0][2], xerr=i[0][1], yerr=i[0][3], color='#1d91c0', fmt='o', capsize=0, markersize=5)
	zpl = plt.errorbar(z[0][0], z[0][2], xerr=z[0][1], yerr=z[0][3], color='#41b6c4', fmt='o', capsize=0, markersize=5)
	jpl = plt.errorbar(j[0][0], j[0][2], xerr=j[0][1], yerr=j[0][3], color='#7fcdbb', fmt='o', capsize=0, markersize=5)
	hpl = plt.errorbar(h[0][0], h[0][2], xerr=h[0][1], yerr=h[0][3], color='#c7e9b4', fmt='o', capsize=0, markersize=5)
	kpl = plt.errorbar(k[0][0], k[0][2], xerr=k[0][1], yerr=k[0][3], color='#c7e9b4', fmt='o', capsize=0, markersize=5)

	#ymin2 = yconv(ymin)
	#ymax2 = yconv(ymax)

	ax3 = twinx()
	xlabel('', fontsize=20)
	ylabel(r'$\rm AB\, Magnitude$', fontsize=18)
	ax3.set_xscale('log')
	ax3.set_xticklabels([])
	xlim([xmin, xmax])
	
	ay3 = twiny()
	ay3.xaxis.tick_top()
	ay3.set_xticklabels([])
	ay3.yaxis.tick_right()
	ay3.set_xscale('log')
	xlim([xmin, xmax])
	# XRT Light-Curve



	ax4 = fig.add_axes([0.12,0.575,0.78,0.378])

	ylabel(r"$\rm Flux_{0.3-10 keV}\, [10^{-13} erg\cdot s^{-1}cm^{-2}]$", fontsize=18)
	ax4.set_xscale('log')
	ax4.set_yscale('log')
	ax4.set_xticklabels([])
	xlim([xmin, xmax])

	fig.savefig("lightcurve.jpeg")
	#show()


if __name__ == "__main__":
	#print __doc__
	PLOT_LC()










