#!/usr/bin/env python

'''
Broadband SED-Plotting

Usage:
\t bb_sed.py <options>

\t \t -m    \t "model" \t \t \t \t \tDefault: "pow"
\t \t -p1   \t "slope 1" \t \t \t \t \tDefault: 0.5
\t \t -p2   \t "slope 2" \t \t \t \t \tDefault: 1.0
\t \t -norm \t "norm" \t \t \t \t \tDefault: 0.005
\t \t -be   \t "break_energy" \t \t \t \tDefault: 0.05
\t \t -f    \t "data_file" \t \t \t \t \tDefault: "sed.txt"
\t \t -nh   \t "neutral hydrogen density / milkyway" \t \tDefault: 0.2
\t \t -nhh  \t "neutral hydrogen density / host_galaxy" \tDefault: 0.2
\t \t -av   \t "galactic foreground extinction/reddening" \tDefault: 0.0
\t \t -avh  \t "host dust extinction" \t \t \tDefault: 0.2
\t \t -red  \t "host_extinction_law" \t \t \t \tDefault: "smc"
\t \t -z    \t "redshift" \t \t \t \t \tDefault: 1.0
\t \t -c    \t "plot color" \t \t \t \t \tDefault: "blue"
\t \t -n    \t "name" \t \t \t \t \tDefault: "grb_sed"
'''

import sys
import math
import scipy
import numpy

import matplotlib as plt
import matplotlib.pyplot as pyplot
from pylab import *

from spec_lib import *
import numpy.ma as ma
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-m", "--model", dest="model", default="pow")
parser.add_argument("-p1", "--slope1", dest="slope1", default=0.5, type=float)
parser.add_argument("-p2", "--sleop2", dest="slope2", default=1.0, type=float)
parser.add_argument("-norm", "--norm", dest="norm", default=0.005, type=float)
parser.add_argument("-be", "--break_energy", dest="break_energy", default=0.05, type=float)
parser.add_argument("-f", "--file", dest="file", default="sed.txt")
parser.add_argument("-nh", "--nh", dest="nh", default=0.2, type=float)
parser.add_argument("-nhh", "--nh_host", dest="nh_host", default=0.2, type=float)
parser.add_argument("-av", "--av", dest="av", default=0.0, type=float)
parser.add_argument("-avh", "--av_host", dest="av_host", default=0.2, type=float)
parser.add_argument("-red", "--reddenig_law", dest="host_extinction_law", default="smc")
parser.add_argument("-z", "--redshift", dest="redshift_z", default=1.0, type=float)
parser.add_argument("-c", "--color", dest="color", default='black')
parser.add_argument("-n", "--name", dest="name", default="grb_sed")
args = parser.parse_args()

class BroadBandSED:

	def PowerLaw(self, energy, p1, norm):
		'''
		Using a single Powerlaw
		'''
		f = []
		for i in energy:
			f.append(norm * ((i)**-p1))
		return f

	def BrokenPowerlaw(self, energy, p1, p2, break_energy, norm):
		'''
		Using a Broken-Powerlaw
		'''
		f = []
		for i in energy:
			if i <= break_energy:
				f.append(norm*(i**(-p1)))
			else:
				f.append(norm*((break_energy**((p2-p1)))*(i**-(p2))))
		return f

	def get_unfolded_model(self, model, energy, norm, p1, p2, break_energy):
		unfolded_model = []
		if model == "pow":
			bb = BroadBandSED()
			unfolded_model.append(bb.PowerLaw(energy, p1, norm))
		elif model == "bknpow":
			bb = BroadBandSED()
			unfolded_model.append(bb.BrokenPowerlaw(energy, p1, p2, break_energy, norm))
		else:
			print " NOTE: bb_sed.py didn't recognize the model, please choose between pow and bknpow \n \n"
		return unfolded_model

	def Model(self, energy, z, nh, nh_host, av, av_host, host_extinction_law):
		wavelength = engy_2_nm(energy)

		phabs_Gal = numpy.array([phabs(i, nh) for i in energy])
		phabs_Host = numpy.array([phabs(i*(1+z), nh_host) for i in energy])
		redd_Gal = numpy.array([ext(i, "mw", av, 0) for i in wavelength])
		redd_Host = numpy.array([ext(i, host_extinction_law, av_host, z) for i in wavelength])

		total_absorption = phabs_Gal*phabs_Host*redd_Gal*redd_Host
		total_absorption = lyman_alpha(energy, total_absorption, z)
		return total_absorption

	def Data(self, file):
		'''
		Reading data out of the given file
		'''
		h = 4.135667516e-15 #eVs
		ev = 1.6021176e-19 	#J
		kevmiJy = ev*h*1e36
		x = []
		y = []
		x_err = []  
		y_err = []
		data_file = open(file, "r")
		for line in data_file:
			s = line.split()
			x.append(float(s[0]))
			x_err.append(float(s[1]))
			y.append(float(s[2])*kevmiJy)
			y_err.append(float(s[3])*kevmiJy)
		return x, y, x_err, y_err

def xconv(x):
	'''
	convert Energy [keV] into Wavelenght [nm]
	'''
	h = 4.135667516e-15 #eVs
	c = 299792458.0
	return ((h*c)/x)*1e+6
def yconv(y):
	'''
	convert Flux [\mu Jy] into AB Magnitude [mag]
	'''
	return (2.5*(29 - math.log(y, 10)) - 48.6)

def SEDPlot():
	'''
	Plotting Broadband-SED with matplotlib
	'''

	h = 4.135667516e-15 #eVs
	ev = 1.6021176e-19 	#J
	kevmiJy = ev*h*1e36

	optical_energy = np.arange(0.0001, 0.0034, 0.0001)
	xray_energy = np.arange(0.31, 10.01, 0.01)
	energy = np.concatenate([optical_energy, xray_energy])
	mask_start = len(optical_energy)
	
	bb = BroadBandSED()
	unfolded_model = bb.get_unfolded_model(args.model, energy, args.norm*kevmiJy, args.slope1, args.slope2, args.break_energy)
	folded_model = unfolded_model*bb.Model(energy, args.redshift_z, args.nh, args.nh_host, args.av, args.av_host, args.host_extinction_law)
	x, y, x_err, y_err = bb.Data(args.file)
	
	mc = ma.array(folded_model[0])
	mc[mask_start] = ma.masked

	fig = figure(figsize=(8,6))
	ax1 = fig.add_axes([0.13,0.13,0.76,0.74])
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_xlim([2e-4, 10.1])
	ax1.set_xlabel(r"$Energy\, [keV]$", fontsize=18)
	ax1.set_ylabel(r"$Flux\, [\mu Jy]$", fontsize=18)
	ax1.errorbar(energy, unfolded_model[0], color=args.color, linestyle="dashed")
	ax1.errorbar(energy, mc, color=args.color)
	ax1.errorbar(x, y, xerr=x_err, yerr=y_err, color=args.color, fmt='o', capsize=0, markersize=3)

	y_min_left, y_max_left = ax1.get_ylim()
	x_min_bottom, x_max_bottom = ax1.get_xlim()

	ax2 = ax1.twinx()
	ax2.set_ylim(yconv(y_min_left), yconv(y_max_left))
	ax2.set_ylabel(r"$AB\, Magnitude\, [mag]$", fontsize=18)

	ax3 = ax1.twiny()
	ax3.set_xlim(xconv(x_min_bottom), xconv(x_max_bottom))
	ax3.set_xscale('log')
	ax3.set_xlabel(r"$Wavelenght_{\lambda}\, [nm]$", fontsize=18)

	print SEDPlot.__doc__
	fig.savefig(args.name+".pdf", format="pdf")
	fig.savefig(args.name+".jpeg", format="jpeg")

if __name__ == "__main__":
	print __doc__
	SEDPlot()




















