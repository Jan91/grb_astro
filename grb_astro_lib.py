from math import *
import os
import time
import datetime
import numpy
import sys
import requests
import numpy as np
import matplotlib as plt
from pylab import *
from lxml import html
from mpl_toolkits.mplot3d import Axes3D


# ---------------------------------------------- #
# -----------  FLux & Magnitudes  -------------- #


def ab2flux(ab):
	'''
	Convert AB Magnitude to Flux
	'''
	return 10**(((23.9-float(ab))/2.5))



def flux2ab(flux):
	'''
	Convert Flux to AB Magnitude
	'''
	return -2.5*log10(float(flux))+23.9


def abfluxerr(ab, err):
	ab, err = float(ab), float(err)
	flux = ab2flux(ab)
	return [ab2flux(ab-err)-flux, flux-ab2flux(ab+err)]

def fluxaberr(flux, err):
	flux, err = float(flux), float(err)
	return (flux2ab(flux-err)-flux2ab(flux+err))/2.



# ---------------------------------------------- #
# --------------  Coordinates  ----------------- #

class Coordinates:

	def __init__(self, ra, dec):
		self.ra = ra
		self.dec = dec


	def addzero(val, n):
		if float(val) < 10:
			if n == 1:
				val = '0%.0f' %val
			elif n == 2:
				val = '0%.2f' %val
			elif n == 3:
				val = '0%.3f' %val
		else:
			if n == 1:
				val = '%.0f' %val
			elif n == 2:
				val = '%.2f' %val
			elif n == 3:
				val = '%.3f' %val
		return val

	def deg2sex(self):
		ra = float(self.ra)
		dec = float(self.dec)	
		hours = int(ra/15)
		minutes = int((ra/15.-hours)*60)
		seconds = float((((ra/15.-hours)*60)-minutes)*60)
		retra = '%s:%s:%s' %(addzero(hours,1), addzero(minutes,1), addzero(seconds,3))

		degree = int(dec)
		minutes = int((dec-degree)*60)
		seconds = float((((dec-degree)*60)-minutes)*60)
		if dec < 0:
			retdec = '-%s:%s:%s' %(addzero(-1*degree,1), addzero(-1*minutes,1), addzero(-1*seconds,2))
		else:
			retdec = '+%s:%s:%s' %(addzero(degree,1), addzero(minutes,1), addzero(seconds,2))

		return retra, retdec

	def sex2deg(self):

		ra = self.ra.split(':')
		dec = self.dec.split(':')

		retra = (float(ra[0])+float(ra[1])/60.+float(ra[2])/3600.)*15

		if float(dec[0]) <= 0:
			retdec = float(dec[0])-float(dec[1])/60.-float(dec[2])/3600.
		else:
			retdec = float(dec[0])+float(dec[1])/60.+float(dec[2])/3600.

		return float(retra), float(retdec)

	def equatorial2plane(self, ra0, dec0, scale):
		ra = radians(self.ra)
		dec = radians(self.dec)
		ra0 = radians(ra0)
		dec0 = radians(dec0)

		delta_ra = ra - ra0

		x = ((cos(dec) * sin(delta_ra)) /
			((sin(dec) * sin(dec0) +
			cos(dec) * cos(dec0) * cos(delta_ra))))
		y = (((sin(dec) * cos(dec0) -
			cos(dec) * sin(dec0) * cos(delta_ra))) /
			((sin(dec) * sin(dec0) +
			cos(dec) * cos(dec0) * cos(delta_ra))))

		x = x * (-1) * (180./pi) * 3600. / scale
		y = y * (180./pi) * 3600. / scale

		return x, y


	def equatorial2galatic(self):
		'''
		Convert equatorial to galatic Coordinates

		still not sure if i do this the right way
		'''
		ra = math.radians(float(self.ra))
		dec = math.radians(float(self.dec))

		alpha_G = math.radians(float(192.85))
		delta_G = math.radians(float(27.13))

		#alpha_G = math.radians(float(266.4))
		#delta_G = math.radians(float(-28.94))

		#b = math.degrees(asin(cos(dec) * cos(delta_G) * cos(ra-alpha_G)) + (sin(dec) * sin(delta_G)))
		#l = 122.9 - math.degrees(asin(cos(dec)*sin(ra-alpha_G)/cos(math.radians(b))))
		
		#return l, b
		l2 = 303 - math.degrees(math.atan(sin(alpha_G-ra)/( (cos(alpha_G-ra)*sin(delta_G)) - (tan(dec)*cos(delta_G)) ) ))
		b2 = math.degrees(math.asin( (sin(dec)*sin(delta_G)) + (cos(dec)*cos(delta_G)*cos(alpha_G-ra)) ))
		 
		#b3 = math.degrees(math.asin(cos(dec)*cos(delta_G)*cos(ra-alpha_G)+(sin(dec)*sin(delta_G))))

		#x = sin(dec)-(sin(b3)*sin(delta_G))
		#y = cos(dec)*cos(delta_G)*sin(ra-alpha_G)
		#l3 = 0
		#if (y>0) and (x>0):
		#	l3 += math.degrees(math.atan(x/y))
		#elif (y>0) and (x<0):
		#	l3 += math.degrees(math.atan(x/y))
		#elif (y<0) and (x<0):
		#	l3 += math.degrees(math.atan(x/y)) - 180
		#elif (y<0) and (x>0):
		#	l3 += math.degrees(math.atan(x/y)) - 380
		#else:
		#	l3 += 0

		return l2, b2

	def distto(self, ra2, dec2):
		ra = radians(float(self.ra))
		dec = radians(float(self.dec))
		ra2 = radians(float(ra2))
		dec2 = radians(float(dec2))

		dra = ra2-ra
		ddec = dec2-dec

		d1 = 2.*asin(sqrt((sin(ddec/2.))**2+cos(dec)*cos(dec2)*(sin(dra/2.))**2))
		d2 = d1 * 380./2./pi

		return d2

	def getebv(self):
		ra = str(self.ra)
		dec = str(self.dec)
		url = "http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=" + ra + "+" + dec
		page = requests.get(url)
		tree = html.fromstring(page.text)

		ebv_schlafly = tree.xpath('//meanvaluesandf/text()')
		ebv_schlegel = tree.xpath('//meanvaluesfd/text()')
	
		split_schlafly = ebv_schlafly[0].split()
		split_schlegel = ebv_schlegel[0].split()
	
		schlafly_corr = float(split_schlafly[0])
		schlegel_corr = float(split_schlegel[0])

		#print "Schlafly & Finkbeiner 2011; Mean Value:", schlafly_corr, "(mag)"
		#print "Schlegel et al. 1998; Mean Value:", schlegel_corr, "(mag)"

		return schlafly_corr, schlegel_corr

	def getnh(self):
		ra = str(self.ra)
		dec = str(self.dec)

		url = ("http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl?Entry=" + 
		ra + "+" + dec + 
		"&NR=GRB%2FSIMBAD%2FNED&CoordSys=Equatorial&equinox=2000&radius=1.00&usemap=0")

		page = requests.get(url)
		tree = html.fromstring(page.text)

		nh_list = tree.xpath('//b/text()')
		nh_string = nh_list[4].split()
		nh = float(nh_string[6])
	
		return nh


#data = open("nh_map/huge_nh_data.txt", "r")
#grbs = open("nh_map/grb_data2010.txt", "r")

#l, b , nh = [], [], []
#grb_l, grb_b, grb_nh, grb_c = [], [], [], []

#for line in data:
#	s = line.split()
#	coords = Coordinates(float(s[0]), float(s[1]))
#	ll, bb = coords.equatorial2galatic()
#	#l.append(float(ll))
#	#b.append(float(bb))
#	nh.append(float(s[2]))
#	l.append(float(s[0]))
#	b.append(float(s[1]))
#data.close()
#
#for line in grbs:
#	s = line.split()
#	coords = Coordinates(float(s[9]), float(s[10]))
#	ll, bb = coords.equatorial2galatic()
#	grb_l.append(float(ll))
#	grb_b.append(float(bb))
#	grb_nh.append(float(s[11]))
#	grb_c.append(str(s[12]))	
#grbs.close()
#
#
#
#fig = figure(figsize=(18,10))
#ax1 = fig.add_subplot(111, projection="hammer")
#ax1 = fig.add_axes([0.1,0.1,0.8,0.8])

#sc = ax1.scatter(l, b, c=nh, s = 6, cmap="PuBu", marker = "o", edgecolor="none", vmax=0.48e22)
#grb_sc = ax1.scatter(grb_l, grb_b, c=grb_c, s = 30, marker = "o", edgecolor="none")

#plt.show()
#fig.savefig("lbtest.jpeg", format = "jpeg")















