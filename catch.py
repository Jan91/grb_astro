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


class Coordinates:

	def __init__(self, ra, dec):
		self.ra = ra
		self.dec = dec

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

		return (x, y)

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
		url = "http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=" + ra + dec
		page = requests.get(url)
		tree = html.fromstring(page.text)

		ebv_schlafly = tree.xpath('//meanvaluesandf/text()')
		ebv_schlegel = tree.xpath('//meanvaluesfd/text()')
	
		split_schlafly = ebv_schlafly[0].split()
		split_schlegel = ebv_schlegel[0].split()
	
		schlafly_corr = float(split_schlafly[0])
		schlegel_corr = float(split_schlegel[0])

		return schlafly_corr

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



for i in np.arange(202.0, 361.0, 1.0):
	for j in np.arange(0.0, 91.0, 1.0):
		j = "+" + str(j)
		i = str(i)
		coords = Coordinates(i, j)
		ebv = coords.getebv()
		print i, j, ebv








