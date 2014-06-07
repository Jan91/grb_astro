#!/usr/bin/env python
'''
Galactic Foregound Reddening for given coordinates
\t Usage:
\t \t get_corr.py <options>

\t Example:
\t \t get_ebv -c 41.17900 -26.15310

\t Options:
\t \t -c "ra_j2000" "dec_j2000" \t \t Default: 41.17900 -26.15310
'''

import requests
from lxml import html
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-c", "--coords", dest="coordinates", nargs=2, default=[str(41.17900), str(-26.15310)])
args = parser.parse_args()
coords = args.coordinates[0] + args.coordinates[1]
ra = args.coordinates[0] 
dec = args.coordinates[1]

def getebv(coords):

	url = "http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=" + coords
	page = requests.get(url)
	tree = html.fromstring(page.text)
	
	ebv_schlafly = tree.xpath('//meanvaluesandf/text()')
	ebv_schlegel = tree.xpath('//meanvaluesfd/text()')
	
	split_schlafly = ebv_schlafly[0].split()
	split_schlegel = ebv_schlegel[0].split()
	
	schlafly_corr = float(split_schlafly[0])
	schlegel_corr = float(split_schlegel[0])

	print "Schlafly & Finkbeiner 2011; Mean Value:", schlafly_corr, "(mag)"
	print "Schlegel et al. 1998; Mean Value:", schlegel_corr, "(mag)"

	return schlafly_corr, schlegel_corr

def getnh(ra, dec):

	url = ("http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl?Entry=" + 
		ra + "+" + dec + 
		"&NR=GRB%2FSIMBAD%2FNED&CoordSys=Equatorial&equinox=2000&radius=1.00&usemap=0")
	page = requests.get(url)
	tree = html.fromstring(page.text)

	nh_list = tree.xpath('//b/text()')
	nh_string = nh_list[4].split()
	nh = float(nh_string[6])
	
	print "Average nH:", nh, "(cm**-2)"
	return nh

if __name__ == "__main__":
	print __doc__
	getebv(coords)
	getnh(ra, dec)
