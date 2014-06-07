#!/usr/bin/env python
'''
Galactic Foregound Reddening for given coordinates
\t Usage:
\t \t get_nh.py <options>

\t Example:
\t \t get_ebv -c 41.17900 -26.15310

\t Options:
\t \t -c "ra_j2000" "dec_j2000" 
'''

import urllib
import requests
from lxml import html
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-c", "--coords", dest="coordinates", nargs=2)
args = parser.parse_args()
ra = args.coordinates[0] 
dec = args.coordinates[1]


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
	getnh(ra, dec)