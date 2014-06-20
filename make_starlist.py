#!/usr/bin/env python

'''
create a STARLIST for gr_vss_lc.py from a given *ana.result file

Usage:
\t make_starlist.py <options>

\t \t -f    \t "resultfile" \t \t \t \t \tDefault: "*ana.result"
\t \t -d    \t "difference" \t \t \t \t \tDefault: "0.008"

'''
import sys
import math
import scipy
import numpy

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", dest="file", default="*ana.result")
parser.add_argument("-d", "--difference", dest="difference", default=0.008, type=float)

args = parser.parse_args()


resultfile = open(args.file, "r")

for line in resultfile:
	s = line.split()
		if "0.98" in line:
			if (math.abs(float(s[2]) - float(s[4]))) < args.difference:
				print line

resultfile.close()