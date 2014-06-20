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

from grb_astro_lib import Coordinates

detec = 0.
notdetec = 0.


data = open("nh_map/grb_data2007.txt", "r")
for line in data:
	s = line.split()
	ra =  s[9]
	dec = s[10]
	coords = Coordinates(ra, dec)
	l, b = coords.equatorial2galatic()
	
	color = s[12]
	if s[0] == "long":
		if color == "green":
			detec += 1.
			t = s[8].split(":")
			delay = ( (float(t[0])*3600) + (float(t[1])*60) + float(t[2]))
			#print delay
	
		if color == "red":
			notdetec += 1.
			t = s[8].split(":")
			delay2 = ( (float(t[0])*3600) + (float(t[1])*60) + float(t[2]))
			#print delay2


	
data.close()

print detec
print notdetec
print (detec/notdetec * 100)