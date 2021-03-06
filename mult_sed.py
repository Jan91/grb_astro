#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-red", "--red", dest="red", default="smc")
parser.add_argument("-zmin", "--zmin", dest="zmin", default=0.0, type=float)
parser.add_argument("-zmax", "--zmax", dest="zmax", default=12.0, type=float)
parser.add_argument("-avmax", "--avmax", dest="avmax", default=0.5, type=float)
parser.add_argument("-ebv", "--ebv", dest="ebv", default=0.00, type=float)
args = parser.parse_args()


seds = open("create_seds.sh", "w")

filter_files = ["GRB_g_relativeLC.txt",
"GRB_r_relativeLC.txt",
"GRB_i_relativeLC.txt",
"GRB_z_relativeLC.txt",
"GRB_J_relativeLC.txt",
"GRB_H_relativeLC.txt",
"GRB_K_relativeLC.txt"]

t, terr = [], []
g, gerr = [], []
r, rerr = [], []
i, ierr = [], []
z, zerr = [], []
J, Jerr = [], []
H, Herr = [], []
K, Kerr = [], []


for filt in filter_files:
	if "_g_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			t.append(float(s[0]))
			terr.append(float(s[1]))
			g.append(float(s[2]))
			gerr.append(float(s[3]))
		f.close()
	elif "_r_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			r.append(float(s[2]))
			rerr.append(float(s[3]))
		f.close()
	elif "_i_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			i.append(float(s[2]))
			ierr.append(float(s[3]))
		f.close()
	elif "_z_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			z.append(float(s[2]))
			zerr.append(float(s[3]))
		f.close()
	elif "_J_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			J.append(float(s[2]))
			Jerr.append(float(s[3]))
		f.close()
	elif "_H_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			H.append(float(s[2]))
			Herr.append(float(s[3]))
		f.close()
	elif "_K_" in filt:
		f = open(filt, "r")
		for line in f:
			s = line.split()
			K.append(float(s[2]))
			Kerr.append(float(s[3]))
		f.close()


for filt in range(1, len(t), 1):
	burst = ""
	if filt < 10:
		burst += "012000"+str(filt)
	else:
		burst += "01200"+str(filt)
	try:
		line = "".join(("grb_z.py -g ", str(g[filt])," ", str(gerr[filt]), 
		" -r ", str(r[filt])," ",str(rerr[filt]),
		" -i ", str(i[filt])," ",str(ierr[filt]),
		" -z ", str(z[filt])," ",str(zerr[filt]),
		" -j ", str(J[filt])," ",str(Jerr[filt]),
		" -h ", str(H[filt])," ",str(Herr[filt]),
		" -k ", str(K[filt])," ",str(Kerr[filt]),
		" -zmin ", str(args.zmin),
		" -zmax ", str(args.zmax),
		" -avmax ", str(args.avmax),
		" -ebv ", str(args.ebv),
		" -red ", args.red,
		" -burst ", burst))
		seds.write(line+"\n")
	except IndexError:
		print "line", filt, "has missing magnitudes"


















