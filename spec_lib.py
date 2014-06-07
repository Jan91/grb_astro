# -*- coding: utf-8 -*-
import numpy, scipy
#from scipy import interpolate, integrate, sqrt

from sedmodels import Avlaws, lyalpha
from astro import abflux

def grond_dict():

  grond_dict = {}
  grond_dict["Bands"] = ["g", "r", "i", "z", "J", "H", "K"]
  grond_dict["Wavelengths"] = {}
  grond_dict["Wavelengths"]["g"] = 4587
  grond_dict["Wavelengths"]["g_err_N"] = 730.6
  grond_dict["Wavelengths"]["g_err_P"] = 760.8
  grond_dict["Wavelengths"]["r"] = 6219.8
  grond_dict["Wavelengths"]["r_err_N"] = 844.9
  grond_dict["Wavelengths"]["r_err_P"] = 847.3
  grond_dict["Wavelengths"]["i"] = 7640.7
  grond_dict["Wavelengths"]["i_err_N"] = 485.8
  grond_dict["Wavelengths"]["i_err_P"] = 515.9
  grond_dict["Wavelengths"]["z"] = 8989.6
  grond_dict["Wavelengths"]["z_err_N"] = 735.6
  grond_dict["Wavelengths"]["z_err_P"] = 538.8
  grond_dict["Wavelengths"]["J"] = 12399.2
  grond_dict["Wavelengths"]["J_err_N"] = 1206.6
  grond_dict["Wavelengths"]["J_err_P"] = 1159
  grond_dict["Wavelengths"]["H"] = 16468.4    
  grond_dict["Wavelengths"]["H_err_N"] = 1366
  grond_dict["Wavelengths"]["H_err_P"] = 1335.5
  grond_dict["Wavelengths"]["K"] = 21705.5
  grond_dict["Wavelengths"]["K_err_N"] = 1561
  grond_dict["Wavelengths"]["K_err_P"] = 1470.6
  grond_dict["Energies"] = {}
 
  h = 6.63e-34
  c = 3.0e8
  keV = 1.6e-16


  for band in grond_dict["Bands"]:

    wav_nm = grond_dict["Wavelengths"][band] * 1e-10
    wav_N_nm = wav_nm - (grond_dict["Wavelengths"]["%s_err_N" % band] * 1e-10)
    wav_P_nm = wav_nm + (grond_dict["Wavelengths"]["%s_err_P" % band] * 1e-10)

#    print "Wavelengths"
#    print wav_nm, wav_N_nm, wav_P_nm
#    print ""

    en = ((h*c)/wav_nm) / (keV)
    en_N = ((h*c)/wav_N_nm) / (keV)
    en_P = ((h*c)/wav_P_nm) / (keV)

    grond_dict["Energies"][band] = en
    grond_dict["Energies"]["%s_err_N" % band] = abs(en-en_N)
    grond_dict["Energies"]["%s_err_P" % band] = abs(en-en_P)
    
  return grond_dict

def idl_hist(array, array_axis, binsize):

  array_x = numpy.array([])
  array_y = numpy.array([])

  for i in range(len(array)):
    array_x = numpy.append(array_x, array_axis[i])
    array_x = numpy.append(array_x, array_axis[i] + binsize)

    array_y = numpy.append(array_y, array[i])
    array_y = numpy.append(array_y, array[i])

  # Ensure it has tails
  array_x = numpy.append(array_x[0], array_x)
  array_x = numpy.append(array_x, array_x[len(array_x)-1])

  array_y = numpy.append(0, array_y)
  array_y = numpy.append(array_y, 0)

  return(array_x, array_y)


def bindata_digitize(x,y,binsize, xerr, yerr):

  x = numpy.array([x[i] for i in range(len(xerr)) if (xerr[i]>0 and yerr[i]>0)])
  y = numpy.array([y[i] for i in range(len(xerr)) if (xerr[i]>0 and yerr[i]>0)])
  xerr = numpy.array([xerr[i] for i in range(len(xerr)) if (xerr[i]>0 and yerr[i]>0)])
  yerr = numpy.array([yerr[i] for i in range(len(xerr)) if (yerr[i]>0)])

  bins = numpy.arange(x[0], x[-1]+binsize, binsize)
  digitized = numpy.digitize(x, bins)

  x_binned = numpy.array([x[digitized == i].mean() for i in range(1, len(bins))])
  y_binned = numpy.array([y[digitized == i].mean() for i in range(1, len(bins))])

  x_err_b, y_err_b = [], []

  print "YERR:", yerr

  for i in range(1,len(bins)):

    temp_x_err = xerr[digitized == i]
    temp_y_err = yerr[digitized == i]

    print "Temp err", temp_y_err

    if len(temp_y_err)<1:
      y_err_b.append(0)
      continue

    val = 0
    for err in temp_y_err:

      if err<0:
        err = 0

      print "Val1:", val
      val = scipy.sqrt(val**2 + err**2)
      print "Val2:", val

    val /= len(temp_y_err)

    print val
#    if i == 3:
#      sys.exit(0)
	


    y_err_b.append(val)

  print y_err_b

  x_err_b = numpy.array(x_err_b)
  y_err_b = numpy.array(y_err_b)

  print y_err_b

  return x_binned, y_binned, x_err_b, y_err_b


def bindata(x,y,binsize, xerr=[], yerr=[]):

  newArrayAvg = []
  newArrayAvgErr = []
  newArrayMed = []
  xArray = []
  xArrayErr = []
  tmpArray = []

  left = x[0]
  right = x[0] + binsize

  val = 0
  val_err = 0
  x_err = 0

  for i in range(len(y)):
    if x[i] >= left and x[i] < right:
      tmpArray.append(y[i])
      if len(yerr)>0 and len(xerr)>0:
        val = (val + y[i])/2.


      if len(yerr)>0:
        val_err = scipy.sqrt(val_err**2 + yerr[i]**2) / 2.0
 
      if len(xerr)>0: 
        x_err = scipy.sqrt(x_err**2 + xerr[i]**2) / 2.0

    else:

      xArray.append((left+right)/2.)
      xArrayErr.append(binsize/2.0)

      left = x[i]
      right = x[i] + binsize
      newArrayAvg.append(val)
      newArrayAvgErr.append(val_err)
      newArrayMed.append(numpy.median(tmpArray))
      tmpArray = []
      val = 0
      val_err = 0
      x_err = 0

  xArray = numpy.array(xArray)
  xArrayErr = numpy.array(xArrayErr)
  newArrayAvg = numpy.array(newArrayAvg)
  newArrayAvgErr = numpy.array(newArrayAvgErr)
  newArrayMed = numpy.array(newArrayMed)
  return (xArray, xArrayErr, newArrayAvg, newArrayMed, newArrayAvgErr)

def absll():
    return {'Lyg_972': [972.536, 0],
            'CIII_977': [977.020, 0],
            'SiII_989': [989.873, 0],
            'NII_989': [989.799, 0],
            'SiII_1012': [1012.502, 0],
           'SiII_1020': [1020.6989, 0],
            'MgII_1026': [1026.1134, 0],
            'Lyb_1025': [1025.728, 0],
            'ArI_1066': [1066.660, 0],
            'FeII_1081': [1081.8748, 0],
            'NII_1083': [1083.990, 0],
            'FeII_1125': [1125.4478, 0],
            'SiII_1190': [1190.254, 0],
            'SiII_1193': [1193.2897, 0],
            'SiII_1206': [1206.500, 0],
            'Lya_1215': [1215.67, 73],
             'NV_1238': [1238.82, 0.14],
             'NV_1242': [1242.80, 0.07],
             'SII_1250': [1250.58, 0.15],
             'SII_1253': [1253.52, 0.24],
             'SII_1259': [1259.52, 0.00],
             'SiII_1260': [1260.42, 1.26],
             'SiII^*_1264': [1264.74, 0.66],
             'CI_1277': [1277.25, 0.09],
             'OI_1302': [1302.17, 0.00],
             'SiII_1304': [1304.37, 2.29],
             'OI^*_1304': [1304.86, 0.00],
             'SiII^*_1309': [1309.28, 0.27],
            'NiII_1317': [1317.22, 0.11],
            'CI_1328': [1328.83, 0.08],
            'CII_1334': [1334.53, 1.73],
            'CII^*_1335': [1335.71, 0.0],
            'ClI_1347': [1347.24, 0.20],
            'NiII_1370': [1370.13, 0.13],
            'SiIV_1393': [1393.76, 0.95],
            'SiIV_1402': [1402.77, 0.68],
            'GaII_1414': [1414.40, 0.16],
            'NiII_1415': [1415.72, 0.16],
            'CO_1419': [1419.0, 0.06],
            'NiII_1454': [1454.84, 0.08],
            'ZnI_1457': [1457.57, 0.08],
            'CO_1477': [1477.5, 0.21],
            'SiII_1526': [1526.71, 0.93],
            'SiII^*_1533': [1533.43, 0.42],
            'CoII_1539': [1539.47, 0.03],
            'CIV_1548': [1548.20, 2.18],
            'CIV_1550': [1550.77, 2.18],
            'CI_1560': [1560.31, 0.09],
            'FeII^*_1570': [1570.25, 0.08],
            'FeII^*_1602': [1602.49, 0.08],
            'FeII_1608': [1608.45, 0.85],
            'FeII_1611': [1611.20, 0.18],
            'FeII*_1618': [1618.47, 0.03],
            'FeII*_1621': [1621.69, 0.11],
            'FeII*_1629': [1629.16, 0.1],
            'FeII*_1631': [1631.13, 0.1],                        
            'FeII*_1634': [1634.35, 0.5],
            'FeII*_1636': [1636.33, 0.5],
            'FeII*_1639': [1639.40, 0.5],                        
            'CI_1656': [1656.93, 0.14],
            'AlII_1670': [1670.79, 1.04],
            'SiI_1693': [1693.29, 0.07],#1693.03 0.07 Â± 0.02
            'NiII_1703': [1703.41, 0.14],#        1702.48        0.14 Â± 0.02        
            'NiII_1709': [1709.60, 0.08],#1709.45        0.08 Â± 0.02        
            'NiII_1741': [1741.55, 0.14],#1742.02        0.14 Â± 0.02        
            'MgI_1747': [1747.79, 0.05],#1748.02        0.05 Â± 0.01        b
            'NiII_1751': [1751.92, 0.09],#        1752.06        0.09 Â± 0.02        
            'NiII_1804': [1804.47, 0.03],#1805.31        0.03 Â± 0.01        a
            'SI_1807': [1807.31, 0.29],#        1808.14        0.29 Â± 0.02        a
            'SiII_1808': [1808.01, 0.29],#                        a
            'SiII*_1816': [1816.93, 0.12],#        1816.95        0.12 Â± 0.01        a, b
            'SiII^*_1817': [1817.45, 0.12],#                        a, b
            'MgI_1827': [1827.93, 0.09],#1827.87        0.09 Â± 0.02        
            'NiII_1842': [1842.89, 0.12],#        1844.27        0.12 Â± 0.02        b
            'AlIII_1854': [1854.72, 0.89],#        1854.83        0.89 Â± 0.02        
            'AlIII_1862': [1862.79, 0.68],#        1863.56        0.68 Â± 0.02        
            'FeI_1875': [1875.16, 0.10],#1874.46        0.10 Â± 0.01        a
            'FeI_1883': [1883.78, 0.17],#1882.67        0.17 Â± 0.02        a
            'SiIII_1892': [1892.03, 0.10],#        1890.58        0.10 Â± 0.02        b
            'FeI_1901': [1901.77, 0.14],#1901.36        0.14 Â± 0.02        a
            'TiII_1905': [1905.77, 0.11],#1905.99        0.11 Â± 0.02        a
            'TiII_1906': [1906.24, 0.11],#                
            'TiII_1910': [1910.75, 0.14],#1910.56        0.14 Â± 0.02        a
            'FeI_1937': [1937.27, 0.09],#1938.04        0.09 Â± 0.02        a
            'CoII_1941': [1941.29, 0.07],#        1941.97        0.07 Â± 0.01        a, b
            'CoII_2012': [2012.17, 0.06],#2010.99        0.06 Â± 0.02        b
            'CrII_2017': [2017.57, 0.08],#2017.61        0.08 Â± 0.02        b
            'ZnII_2026': [2026.14, 0.60],#                
            'CrII_2026': [2026.27, 0.60],#2025.97        0.60 Â± 0.02        a
            'MgI_2026': [2026.48, 0.00],#                
            'CrII_2056': [2056.26, 0.19],#        2055.51        0.19 Â± 0.02        
            'CrII_2062': [2062.23, 0.53],#2062.67        0.53 Â± 0.02        a
            'ZnII_2062': [2062.66, 0.53],#                a
            'CrII_2066': [2066.16, 0.12],#2066.29        0.12 Â± 0.02        a
            'NiII*_2166': [2166.23, 0.26],#        2167.65        0.26 Â± 0.02        a
            'FeI_2167': [2167.45, 0.00],#                
            'NiII^*_2175': [2175.22, 0.07],#        2175.83        0.07 Â± 0.01        
            'MnI_2185': [2185.59, 0.51],#2186.96        0.51 Â± 0.02        
            'NiII^*_2217': [2217.2, 0.24],#2217.02        0.24 Â± 0.01        
            'NiII_2223': [2223.0, 0.07],#2224.98        0.07 Â± 0.01        
            'FeII_2249': [2249.88, 0.31],#        2250.02        0.31 Â± 0.02        
            'FeII_2260': [2260.78, 0.38],#2261.00        0.38 Â± 0.02        
            'NII_2316': [2316.7, 0.19],#2316.36        0.19 Â± 0.02        
            'FeII^{**}_2328': [2328.11, 0.10],#        2327.65        0.10 Â± 0.02        
            'FeII^*_2333': [2333.52, 0.34],#        2333.30        0.34 Â± 0.02        
            'FeII_2344': [2344.21, 1.74],#2346.01        1.74 Â± 0.02        a
            'FeII^*_2345': [2345.00, 0.00],#                        
            'FeII^*_2349': [2349.02, 0.26],#        2349.60        0.26 Â± 0.02        a
            'FeII^*_2359': [2359.83, 0.28],#        2360.01        0.28 Â± 0.01        
            'FeII^*_2365': [2365.55, 0.18],#        2365.40        0.18 Â± 0.01        
            'FeII_2374': [2374.46, 1.00],#2374.46        1.00 Â± 0.02        
            'FeII^*_2381': [2381.49, 0.00],#                        a
            'FeII_2382': [2382.77, 1.65],#2383.59        1.65 Â± 0.02        a
            'FeII^*_2383': [2383.79, 0.00],#                        a
            'FeII^*_2389': [2389.36, 0.18],#        2389.40        0.18 Â± 0.02        a
            'FeII^*_23961': [2396.15, 0.00],#                        a
            'FeII^*_23963': [2396.36, 0.71],#        2396.50        0.71 Â± 0.02        a
            'FeII^*_2399': [2399.98, 0.00],#                        a
            'FeII^*_2405': [2405.16, 0.40],#        2406.04        0.40 Â± 0.02        a
            'FeII^*_2407': [2407.39, 0.00],#                        a
            'FeII^*_24112': [2411.25, 0.00],#                        a
            'FeII^*_24118': [2411.80, 0.51],#        2411.68        0.51 Â± 0.03        a
            'FeII^*_2414': [2414.05, 0.00],#                        a
            'ScII_2561': [2561.00, 0.18],#2562.52        0.18 Â± 0.01        a
            'ScII_2563': [2563.97, 0.00],#                a
            'MnII_2576': [2576.88, 0.45],#2577.04        0.45 Â± 0.02        
            'FeII_2586': [2586.65, 1.33],#2586.49        1.33 Â± 0.02        a
            'MnII_2594': [2594.50, 0.45],#2594.36        0.45 Â± 0.02        
            'FeII*_2599': [2599.15, 0.00],#        2600.23        1.85 Â± 0.03        a
            'FeII_2600': [2600.17, 1.85],#                        a
            'MnII_2606': [2606.46, 0.56],#2607.53        0.56 Â± 0.01        a
            'FeII^*_2607': [2607.87, 0.00],#                        a
            'FeII^*_2612': [2612.65, 0.51],#        2613.11        0.51 Â± 0.01        a
            'FeII^*_2614': [2614.61, 0.00],#                        
            'FeII^*_2618': [2618.40, 0.06],#        2618.51        0.06 Â± 0.01        a
            'FeII^*_2621': [2621.19, 0.00],#                        a
            'FeII^*_2622': [2622.45, 0.00],#                        a
            'FeII^*_2626': [2626.45, 0.00],#        2629.82        0.90 Â± 0.02        a
            'FeII^*_2629': [2629.08, 0.90],#                        a
            'FeII^*_2631': [2631.83, 0.00],#                        a
            'FeII^*_2632': [2632.11, 0.00],#                        a
            'FeII^*_2740': [2740.4, 0.07],#739.45        0.07 Â± 0.01        b
            'FeII^*_2747': [2747.9, 0.16],#749.50        0.16 Â± 0.01        b
            'FeII^*_2756': [2756.28, 0.08],#        2756.50        0.08 Â± 0.01        b
            'MgII_2796': [2796.35, 1.71],#2796.21        1.71 Â± 0.02        
            'MgII_2803': [2803.53, 1.47],#2803.50        1.47 Â± 0.02        
            'MgI_2852': [2852.96, 0.78],#        2852.97        0.78 Â± 0.01        
            'FeI_2967': [2967.77, 0.05],#        2966.24        0.05 Â± 0.01        b
            'FeI_2984': [2984.44, 0.05]        ,#2983.49        0.05 Â± 0.01        
            'TiII_3073': [3073.88, 0.08],#3076.01        0.08 Â± 0.01        
            'TiII_3242': [3242.93, 0.11],#3239.56        0.11 Â± 0.01        
            'TiII_3384': [3384.74, 0.15],#3385.85        0.15 Â± 0.02        b
            'CaII_3934': [3934.78, 0.76],#3933.97        0.76 Â± 0.02        
            'CaII_3969': [3969.59, 0.66],#3969.98        0.66 Â± 0.02        
            'CaI_4227': [4227.92, 0.11],#        4226.93        0.11 Â± 0.02        b
            'MgH_5209': [5209.45, 0.09] }#5210.75        0.09 Â± 0.02        b

def plotspectra(spectratxt):
	
	# Load file
	spectraf = open(spectratxt, "r")
	spectrad = spectraf.readlines()
	spectraf.close()

	# Parse text
	spectradata_x = numpy.array([])
	spectradata_y = numpy.array([])

	for i in range(len(spectrad)):

		tmp = spectrad[i].replace(".","").split(" ")
		tmpnorm = spectrad[i].replace("\n","").split(" ")
		tmpnorm = numpy.array([i for i in tmpnorm if i !=""])
		if tmp[0].isdigit():
	
			spectradata_x = numpy.append(spectradata_x, float(tmpnorm[0]))
			spectradata_y = numpy.append(spectradata_y, float(tmpnorm[1]))

	return spectradata_x, spectradata_y

def parse_filters(inf="filters.txt"):

  filtf = open(inf, "r")
  filtl = filtf.readlines()
  filtf.close()

  grdict = {}
  for key in ["X", "g", "r", "i", "z", "J", "H", "K"]:
    grdict[key] = []

  for line in filtl:
    line_parse = [i for i in line.replace("\n","").split("\t") if i!=""]
    
    if line_parse[0] == "A":	continue
    
    grdict["X"].append(float(line_parse[0]))
    grdict["g"].append(float(line_parse[1]))
    grdict["r"].append(float(line_parse[2]))
    grdict["i"].append(float(line_parse[3]))
    grdict["z"].append(float(line_parse[4]))
    grdict["J"].append(float(line_parse[5]))
    grdict["H"].append(float(line_parse[6]))
    grdict["K"].append(float(line_parse[7]))
    
    
  for key in grdict:
     grdict[key] = numpy.array(grdict[key])
     
  return grdict

def power_law_kev(x, f0, b):
    f = []
    for i in x:
      f.append(f0 * ((i)**-b))

    return f

def broken_power_law_kev(x,f0,b1,b2,eb):
    f = []
    for i in x:
      if i <= eb:
          f.append(f0*(i**(-b1)))
      else:
          f.append(f0*((eb**((b2-b1)))*(i**-(b2))))

    return f

def ext(y, law='smc', av=0.2148, z=0.9840):
    # y in nm
    #av = 0.2148
    #z = 0.9840
    
    smc=Avlaws(y/(1+z), law)
    ext = scipy.exp(-1/1.086*av*smc)
    
    return ext

def xrt_2_mujy(x):
  ev_2_J = 1.60217656535E-19
  ev_2_hz = 2.417989348E14
  cm2_2_m2 = 1.0e4
  J_2_erg = 1.0e7

  flux_J_cm2_s_hz = (ev_2_J / ev_2_hz) * x

  flux_w_m2_hz = flux_J_cm2_s_hz * cm2_2_m2

  flux_jy = 1e26 * flux_w_m2_hz
  flux_ujy = 1e6 * flux_jy
  
  return flux_ujy
 
def xrt_2_erg(flux_kev_cm2_s_kev, wavelength_angstrom):
  ev_2_J = 1.60217656535E-19
  ev_2_hz = 2.417989348E14
  cm2_2_m2 = 1.0e4
  J_2_erg = 1.0e7

  flux_J_cm2_s_hz = (ev_2_J / ev_2_hz) * flux_kev_cm2_s_kev

  flux_w_m2_hz = flux_J_cm2_s_hz * cm2_2_m2

  flux_jy = 1e26 * flux_w_m2_hz
  flux_ujy = 1e6 * flux_jy

  flux_erg_s_cm2_Hz = flux_J_cm2_s_hz * J_2_erg

  hz_2_m = 3.335640951E-9
  hz_2_ang = (1e10)/hz_2_m

  flux_erg_s_cm2__A = flux_erg_s_cm2_Hz * hz_2_ang

  flux_erg_s_cm2_A = flux_erg_s_cm2__A / (wavelength_angstrom**2)
   
  return flux_erg_s_cm2_A
 
def nm_2_engy(x):
  
  # give wavelength in nm
  c = 2.99792458E8
  h = 6.62606957E-34
  keV = 1.60217656535E-16

  ### 
  wavelength_m = x*1e-9

  ###
  frequency_hz = c / wavelength_m

  ###
  energy_ev = frequency_hz * h
  energy_kev = energy_ev / keV
  
  return energy_kev

def engy_2_nm(energy_keV):

  # give wavelength in keV
  c = 2.99792458E8
  h = 6.62606957E-34
  keV = 1.60217656535E-16

  ###
  energy_ev = energy_keV * keV

  ###
  frequency_hz = energy_ev / h

  ###
  wavelength_m = c / frequency_hz
 
  ###
  wavelength_nm = wavelength_m / 1e-9
 
  return wavelength_nm

def lyman_alpha(opten, optmo, z):
   da, db, dc, dd = lyalpha(z)
   i = 0

   for energy in opten:

       if (energy*(1+z) > 10.14E-3) and (energy*(1+z)) < 0.1:
           optmo[i] = optmo[i]*da

       if (energy*(1+z) > 11.9E-3) and (energy*(1+z)) < 0.1:
           optmo[i] = optmo[i]*db

       if (energy*(1+z) > 12.59E-3) and (energy*(1+z)) < 0.1:
           optmo[i] = optmo[i]*dc

       if (energy*(1+z) > 12.98E-3) and (energy*(1+z)) < 0.1:
           optmo[i] = optmo[i]*dd

       if (energy*(1+z) > 13.6E-3) and (energy*(1+z)) < 0.1:
           optmo[i] = optmo[i]*1E-30
       i+=1
   return optmo


def rmlybr(model, energy, lim=1E-3):
   newmodel, newenergy = [], []
   for value, enval in zip(model,energy):
      if value > lim:
         newmodel.append(value)
         newenergy.append(enval)
   return numpy.array(newmodel), numpy.array(newenergy)

def sigma_E(energy_keV):
   """ Original fortran code:
      REAL EMAX(14)
C          (edge energies where polynomial changes -- in keV)
      REAL C0(14)
C          (zero-order polynomial coefficients for 14 energy intervals)
      REAL C1(14)
C          (1st-order polynomial coefficients)
      REAL C2(14)
C          (2nd-order polynomial coefficients)
C
      data emax/.100,.284,.400,.532,.707,.867,1.303,
     #    1.840,2.471,3.210,4.038,7.111,8.331,10.0/
      data c0/17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,
     #    202.7,342.7,352.2,433.9,629.0,701.2/
      data c1/608.1,267.9,18.8,66.8,145.8,-380.6,169.3,
     #    146.8,104.7,18.7,18.7,-2.4,30.9,25.2/
      data c2/-2150.,-476.1,4.3,-51.4,-61.1,294.0,-47.7,
     #    -31.5,-17.0,0.0,0.0,0.75,0.0,0.0/
C
C     Start:
C
      E=energy/1.e3
C          (convert to keV)
      do 100 i=1,14
        if (E .lt. Emax(i)) goto 200
100   continue
      i=14
200   sigism=(c0(i)+c1(i)*E+c2(i)*E*E)/E**3 * 1.E-24
      return
      end
   """
   # Input is keV
   emax = numpy.array([.100,.284,.400,.532,.707,.867,1.303,1.840,2.471,3.210,4.038,7.111,8.331,10.0])
   c0 = numpy.array([17.3,34.6,78.1,71.4,95.5,308.9,120.6,141.3,202.7,342.7,352.2,433.9,629.0,701.2])
   c1 = numpy.array([608.1,267.9,18.8,66.8,145.8,-380.6,169.3,146.8,104.7,18.7,18.7,-2.4,30.9,25.2])
   c2 = numpy.array([-2150.,-476.1,4.3,-51.4,-61.1,294.0,-47.7,-31.5,-17.0,0.0,0.0,0.75,0.0,0.0])
   sig_ism = 0.0

   for i in range(14):
      if energy_keV < emax[i] and energy_keV > 0.03:
	sig_ism = (c0[i]+c1[i]*energy_keV+c2[i]*energy_keV*energy_keV)/energy_keV**3 * 1.0e-24

        break    
      else:
        continue

   return sig_ism

def phabs(energy_keV, n_H):

   # n_H: 10^22 atoms/cm^2
   # energy_keV: keV

   sig_ism = sigma_E(energy_keV)
   M_E = scipy.exp(-(n_H*1e22) * sig_ism)
   return M_E

def highecut(energy, e_cut, e_folding):

  if energy >= e_cut:
    return scipy.exp((e_cut-energy)/e_folding)
  else:
    return 1.0
    
def parse_opt(inf="grb121217Aobs.sed"):

	f = open(inf, "r")
	l = f.readlines()
	f.close()

	OptDict = {\
		"Wavelength": [], \
		"Wavelength_Err_P": [], \
		"Wavelength_Err_N": [], \
		"Magnitude": [], \
		"Magnitude_Err": [], \
	}

	for i in range(len(l)):

		line = [j for j in l[i].strip().split(" ") if j != ""]

		OptDict["Wavelength"].append(float(line[0]))
		OptDict["Wavelength_Err_P"].append(float(line[1]))
		OptDict["Wavelength_Err_N"].append(abs(float(line[2])))
		OptDict["Magnitude"].append(float(line[3]))
		OptDict["Magnitude_Err"].append(float(line[4]))

	for item in OptDict:
		OptDict[item] = numpy.array(OptDict[item])

	return OptDict

def parse_xrt(inf="grb121217Axrt.spec", ignore=False):

        f = open(inf, "r")
        l = f.readlines()
        f.close()

        XRTDict = {\
                "Counts": [], \
                "Counts_Err": [], \
                "Energy": [], \
                "Energy_Err": [], \
                "Model": [], \
        }

        for i in range(len(l)):

                line = [j for j in l[i].strip().split(" ") if j != ""]

		if line[0][0] in ["R", "@", "!"]: continue
		if "NO" in line: continue

		XRTDict["Energy"].append(float(line[0]))
		XRTDict["Energy_Err"].append(float(line[1]))
		XRTDict["Counts"].append(float(line[2]))
		XRTDict["Counts_Err"].append(float(line[3]))

		XRTDict["Model"].append(float(line[4]))

        for item in XRTDict:
                XRTDict[item] = numpy.array(XRTDict[item])

	if ignore:
		indx = []
		for ignore_l, ignore_r in ignore:
			indx_l = XRTDict["Energy"] > ignore_l

			indx_r = XRTDict["Energy"] < ignore_r

			indx = indx_l-indx_r
			for item in XRTDict:
				XRTDict[item] = XRTDict[item][indx]

        return XRTDict


def band_function_kev(E, E_c, norm, alpha, beta):
  
  """  Approximation to the Gamma-Ray Burst model by D. Band, et al.,
   1993 (ApJ 413, 281) for XSPEC programs:
   N(E) = n  * (E/E_p)**alfa * exp(-E/tem), if E < (alfa-beta)tem
        = n' * (E/E_p)**beta              , if E > (alfa-beta)tem"""

  E_b = (alpha-beta)*E_c


  if E < E_b:
    f = ( (E/100.)**(alpha) ) * scipy.exp(-1.0*E/E_c)
  elif E >= E_b:
    f = ( ( (alpha-beta)*(E_c/100.) )**(alpha-beta) ) * ( (E/100.)**(beta) ) * ( scipy.exp(-1.0*(alpha-beta)) )

  f *= (norm)

  return f

def parseModel(model):

  """You need to give the XRT spectrum class. It will then extract all the models and their best fit parameters. This is useful for such cases as:
  
  modelDict["zdust"]["E(B-V)"]["Value"]

  but you can do it yourself by using the below code as well.
  """

  modelDict = {}

  for component in model.componentNames:
    comp = getattr(model, component)
    modelDict[component] = {}

    print ""
    print "Component name: %s" % component

    for parameter in comp.parameterNames:

      print "paramater: %s" % parameter

      par = getattr(comp, parameter)
      value = par.values

      print "values: %s" % value
      print ""
      print ""

      modelDict[component][parameter] = {}
      if len(value) == 1:
        modelDict[component][parameter]["Value"] = value
        modelDict[component][parameter]["Delta"] = 0.0
        modelDict[component][parameter]["min"] = 0.0
        modelDict[component][parameter]["bot"] = 0.0
        modelDict[component][parameter]["top"] = 0.0
        modelDict[component][parameter]["max"] = 0.0
      else:
        modelDict[component][parameter]["Value"] = value[0]
        modelDict[component][parameter]["Delta"] = value[1]
        modelDict[component][parameter]["min"] = value[2]
        modelDict[component][parameter]["bot"] = value[3]
        modelDict[component][parameter]["top"] = value[4]
        modelDict[component][parameter]["max"] = value[5]

  return modelDict

def twobrpl(t, par):
  # Parameter 0: Power law slope 1
  # Parameter 1: Breaktime 1
  # Paramter 2: Smoothness 1
  # Paramter 3: Power law slope 2
  # Paramter 4: Breaktime 2
  # Paramter 5: Smoothness 2
  # Paramter 6: Power law slope 3
  # Paramter 7: Power law slope 4
  # Parameter 8: Relative Normalisation

  brkpl1 = ((t/par["Time_1"])**(-par["smooth_1"]*par["alpha_1"])+(t/par["Time_1"])**(-par["smooth_1"]*par["alpha_2"]))**(-1/par["smooth_1"])
  brkpl2 = par["Relative_Normalisation"]*((t/par["Time_2"])**(-par["smooth_2"]*par["alpha_3"])+(t/par["Time_2"])**(-par["smooth_2"]*par["alpha_4"]))**(-1/par["smooth_2"])

  return (brkpl1+brkpl2)
