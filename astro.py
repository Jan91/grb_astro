# -*- coding: utf-8 -*-
"""
Miscellaneous functions for astronomical use

"""
__author__ = "Thomas Krühler <kruehler@mpe.mpg.de>"

__version__ = "0.1.1"

from math import *
import os
import time
import datetime
import numpy
import sys
from scipy import optimize
from pylab import array, sqrt
from xml.dom import minidom
from urllib2 import urlopen
from urllib import urlretrieve
from socket import setdefaulttimeout
from matplotlib.patches import Polygon


timeout = 90
setdefaulttimeout(timeout)

#####################################################

def airtovac(wl):
    """
    This corrects wavelengthes from air to vacuum system
    Wavelegnth in Angstrom
    Based in idlastro library
    """
    if isinstance(wl, list):
        wlv = []
        for swl in wl:
            if swl > 2000:
            # Convert to wavenumber squared
                sigma2 = (1.E4/swl)**2
                fact = 1. + 5.792105E-2/(238.0185-sigma2) + 1.67917E-3/(57.362-sigma2)
                wlv.append(swl*fact)
    elif isinstance(wl, float):
        sigma2 = (1.E4/wl)**2
        fact = 1. + 5.792105E-2/(238.0185-sigma2) + 1.67917E-3/(57.362-sigma2)
        wlv = wl*fact
        
    return wlv
#####################################################

def LDMP(z, H0=71, WM=0.27, WV=0.73):

    c = 299792.458 # velocity of light in km/sec
    Tyr = 977.8    # coefficent for converting 1/H into Gyr
    DTT = 0.5      # time from z to now in units of 1/H0
    DTT_Gyr = 0.0  # value of DTT in Gyr
    age = 0.5      # age of Universe in units of 1/H0
    age_Gyr = 0.0  # value of age in Gyr
    zage = 0.1     # age of Universe at redshift z in units of 1/H0
    zage_Gyr = 0.0 # value of zage in Gyr
    DCMR = 0.0     # comoving radial distance in units of c/H0
    DCMR_Mpc = 0.0 
    DCMR_Gyr = 0.0
    DA = 0.0       # angular size distance
    DA_Mpc = 0.0
    DA_Gyr = 0.0
    kpc_DA = 0.0
    DL = 0.0       # luminosity distance
    DL_Mpc = 0.0
    DL_Gyr = 0.0   # DL in units of billions of light years
    V_Gpc = 0.0
    a = 1.0        # 1/(1+z), the scale factor of the Universe
    z = float(z)


    

    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    age = 0.
    n=5000         # number of points in integrals
    for i in range(n):
        a = az*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot

    zage = az*age/n
    zage_Gyr = (Tyr/H0)*zage
    DTT = 0.0
    DCMR = 0.0

# do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    age = DTT+zage
    age_Gyr = age*(Tyr/H0)
    DTT_Gyr = (Tyr/H0)*DTT
    DCMR_Gyr = (Tyr/H0)*DCMR
    DCMR_Mpc = (c/H0)*DCMR

# tangential comoving distance

    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5*(exp(x)-exp(-x))/x 
        else:
            ratio = sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    DCMT = ratio*DCMR
    DA = az*DCMT
    DA_Mpc = (c/H0)*DA
    kpc_DA = DA_Mpc/206.264806
    DA_Gyr = (Tyr/H0)*DA
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    DL_Gyr = (Tyr/H0)*DL

# comoving volume computation

    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio = (0.125*(exp(2.*x)-exp(-2.*x))-x/2.)/(x*x*x/3.)
        else:
            ratio = (x/2. - sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/5. + (2./105.)*y*y
    VCM = ratio*DCMR*DCMR*DCMR/3.
    V_Gpc = 4.*pi*((0.001*c/H0)**3)*VCM

    if False:
        print 'For H_0 = ' + '%1.1f' % H0 + ', Omega_M = ' + '%1.2f' % WM + ', Omega_vac = ',
        print '%1.2f' % WV + ', z = ' + '%1.3f' % z
        print 'It is now ' + '%1.3f' % age_Gyr + ' Gyr since the Big Bang.'
        print 'The age at redshift z was ' + '%1.3f' % zage_Gyr + ' Gyr.'
        print 'The light travel time was ' + '%1.3f' % DTT_Gyr + ' Gyr.'
        print 'The comoving radial distance, which goes into Hubbles law, is',
        print '%1.1f' % DCMR_Mpc + ' Mpc or ' + '%1.1f' % DCMR_Gyr + ' Gly.'
        print 'The comoving volume within redshift z is ' + '%1.1f' % V_Gpc + ' Gpc^3.'
        print 'The angular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc or',
        print '%1.1f' % DA_Gyr + ' Gly.'
        print 'This gives a scale of ' + '%.2f' % kpc_DA + ' kpc/".'
        print 'The luminosity distance D_L is ' + '%1.1f' % DL_Mpc + ' Mpc or ' + '%1.1f' % DL_Gyr + ' Gly.'
        print 'The distance modulus, m-M, is '+'%1.2f' % (5*log10(DL_Mpc*1e6)-5)
     
    return(DL_Mpc)

##################################################################

def plotlayout(keyword):
    if keyword == 'pub':
       fs, labs, ps1, ps2, lw1, lw2 = 32, 24, 9, 6, 2.0, 1.5
    else:
       fs, labs, ps1, ps2, lw1, lw2 = 26, 22, 6, 3, 1.5, 1.0
    return fs, labs, ps1, ps2, lw1, lw2

############################################
def writecmd(filename, clobber = 0):
    if clobber == 0:
        fileopen = 'a'
    else:
        fileopen = 'w'
    now = datetime.datetime.now()
    f = open(filename, fileopen)   
    f.write(now.strftime("%Y-%m-%d %H:%M:%S"))
    f.write('\t')
    for arg in sys.argv:
        f.write(arg+' ')
    f.write('\n')
    f.close()
    
############################################

def isnumber(num):
    try:
        float(num)
        return True
    except:
        return False

###########################################

def abflux(ab):
    return 10**(((23.9-float(ab))/2.5))

###########################################

def fluxab(flux):
    return -2.5*log10(float(flux))+23.9

###########################################

def abfluxerr(ab, err):
    ab, err = float(ab), float(err)
    flux = abflux(ab)
    return [abflux(ab-err)-flux, flux-abflux(ab+err)]

###########################################

def fluxaberr(flux, err):
    flux, err = float(flux), float(err)
    return (fluxab(flux-err)-fluxab(flux+err))/2.
    
###########################################

def addzero(val, n):
    if isnumber(val):
        if float(val) < 10:
            if n == 1:
                val = '0%.0f' %val
            if n == 2:
                val = '0%.2f' %val
            if n == 3:
                val = '0%.3f' %val
        else:
            if n == 1:
                val = '%.0f' %val
            if n == 2:
                val = '%.2f' %val
            if n == 3:
                val = '%.3f' %val
    return val

###########################################

def deg2sexa(ra, dec):
    if isnumber(ra):
        ra = float(ra)
        hours = int(ra/15)
        minu = int((ra/15.-hours)*60)
        seco = float((((ra/15.-hours)*60)-minu)*60)
        retra = '%s:%s:%s' %(addzero(hours,1), addzero(minu,1), addzero(seco,3))
    else:
        retra = ra

    if isnumber(dec):
        dec = float(dec)
        degree = int(dec)
        minutes = int((dec-degree)*60)
        seconds = float((((dec-degree)*60)-minutes)*60)
        if dec < 0:
            retdec = '-%s:%s:%s' %(addzero(-1*degree,1), addzero(-1*minutes,1), addzero(-1*seconds,2))
        else:
            retdec = '+%s:%s:%s' %(addzero(degree,1), addzero(minutes,1), addzero(seconds,2))
    else:
        retdec = dec

    return retra, retdec

###########################################

def sexa2deg(ra, dec):

    if isnumber(ra):
        retra = ra
    else:
        ra = ra.split(':')
        retra = (float(ra[0])+float(ra[1])/60.+float(ra[2])/3600.)*15

    if isnumber(dec):
        retdec = dec
    else:
        dec = dec.split(':')
        if float(dec[0]) <= 0:
            retdec = float(dec[0])-float(dec[1])/60.-float(dec[2])/3600.
        else:
            retdec = float(dec[0])+float(dec[1])/60.+float(dec[2])/3600.
    return float(retra), float(retdec)

###########################################

def fill_between(ax, x, y1, y2, **kwargs):
    # add x,y2 in reverse order for proper polygon filling
    verts = zip(x,y1) + [(x[i], y2[i]) for i in range(len(x)-1,-1,-1)]
    poly = Polygon(verts, **kwargs)
    ax.add_patch(poly)
    ax.autoscale_view()
    return poly

###########################################

def equatorial_to_plane(ra, dec, ra0, dec0, scale):
    ra = radians(ra)
    dec = radians(dec)
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

###########################################

def mjdut(mjd):
    jd = float(mjd)+2400000.5
    jd=jd+0.5
    Z=int(jd)
    F=jd-Z
    alpha=int((Z-1867216.25)/36524.25)
    A=Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B-122.1)/365.25)
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    dd = B - D - int(30.6001*E) + F

    if E<13.5:
        mm=E-1
    if E>13.5:
        mm=E-13
    if mm>2.5:
        yyyy=C-4716
    if mm<2.5:
        yyyy=C-4715

    months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]

    h=int((dd-int(dd))*24)
    minute=int((((dd-int(dd))*24)-h)*60)
    sec=86400*(dd-int(dd))-h*3600-minute*60
    if h < 10:
        h= '0'+str(h)
    if minute < 10:
        minute = '0'+str(minute)
    if sec < 10:
        sec = '0%.2f'%sec
    else:
        sec = '%.2f'%sec

    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2

    daysum=0
    for y in range(mm-1):
        daysum=daysum+days[y]
    daysum=daysum+dd-1

    if days[1]==29:
        fracyear=yyyy+daysum/366.
    else:
        fracyear=yyyy+daysum/365.
    gregdate = months[mm-1]+" %i, %i, " % (dd, yyyy)+'%s'%h+':%s'%minute+':%s'%sec+' UTC'
    return fracyear, gregdate

############################################

def dist(ra1, dec1, ra2, dec2):
    ra1 = radians(float(ra1))
    dec1 = radians(float(dec1))
    ra2 = radians(float(ra2))
    dec2 = radians(float(dec2))
    dra = ra2-ra1
    ddec = dec2-dec1
    d1 = 2.*asin(sqrt((sin(ddec/2.))**2+cos(dec1)*cos(dec2)*(sin(dra/2.))**2))
    d2 = d1 * 380./2./pi
    return d2

############################################

def linfit2(x, y, yerr):
    x = array(x)
    y = array(y)
    yerr = array(yerr)

    fitfunc = lambda para, xvar: para[0] + para[1]*xvar + para[2]*xvar**2
    errfunc = lambda para, xvar, yvar, errvar: (yvar - fitfunc(para, xvar)) / errvar
    pinit = [0.0, 1.0, 1.0]
    out, covar, inf, mesg, ier = optimize.leastsq(errfunc, pinit, args=(x, y, yerr), full_output=1)

    med = out[0]
    mult1= out[1]
    mult2= out[2]

    chi2 = 0
    for (xval,yval,yvalerr) in zip(x,y,yerr):
        chi2+=((med+mult1*xval+mult2*xval**2-yval)/yvalerr)**2

    dof = len(y)-1
    redchi2 = chi2/dof
    sd1 = sqrt(covar[0][0])*sqrt(redchi2)
    sd2 = sqrt(covar[1][1])*sqrt(redchi2)
    sd3 = sqrt(covar[2][2])*sqrt(redchi2)

    return med, mult1,mult2, sd1, sd2, sd3, chi2, dof, redchi2

###########################################

def linfit1(x, y, yerr):
    x = array(x)
    y = array(y)
    yerr = array(yerr)

    fitfunc = lambda para, xvar: para[0] + para[1]*xvar
    errfunc = lambda para, xvar, yvar, errvar: (yvar - fitfunc(para, xvar)) / errvar
    pinit = [0.0, 1.0]
    out, covar, inf, mesg, ier = optimize.leastsq(errfunc, pinit, args=(x, y, yerr), full_output=1)

    med = out[0]
    mult = out[1]

    chi2 = 0
    for (xval,yval,yvalerr) in zip(x,y,yerr):
        chi2+=((med+mult*xval-yval)/yvalerr)**2

    dof = len(y)-1
    redchi2 = chi2/dof
    sd = sqrt(covar[0][0])*sqrt(redchi2)

    return med, mult, sd, chi2, dof, redchi2

############################################

def linfit(x, y, yerr):
    x = array(x)
    y = array(y)
    yerr = array(yerr)

    fitfunc = lambda para, xvar: para[0] + xvar
    errfunc = lambda para, xvar, yvar, errvar: (yvar - fitfunc(para, xvar)) / errvar
    pinit = [0.0]
    out, covar, inf, mesg, ier = optimize.leastsq(errfunc, pinit, args=(x, y, yerr), full_output=1)

    med = out

    chi2 = 0
    for (xval,yval,yvalerr) in zip(x,y,yerr):
        chi2+=((med+xval-yval)/yvalerr)**2

    dof = len(y)-1
    redchi2 = chi2/dof
    sd = sqrt(covar[0][0])*sqrt(redchi2)

    return med, sd, chi2, dof, redchi2

###########################################

def getebv(ra, dec):
    if dec > 0:
        dec = '++'+str(dec)
    else:
        dec = '+'+str(dec)
    ra = str(ra)
    urlstring = 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr='+ra+'+'+dec
    ebvall = minidom.parse(urlopen(urlstring))
    ebv = float(str(ebvall.getElementsByTagName('meanValue')[0].firstChild.data).split()[0])
    std = float(str(ebvall.getElementsByTagName('std')[0].firstChild.data).split()[0])
    ref = float(str(ebvall.getElementsByTagName('refPixelValue')[0].firstChild.data).split()[0])
    av = ebv*3.08
    return ebv, std, ref, av

###########################################

def getfc(ra, dec, radius, band):
    if dec > 0:
        dec = '++'+str(dec)
    else:
        dec = '+'+str(dec)
    ra = str(ra)
    radius = str(radius)
    band = str(band)
    if band in 'ugrizUBVRI':
        survey = 'survey=sdss&survey=dss'
    elif band in 'JHK':
        survey = 'survey=2mass'
    urlstring = 'http://irsa.ipac.caltech.edu/cgi-bin/FinderChart/nph-finder?locstr='+ra+dec+'&mode=prog&subsetsize='+radius+'&'+survey
    fc = minidom.parse(urlopen(urlstring))
    fitsurls = fc.getElementsByTagName('fitsurl')
    for fitsurl in fitsurls:
        chart = str(fitsurl.firstChild.data)
        head, tail = os.path.split(chart)
        # DSS2
        if 'poss2' in tail:
            cat = 'dss2'
            if band in ('ugUBV'):
                if 'blue' in tail:
                    downurl = chart
            if band in ('rR'):
                if 'red' in tail:
                    downurl = chart
            if band in ('izI'):
                if 'ir' in tail:
                    downurl = chart
        #SDSS
        if 'fpC' in tail:
            cat = 'sdss'
            if band in ('u'):
                if '-u' in tail:
                    downurl = chart
            if band in ('g'):
                if '-g' in tail:
                    downurl = chart
            if band in ('r'):
                if '-r' in tail:
                    downurl = chart
            if band in ('i'):
                if '-i' in tail:
                    downurl = chart
            if band in ('z'):
                if '-z' in tail:
                    downurl = chart
        #2MASS
        if '2mass' in tail:
            cat = '2mass'
            if band in ('J'):
                if '-j' in tail:
                    downurl = chart
            if band in ('H'):
                if '-h' in tail:
                    downurl = chart
            if band in ('K'):
                if '-k' in tail:
                    downurl = chart
    image = os.path.join(os.getcwd(),cat+'_'+ra+'_'+dec+'_'+band+'.fits')
    urlretrieve(downurl, image)

    return image

###########################################

def beta_ox(fluxR, fluxX, wl=625.):

    h = 4.135667E-15
    c = 2.99792458E8
    EX = 3.0E3 #X-ray energy in eV
    nuX = EX/h
    nuR = c/(wl*1E-9)

    beta_ox = log10(fluxR/fluxX)/log10(nuX/nuR)
    return beta_ox

###########################################

def readx(xfilename):

   ev = 1.602176E-19 #J
   h = 4.135667E-15 #eVs
   kevmuJy = ev*h*1E36

   xx, xy, xxerr, xyerr, xmodel = [],[],[],[],[]
   fin = file(xfilename)
   lines = [line.strip().split() for line in fin.readlines() if not line.strip().startswith('!')]
   fin.close()
   for line in lines:
      if isnumber(line[0]):
               xx.append(float(line[0]))
               xxerr.append(float(line[1]))
               keV=float(line[2])
               keVe=float(line[3])
               model = float(line[4])
               r = 1.
               flux = keV*kevmuJy/r
               fluxerr = keVe*kevmuJy/r
               modelflux = model*kevmuJy/r
               xy.append(flux)
               xyerr.append(fluxerr)
               xmodel.append(modelflux)
      elif line[0] == 'NO':
         break
      else:
         pass

   return array(xx), array(xy), array(xxerr), array(xyerr), array(xmodel)

###########################################

def readxm(xfilename):

   ev = 1.602176E-19 #J
   h = 4.135667E-15 #eVs
   kevmuJy = ev*h*1E36

   xx, xmodel = [],[]
   fin = file(xfilename)
   lines = [line.strip().split() for line in fin.readlines() if not line.strip().startswith('!')]
   fin.close()
   for line in lines:
      if isnumber(line[0]):
               xx.append(float(line[0]))
               model = float(line[4])
               r = 1.
               modelflux = model*kevmuJy/r
               xmodel.append(modelflux)
      elif line[0] == 'NO':
         break
      else:
         pass

   return array(xx), array(xmodel)

###########################################

def reado(ofilename, en = 1):

   keVnm = 1.23983

   ox, oy, xerr1, xerr2, yerr1, yerr2 = [],[],[],[],[],[]
   ux, uy, uerr1, uerr2, uyerr1, uyerr2 = [],[],[],[],[],[]
   fin = file(ofilename)
   lines = [line.strip().split() for line in fin.readlines() if not line.strip().startswith('#')]
   fin.close()
   for line in lines:

      wl = float(line[0])
      wlmin = float(line[0])+float(line[2])
      wlmax = float(line[0])+float(line[1])

      energy = keVnm*10/wl

      def energy(wl):
         return keVnm*10/wl

      if float(line[4]) != 0:
         if en == 1:
             ox.append(energy(wl))
             xerr2.append(energy(wlmin)-energy(wl))
             xerr1.append(energy(wl)-energy(wlmax))

             flux = 10**((23.9-float(line[3]))/2.5)
             oy.append(flux)
             ye2 = -10**((23.9-float(line[3])-float(line[4]))/2.5)+flux
             ye1 = +10**((23.9-float(line[3])+float(line[4]))/2.5)-flux
             yerr1.append(ye1)
             yerr2.append(ye2)
             
         elif en == 0:
             ox.append(wl)
             xerr2.append(wlmax-wl)
             xerr1.append(wl-wlmin)
             oy.append(float(line[3]))
             yerr1.append(float(line[4]))
             yerr2.append(float(line[4]))
      else:
         if en == 1:
              ux.append(energy(wl))
              uerr2.append(energy(wlmin)-energy(wl))
              uerr1.append(energy(wl)-energy(wlmax))
              uy.append(10**((23.9-float(line[3]))/2.5))
              uyerr1.append(0)
              uyerr2.append(0)
              
         elif en == 0:
             ux.append(wl)
             uerr2.append(wlmax-wl)
             uerr1.append(wl-wlmin)
             uy.append(float(line[3]))
             uyerr1.append(0)
             uyerr2.append(0)

   return (array(ox), array(oy), array(xerr1), array(xerr2), array(yerr1), 
           array(yerr2), array(ux), array(uy), array(uerr1), array(uerr2), 
           array(uyerr1), array(uyerr2))

##########################################################
    
def smooth(x,window_len=11,window='hanning'):
       if x.ndim != 1:
           raise ValueError, "smooth only accepts 1 dimension arrays."
   
       if x.size < window_len:
           raise ValueError, "Input vector needs to be bigger than window size."
   
   
       if window_len<3:
           return x
   
   
       if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
           raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
   
   
       s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
       #print(len(s))
       if window == 'flat': #moving average
           w=numpy.ones(window_len,'d')
       else:
           w=eval('numpy.'+window+'(window_len)')
   
       y=numpy.convolve(w/w.sum(),s,mode='same')
       return y[window_len:-window_len+1]
   
