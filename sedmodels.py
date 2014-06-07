#!/usr/bin/env python

import sys
import math
import scipy

from pylab import exp, array
from scipy import interpolate
from matplotlib.patches import Polygon

def fill_between(ax, x, y1, y2, **kwargs):
    # add x,y2 in reverse order for proper polygon filling
    verts = zip(x,y1) + [(x[i], y2[i]) for i in range(len(x)-1,-1,-1)]
    poly = Polygon(verts, **kwargs)
    ax.add_patch(poly)
    ax.autoscale_view()
    return poly


def dla(z, wl1 ,wl2, NH):
    
    c = 299792458.
    NH = NH*1E4
    lambda_a = 121.6*1E-9
    nu_a = c/lambda_a
    
    f_a = 0.4162
    gu = 3.
    gl = 1.
    delta_cla = 1.503E9
    delta_a = 3*(gl/gu)*f_a*delta_cla
    pi = 3.1416

    const = 3*lambda_a**2*f_a*delta_cla/8/pi
    tau1 = []
    sigma1 = []
    tau2 = []
    sigma2 = []

    for wave in wl1:
        nu = c/wave/1E-9*(1+z)
        sigma_a1 = const *delta_a*((nu/nu_a)**4/(4*pi**2*(nu-nu_a)**2+delta_a**2*(nu/nu_a)**6/4.))
        sigma1.append(sigma_a1)
        tau_a1 = NH*sigma_a1
        tau1.append(tau_a1)

    for wave in wl2:
        nu = c/wave/1E-9*(1+z)
        sigma_a2 = const *delta_a*((nu/nu_a)**4/(4*pi**2*(nu-nu_a)**2+delta_a**2*(nu/nu_a)**6/4.))
        sigma2.append(sigma_a2)
        tau_a2 = NH*sigma_a2
        tau2.append(tau_a2)
    return tau1, tau2


def lyalpha(z):
   
    nw = 1000.
    w1 = 1030.*(1+z)
    w2 = 1210.*(1+z)
    wstep = (w2-w1)/nw
    
    wtemp=[]
    ptau=[]

    for i in range(1, int(nw)):
       
        wavetemp = w1+(i-1)*wstep
        wtemp.append(wavetemp)
        ptau.append(math.exp(-1.8E-3*(wavetemp/1216.)**3.92))      

    da = (1./(180.*(1.+z)))*trapz1(wtemp, ptau, nw)

    w1 = 980.*(1+z)
    w2 = 1025.*(1+z)

    wstep = (w2-w1)/nw
    
    wtemp=[]
    ptau=[]

    for i in range(1, int(nw)):
        wavetemp = w1+(i-1)*wstep
        wtemp.append(wavetemp)
        ptau.append(exp(-0.9E-3*(wavetemp/1026.)**3.92))

    
    db = (1./(45.*(1.+z)))*trapz1(wtemp, ptau, nw)


    w1 = 951.*(1+z)
    w2 = 972.*(1+z)
    wstep = (w2-w1)/nw
    
    wtemp=[]
    ptau=[]

    for i in range(1, int(nw)):
        wavetemp = w1+(i-1)*wstep
        wtemp.append(wavetemp)
        ptau.append(exp(-0.6E-3*(wavetemp/973.)**3.92))

    dc = (1./(21.*(1.+z)))*trapz1(wtemp, ptau, nw)

    w1 = 913.*(1+z)
    w2 = 949.*(1+z)
    wstep = (w2-w1)/nw
    
    wtemp=[]
    ptau=[]

    for i in range(1, int(nw)):
        wavetemp = w1+(i-1)*wstep
        wtemp.append(wavetemp)
        ptau.append(exp(-4.5E-5*(wavetemp/950.)**3.92))

    dd = (1./(36.*(1.+z)))*trapz1(wtemp, ptau, nw)
   
    if da > 1:
            da = 1
    if da < 0:
            da = 0
    if db > 1:
            db = 1
    if db < 0:
            db = 0
    if dc > 1:
            dc = 1
    if dc < 0:
            dc = 0
    if dd > 1:
            dd = 1
    if dd < 0:
            dd = 0
         
    return da, db, dc, dd


def trapz1(x,y,n):
    trapz1=0.

    if n == 1:
        trapz1=y[0]*x[0]/2.
    else:
        for j in range(1, int(n)-1):
            trapz1=trapz1+abs(((x[j]-x[j-1])*(y[j]+y[j-1]))/2.)
    return trapz1

def PLX(x, xnorm, pl):
    return xnorm*x**(-1*pl)

def BRPLX(x, xnorm, pl1, pl2, eb):
    brpl = []
    for en in x:
        if en < eb:
            brpl.append(xnorm*en**(-1*pl1))
        else:
            brpl.append(xnorm*eb**(pl2-pl1)*en**(-1*pl2))
    return array(brpl) 

def red_sed(en, flux, red, av, z):
    y = 1.23983/en/(1+z)
    if red == 'mw':   
        mw=Avlaws(y, 'mw')
        return flux*exp(-1/1.086*av*mw)
    if red == 'smc':
        smc=Avlaws(y, 'smc')
        return flux*exp(-1/1.086*av*smc)
    if red == 'lmc':  
        lmc=Avlaws(y, 'lmc')
        return flux*exp(-1/1.086*av*lmc)

def fmext(y, av = 1, rv = 3.08, gamma = 0.99, x0 = 4.596, c3 = 3.23, c4 = 0.41):
    ebv = av/rv
    c2 = -0.824 + 4.717/rv
    c1 = 2.030 - 3.007*c2
    xcutuv = 1000.0/270.0
    x = 1000./y
    ancpointsx = [0, 0.377, 0.820, 1.667, 1.828, 2.141, 2.4333]
    ancpointsy = [0, 0.265/3.1, 0.829/3.1, #, 2.688/3.1, 3.055/3.1, 3.806/3.1, 4.315/3.1]
                  -4.22809e-01/rv + 1.00270 + 2.13572e-04*rv**1,
                  -5.13540e-02/rv + 1.00216 - 7.35778e-05*rv**1,  
                  +7.00127e-01/rv + 1.00184 - 3.32598e-05*rv**1,
                  +1.19456/rv + 1.01707 - 5.46959e-03*rv**1 - 4.45636e-05*rv**2]

    # Fitzpatrick and Massa 90, full parametrization
    # A_V = E_B-V * R_V
    # yuv = k(lambda)
    # k(lambda) = E(lambda)/E(B-V), A(lambda)=E_B-V*R_V*k(lambda)
    # Convert to microns
    if x >= xcutuv:
        if x >= 5.9:
            fuv = 0.5392*(x-5.9)**2 + 0.05644*(x-5.9)**3
        else:
            fuv = 0
        yuv = c1 + c2*x + c3*(x**2/((x**2-x0**2)**2+x**2*gamma**2))+c4*fuv
    else:
        for uvwl in range(270,240,-10):
            x2 = 1000./uvwl
            yuv = c1 + c2*x2 + c3*(x2**2/((x2**2-x0**2)**2+x2**2*gamma**2))
            #print yuv
            ancpointsx.append(x2)
            ancpointsy.append(yuv) 
        fmspline = scipy.interpolate.splrep(ancpointsx, ancpointsy, s=0)
        yuv = scipy.interpolate.splev(array(x), fmspline)
        
    return yuv*av, ancpointsx, ancpointsy

def Avlaws(y, law):
    # y is wavelength in nm
    if law == 'mw':
        rv = 3.08
        a1, l1, b1, n1 =165.,47.,90.,2.
        a2, l2, b2, n2 =14.,80.,4.,6.5
        a3, l3, b3, n3 =0.045, 220.,-1.95,2.
        a4, l4, b4, n4 =0.002,9700.,-1.95,2.
        a5, l5, b5, n5 =0.002,18000.,-1.8,2.
        a6, l6, b6, n6 =0.012,25000.,0.,2.

        mw1=a1/((y/l1)**n1+(l1/y)**n1+b1)
        mw2=a2/((y/l2)**n2+(l2/y)**n2+b2)
        mw3=a3/((y/l3)**n3+(l3/y)**n3+b3)
        mw4=a4/((y/l4)**n4+(l4/y)**n4+b4)
        mw5=a5/((y/l5)**n5+(l5/y)**n5+b5)
        mw6=a6/((y/l6)**n6+(l6/y)**n6+b6)

        return rv*(mw1+mw2+mw3+mw4+mw5+mw6)*1.35649/3.08

    if law == 'smc':
        rv = 2.93
        a12, l12, b12, n12 =185.,42.,90.,2.
        a22, l22, b22, n22 =27.,80.,5.5,4.0
        a32, l32, b32, n32 =0.005, 220.,-1.95,2.
        a42, l42, b42, n42 =0.010,9700.,-1.95,2.
        a52, l52, b52, n52 =0.012,18000.,-1.8,2.
        a62, l62, b62, n62 =0.030,25000.,0.,2.

        smc1=a12/((y/l12)**n12+(l12/y)**n12+b12)
        smc2=a22/((y/l22)**n22+(l22/y)**n22+b22)
        smc3=a32/((y/l32)**n32+(l32/y)**n32+b32)
        smc4=a42/((y/l42)**n42+(l42/y)**n42+b42)
        smc5=a52/((y/l52)**n52+(l52/y)**n52+b52)
        smc6=a62/((y/l62)**n62+(l62/y)**n62+b62)

        return rv*(smc1+smc2+smc3+smc4+smc5+smc6)*1.3875/2.93

    if law == 'lmc':
        rv = 3.16
        a11, l11, b11, n11 =175.,46.,90.,2.
        a21, l21, b21, n21 =19.,80.,5.5,4.5
        a31, l31, b31, n31 =0.028, 220.,-1.95,2.
        a41, l41, b41, n41 =0.005,9700.,-1.95,2.
        a51, l51, b51, n51 =0.006,18000.,-1.8,2.
        a61, l61, b61, n61 =0.020,25000.,0.,2.

        lmc1=a11/((y/l11)**n11+(l11/y)**n11+b11)
        lmc2=a21/((y/l21)**n21+(l21/y)**n21+b21)
        lmc3=a31/((y/l31)**n31+(l31/y)**n31+b31)
        lmc4=a41/((y/l41)**n41+(l41/y)**n41+b41)
        lmc5=a51/((y/l51)**n51+(l51/y)**n51+b51)
        lmc6=a61/((y/l61)**n61+(l61/y)**n61+b61)

        return rv*(lmc1+lmc2+lmc3+lmc4+lmc5+lmc6)*1.3165/3.16

    if law == 'sn':
        a11, l11, b11, n11 =157.8,67.7,136.97,2.56
        a21, l21, b21, n21 =157.4,79.73,23.59,5.87
        a31, l31, b31, n31 =0.2778, 272.8,-1.87,0.87
        a41, l41, b41, n41 =0.005,9700.,-1.95,2.
        a51, l51, b51, n51 =0.006,18000.,-1.8,2.
        a61, l61, b61, n61 =0.020,25000.,0.,2.

        sn1=a11/((y/l11)**n11+(l11/y)**n11+b11)
        sn2=a21/((y/l21)**n21+(l21/y)**n21+b21)
        sn3=a31/((y/l31)**n31+(l31/y)**n31+b31)
        sn4=a41/((y/l41)**n41+(l41/y)**n41+b41)
        sn5=a51/((y/l51)**n51+(l51/y)**n51+b51)
        sn6=a61/((y/l61)**n61+(l61/y)**n61+b61)

        return (sn1+sn2+sn3+sn4+sn5+sn6)*1.0000

    if law == 'cal':
        #ff = k(lambda); A(lambda) = A_V*k(lambda)/R ; R = 4.05
        rv=4.05
        p11=1/0.11
        ff11=2.659*(-2.156+1.509*p11-0.198*p11**2+0.011*p11**3)+rv
        p12=1/0.12
        ff12=2.659*(-2.156+1.509*p12-0.198*p12**2+0.011*p12**3)+rv
        slope1=(ff12-ff11)/100.
        ff99=2.659*(-1.857+1.040/2.19)+rv
        ff100=2.659*(-1.857+1.040/2.2)+rv
        slope2=(ff100-ff99)/100.
        p = 1./(y*1E-3)
        if y > 2200:
            ff=(ff99+(y-2190.)*slope2)/rv
        elif (y <= 2200) and (y>630):
            ff=(2.659*(-1.857+1.040*p)+rv)/rv
        elif (y <= 630) and (y>90):
            ff=(2.659*(-2.156+1.509*p-0.198*p**2+0.011*p**3)+rv)/rv
        else:
            ff = 0
        return ff
        

def PL(en, red, pl1, av, z):
    y = 1.23983/en/(1+z)
    wlnorm = 1000/(1+z)
    pl = (y/wlnorm)**(pl1)
    if red == 0:      
        return pl
    else:
        redlaw=Avlaws(y, red)
        return pl*exp(-1/1.086*av*redlaw)

def BRPL(en, red, pl1, eb, s, pl2, av, z):
    y = 1.23983/en/(1+z)
    wlb = 1.23983/eb/(1+z)
    brpl = ((y/wlb)**(-pl1*s)+(y/wlb)**(-pl2*s))**(-1/s)
    if red == 0:      
        return brpl
    else:
        redlaw=Avlaws(y, red)
        return brpl*exp(-1/1.086*av*redlaw)

def TRBRPL(en, red, pl1, eb1, s1, pl2, eb2, s2, pl3, av, z):
    y = 1.23983/en/(1+z)
    wlb1 = 1.23983/eb1/(1+z)
    wlb2 = 1.23983/eb2/(1+z)
    brpl1 = ((y/wlb1)**(-pl1*s1)+(y/wlb1)**(-pl2*s1))**(-1/s1)
    brpl2 = ((wlb2/wlb1)**(-pl1*s1)+(wlb2/wlb1)**(-pl2*s1))**(-1/s1)*(y/wlb2)**pl3
    trbrpl = (brpl1**(-s2)+brpl2**(-s2))**(-1/s2)
    if red == 0:      
        return trbrpl
    else:
        redlaw=Avlaws(y, red)
        return trbrpl*exp(-1/1.086*av*redlaw)
        
def PLwl(x, x2, red, pl, av, z, c1=0, c2=0, c3=0, c4=0):

    y = x/(1+z)
    y2 = x2/(1+z)
    wlnorm = 200*(1+z)
  
    if red == 0:      
        return ((y/wlnorm)**(pl)), ((y2/wlnorm)**(pl))

    if red == 'mw':   
        mw=Avlaws(y, 'mw')        
        return ((y/wlnorm)**(pl))*exp(-1/1.086*av*mw), ((y2/wlnorm)**(pl))

    if red == 'smc':
        smc=Avlaws(y, 'smc')
        return ((y/wlnorm)**(pl))*exp(-1/1.086*av*smc), ((y2/wlnorm)**(pl))

    if red == 'lmc':   
        lmc = Avlaws(y, 'lmc')        
        return ((y/wlnorm)**(pl))*exp(-1/1.086*av*lmc), ((y2/wlnorm)**(pl))
        
    if red == 'drude':   

        drude1 = c1/((y/80)**c2+(80/y)**c2+c3)
        drude2 = (233*(1-c1/(6.88**c2+0.145**c2+c3)-c4/4.60))/((y/46)**2+(46/y)**2+90)
        drude3 = c4/((y/217.5)**2+(217.5/y)**2-1.95)
        drude=(drude1+drude2+drude3)
        return ((y/wlnorm)**(pl))*exp(-1/1.086*av*drude), ((y2/wlnorm)**(pl))


def BRPLwl(x, x2, red, pl1, eb, s, pl2, av, z, c1=0, c2=0, c3=0, c4=0):

    y = x/(1+z)
    y2 = x2/(1+z)
    brpl = ((y/eb)**(-pl1*s)+(y/eb)**(-pl2*s))**(-1/s)
    if red == 0:      
        return brpl
    else:
        redlaw=Avlaws(y, red)
        return brpl*exp(-1/1.086*av*redlaw), brpl


def PLerr(en, plmin, plmax, normmin, normmax):

    step = 20
    optplmin = []
    optplmax = []
    for energy in en:
        fluxmin = 1E10
        fluxmax = 0.
        norminc = abs(normmax-normmin)/float(step)
        pl = plmin
        plinc = abs(plmax-plmin)/float(step)
        for i in range(step):
            norm = normmin
            for j in range(step):
                flux = norm*(energy)**(-1*pl)
                #print flux
                if flux < fluxmin:
                    fluxmin = flux
                if flux > fluxmax:
                    fluxmax = flux
                norm+= norminc
            pl+= plinc
        optplmin.append(fluxmin)
        optplmax.append(fluxmax)

    return array(optplmin), array(optplmax)


def BRPLerr(en, pl1min, pl1max, normmin, normmax, pl2min, pl2max, ebmin, ebmax):
    step = 10
    optplmin = []
    optplmax = []
    for energy in en:
        fluxmin = 1E20
        fluxmax = 0.
        norminc = abs(normmax-normmin)/float(step)
        pl1 = pl1min
        pl1inc = abs(pl1max-pl1min)/float(step)
        pl2 = pl2min
        pl2inc = abs(pl2max-pl2min)/float(step)
        ebinc = abs(ebmax-ebmin)/float(step)
        norm = normmin
        for i in range(step):
            pl1 = pl1min
            for j in range(step):
                pl2=pl2min
                for k in range(step):
                    eb = ebmin
                    for l in range(step):
                        if energy < eb:
                            flux = norm*(energy)**(-1*pl1)
                        else:
                            flux = norm*eb**(pl2-pl1)*(energy)**(-1*pl2)
                        if flux < fluxmin:
                            fluxmin = flux
                        if flux > fluxmax:
                            fluxmax = flux
                        eb+=ebinc
                    pl2+=pl2inc
                pl1+=pl1inc
            norm+=norminc
        optplmin.append(fluxmin)
        optplmax.append(fluxmax)
  
    return array(optplmin), array(optplmax)
