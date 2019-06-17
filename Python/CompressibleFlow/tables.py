#!/usr/bin/env python

""" Tables.py - Make tables  of atmospheric properties

Adapted by
    Richard J. Kwan, Lightsaber Computing
from original programs by
    Ralph L. Carmichael, Public Domain Aeronautical Software

Revision History
Date	     Vers Person Statement of Changes
2004 Oct 04  1.0  RJK    Initial program
"""

import sys, math

version = "1.0 (2004 Oct 04)"
greeting = "Tables - A Python program to compute atmosphere tables"
author = "Ralph L. Carmichael, Public Domain Aeronautical Software"
modifier = "Thomas"
farewell = "Four files added to your directory."
finalmess = "Normal termination of tables."

#   P H Y S I C A L   C O N S T A N T S

FT2METERS = 0.3048		# mult. ft. to get meters (exact)
KELVIN2RANKINE = 1.8		# mult deg K to get deg R
PSF2NSM	= 47.880258		# mult lb/sq.ft to get sq.m
SCF2KCM = 515.379		# mult slugs/cu.ft to get kg/cu.m
TZERO	= 288.15		# sea-level temperature, kelvins
PZERO	= 101325.0		# @param PZERO  sea-level pressure, N/sq.m     
RHOZERO	= 1.225			# sea-level density, kg/cu.m
AZERO	= 340.294		# speed of sound at S.L.  m/sec
BETAVISC = 1.458E-6		# viscosity constant
SUTH	= 110.4			# Sutherland's constant, kelvins

def LongUSTable():
    """
    Atmosphere table US format long
    """
    Itxt = open('us1py.prt', 'w')
    Itxt.write(' alt   sigma     delta    theta ')
    Itxt.write(' temp   press     dens       a    visc  k.visc\n')
    Itxt.write(' Kft                            ')
    Itxt.write(' degR lb/sq.ft  s/cu.ft     fps s/ft-s sq.ft/s\n')

    for i in range(-1, 57):
        altKm=5*i*FT2METERS
        (sigma, delta, theta) = Atmosphere(altKm)
        Itxt.write("%4i %9.3E %9.3E %6.4f " % 
        (5*i, sigma, delta, theta))
        temp=(KELVIN2RANKINE*TZERO)*theta
        pressure=(PZERO/PSF2NSM)*delta
        density=(RHOZERO/SCF2KCM)*sigma
        asound=(AZERO/FT2METERS)*math.sqrt(theta)
        Itxt.write("%5.1f %8.3E %9.3E %6.1f" %
        (temp, pressure, density, asound))
        viscosity=(1.0/PSF2NSM)*MetricViscosity(theta)
        kinematicViscosity=viscosity/density
        Itxt.write("%6.3f %8.2E\n" %
        (1.0E6*viscosity, kinematicViscosity))

    Itxt.close()

def ShortUSTable():
    """
    Atmosphere table US short long
    """
    Itxt = open('us2py.prt', 'w')
    Itxt.write(' alt  sigma  delta  theta ')
    Itxt.write(' temp  press    dens     a     visc   k.visc ratio\n')
    Itxt.write(' Kft                      ')
    Itxt.write(' degR   psf   s/cu.ft   fps s/ft-sec sq.ft/s   1/ft\n')

    for i in range(-1, 66):
        altKm=i*FT2METERS
        (sigma, delta, theta) = SimpleAtmosphere(altKm)
        Itxt.write("%4i %6.4f %6.4f %6.4f" % 
        (i, sigma, delta, theta))
        temp=KELVIN2RANKINE*TZERO*theta
        pressure=PZERO*delta/47.88
        density=RHOZERO*sigma/515.379
        asound=(AZERO/FT2METERS)*math.sqrt(theta);
        Itxt.write("%6.1f %6.1f %9.7f %6.1f" %
        (temp, pressure, density, asound))
        viscosity=(1.0/PSF2NSM)*MetricViscosity(theta)
        kinematicViscosity=viscosity/density
        vratio=asound/kinematicViscosity
        Itxt.write("%6.3f %8.2E %4.2f\n"%
        (1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio))

    Itxt.close()

def LongSITable():
    """
    Atmosphere table SI format long
    """
    Itxt = open('tableSI.txt', 'w')
    Itxt.write(' alt    sigma      delta    theta ')
    Itxt.write(' temp   press      dens     a   visc   k.visc\n')
    Itxt.write('  Km                              ')
    Itxt.write('   K   N/sq.m    kg/cu.m  m/sec kg/m-s sq.m/s\n')

    for i in range(-1,44):
        altKm=2*i
        (sigma, delta, theta) = Atmosphere(altKm)
        Itxt.write("%4i %8.4E %8.4E %6.4f" %
        (2*i, sigma,delta,theta))
        temp=TZERO*theta
        pressure=PZERO*delta
        density=RHOZERO*sigma
        asound=AZERO*math.sqrt(theta)
        Itxt.write("%6.1f %8.3E %8.3E %5.1f" %
        (temp,pressure,density,asound))
        viscosity=MetricViscosity(theta)
        kinematicViscosity=viscosity/density
        Itxt.write("%6.2f %8.2E\n" %
        (1.0E6*viscosity, kinematicViscosity))

    Itxt.close()

def ShortSITable():
    """
    Atmosphere table SI format short
    """
    Itxt = open('si2py.prt', 'w')
    Itxt.write(' alt  sigma  delta  theta ')
    Itxt.write(' temp  press  dens   a    visc   k.visc ratio\n')
    Itxt.write('  Km                      ')
    Itxt.write('   K  N/sq.m  kcm  m/sec kg/m-s  sq.m/s  1/m\n')

    for i in range(-1, 41):
        altKm=0.5*i
        (sigma, delta, theta) = SimpleAtmosphere(altKm)
        Itxt.write("%4.1f %6.4f %6.4f %6.4f" %
        (altKm, sigma, delta, theta))
        temp=TZERO*theta
        pressure=PZERO*delta
        density=RHOZERO*sigma
        asound=AZERO*math.sqrt(theta)
        Itxt.write("%6.1f %6.0f %5.3f %5.1f" %
        (temp, pressure,density,asound))
        viscosity=MetricViscosity(theta)
        kinematicViscosity=viscosity/density
        vratio=asound/kinematicViscosity
        Itxt.write("%6.2f %8.2E %5.2f\n" %
        (1.0E6*viscosity, kinematicViscosity, 1.0E-6*vratio))

    Itxt.close()

def Atmosphere(alt):
    """ Compute temperature, density, and pressure in standard atmosphere.
    Correct to 86 km.  Only approximate thereafter.
    Input:
	alt	geometric altitude, km.
    Return: (sigma, delta, theta)
	sigma	density/sea-level standard density
	delta	pressure/sea-level standard pressure
	theta	temperature/sea-level std. temperature
    """

    REARTH = 6369.0		# radius of the Earth (km)
    GMR = 34.163195
    NTAB = 8			# length of tables

    htab = [ 0.0,  11.0, 20.0, 32.0, 47.0,
	51.0, 71.0, 84.852 ]
    ttab = [ 288.15, 216.65, 216.65, 228.65, 270.65,
	270.65, 214.65, 186.946 ]
    ptab = [ 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3, 1.0945601E-3,
	6.6063531E-4, 3.9046834E-5, 3.68501E-6 ]
    gtab = [ -6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0, 0.0 ]

    h = alt*REARTH/(alt+REARTH)	# geometric to geopotential altitude

    i=0; j=len(htab)
    while (j > i+1):
        k = int((i+j)/2)
        if h < htab[k]:
            j = k
        else:
            i = k
    tgrad = gtab[i]		# temp. gradient of local layer
    tbase = ttab[i]		# base  temp. of local layer
    deltah=h-htab[i]		# height above local base
    tlocal=tbase+tgrad*deltah	# local temperature
    theta = tlocal/ttab[0]	# temperature ratio

    if 0.0 == tgrad:
        delta=ptab[i]*math.exp(-GMR*deltah/tbase)
        sigma=delta/theta
    else:
        delta=ptab[i]*math.pow(tbase/tlocal, GMR/tgrad)
        sigma = delta/theta
    return ( sigma, delta, theta )

def SimpleAtmosphere(alt):
    """ Compute temperature, density, and pressure in simplified
    standard atmosphere.

    Correct to 20 km.  Only approximate thereafter.

    Input:
	alt	geometric altitude, km.
    Return: (sigma, delta, theta)
	sigma	density/sea-level standard density
	delta	pressure/sea-level standard pressure
	theta	temperature/sea-level std. temperature
    """
    REARTH = 6369.0		# radius of the Earth (km)
    GMR = 34.163195		# gas constant

    h = alt*REARTH/(alt+REARTH)	# geometric to geopotential altitude

    if h<11.0:		# troposphere
        theta=(288.15-6.5*h)/288.15
        delta=math.pow(theta, GMR/6.5)
        sigma=delta/theta
    else:		# stratosphere
        theta=216.65/288.15
        delta=0.2233611*math.exp(-GMR*(h-11.0)/216.65)
        sigma = delta/theta
    return ( sigma, delta, theta )

def MetricViscosity(theta):
    """
    Metric viscosity for Air
    """
    t=theta*TZERO
    return BETAVISC*math.sqrt(t*t*t)/(t+SUTH)

def main():
    """
    Main to calculate atmosphere tables
    """
    print("Executing ", sys.argv[0])
    print(greeting)
    print(author)
    if modifier != "":
        print("Modified by ", modifier)
        print("    version ", version)
        LongUSTable()
        ShortUSTable()
        LongSITable()
        ShortSITable()
    print(farewell)
    print(finalmess)


