import numpy as np
from __future__ import division

def HMS(hour, minute, second):
    """
    Converts RA from hour-minute-second notation to decimal notation.
    
    Arguments:
    hour: RA value of degree. 
    minute: RA value of minute.
    second: RA value of second.
    
    Returns:
    Decimal representation of the degree-minute-second declination.
    """
    minute = minute/60
    second = second/3600
    decimal = 15*(hour + minute + second)
    return decimal
  
    
def DMS(degree, arcmin, arcsec):
    """
    Converts dec from degree-minute-second notation to decimal notation.  
    
    Arguments:
    degree: dec value of degree. 
    arcmin: dec value of arcminute.
    arcsec: dec value of arcsecond.
    
    Returns:
    Decimal representation of the degree-minute-second declination.
    """
    arcmin = arcmin/60
    arcsec = arcsec/3600
    if degree < 0:
        decimal = -(np.abs(degree) + arcmin + arcsec)
    else:
        decimal = np.abs(degree) + arcmin + arcsec
    return decimal

def gal2equ(l, b):
    """
    Converts galactic coordinates to equatorial coordinates, according to the epoch J2000. 
    
    Arguments:
    l, Galactic longitude (degrees).
    b, Galactic latitude (degrees).
    
    Returns: 
    RA, Equatorial right ascension (degrees).
    dec, Equatorial declination (degrees).
    """
    l = np.radians(l)
    b = np.radians(b) 
    x = np.array([0.,0,0])
    x[0] = np.cos(b) * np.cos(l)
    x[1] = np.cos(b) * np.sin(l)
    x[2] = np.sin(b)
    R = np.array([-0.054876, -0.873437, -0.483835, 0.494109, -0.444830, 0.746982, -0.867666, -0.198076, 0.455984])
    R = R.reshape(3, 3)
    iR = np.linalg.inv(R)
    xp = np.dot(iR, x)
    RA = np.arctan2(xp[1], xp[0])
    dec = np.arcsin(xp[2])
    RA = np.degrees(RA)
    dec = np.degrees(dec)
    return (RA, dec)

def equ2gal(RA, dec):
    """
    Converts equatorial coordinates to galactic coordinates, according to the epoch J2000. 
    
    Arguments:
    RA, Equatorial right ascension (degrees).
    dec, Equatorial declination (degrees).
    
    Returns: 
    l, Galactic longitude (degrees).
    b, Galactic latitude (degrees).
    """
    RA = np.radians(RA)
    dec = np.radians(dec) 
    x = np.array([0.,0,0])
    x[0] = np.cos(dec) * np.cos(RA)
    x[1] = np.cos(dec) * np.sin(RA)
    x[2] = np.sin(dec)
    R = np.array([-0.054876, -0.873437, -0.483835, 0.494109, -0.444830, 0.746982, -0.867666, -0.198076, 0.455984])
    R = R.reshape(3, 3)
    xp = np.dot(R, x)
    l = np.arctan2(xp[1], xp[0])
    b = np.arcsin(xp[2])
    l = np.degrees(l)
    b = np.degrees(b)
    if l < 0:
        l = l + 360
    return (l, b)

