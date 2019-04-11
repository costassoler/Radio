import ugradio
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

def gal2hor(l, b, jd=None):
    """
    Converts galactic coordinates to horizontal coordinates at Leuschner. 
    
    Arguments:
    l, Galactic longitude (degrees).
    b, Galactic latitude (degrees).
    jd, Julian date (default=now).
    
    Returns: 
    alt, Horizontal altitude (degrees)
    az, Topocentric azimuth (degrees)
    """
    equ = equ2gal(l,b)
    leo = ugradio.leo
    hor = ugradio.coord.get_altaz(equ[0], equ[1], jd=jd, lat=leo.lat, lon=leo.lon, alt=leo.alt)
    alt = hor[0]
    if alt < 0:
        alt = alt + 90
    if alt > 90:
        alt = alt - 90
    azi = hor[1]
    if azi < 0:
        azi = alt + 360
    if azi > 360:
        azi = azi - 360
    return (alt, azi)

def danger(l, b):
    """
    Determines if given galtctic coordinates are safe or out of range for Leuschner. 
    
    Arguments:
    l, Galactic longitude (degrees).
    b, Galactic latitude (degrees).
    
    Returns:
    alt, Horizontal altitude (degrees).
    azi, Horizontal azimuth (degrees).
    """
    min_alt, max_alt = 15, 85
    min_azi, max_azi = 5, 350
    hor = gal2hor(l, b)
    alt, azi = hor[0], hor[1]
    print '(alt, azi) =',hor
    if alt < min_alt or alt > max_alt:
        print 'Coord Lord says: ALT OUT OF RANGE. DO NOT ATTEMPT POINTING.'
    if azi < min_azi or azi > max_azi:
        print 'Coord Lord says: AZI OUR OF RANGE. DO NOT ATTEMPT POINTING.'
    if alt > min_alt and alt < max_alt and azi > min_azi and azi < max_azi:
        print 'Coord Lord says: SAFE. GO FORTH AND POINT.'

