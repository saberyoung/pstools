"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
import numpy as np

def moon_phase(month, day, year):
    ages = [18, 0, 11, 22, 3, 14, 25, 6, 17, 28, 9, 20, 1, 12, 23, 4, 15, 26, 7]
    offsets = [-1, 1, 0, 1, 2, 3, 4, 5, 7, 7, 9, 9]
    description = ["new (totally dark)",
      "waxing crescent (increasing to full)",
      "in its first quarter (increasing to full)",
      "waxing gibbous (increasing to full)",
      "full (full light)",
      "waning gibbous (decreasing from full)",
      "in its last quarter (decreasing from full)",
      "waning crescent (decreasing from full)"]
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    if day == 31:day = 1
    days_into_phase = ((ages[(year + 1) % 19] +
                        ((day + offsets[month-1]) % 30) +
                        (year < 1900)) % 30)
    index = int((days_into_phase + 2) * 16/59.0)    
    if index > 7:index = 7
    status = description[index]
    # light should be 100% 15 days into phase
    light = int(2 * days_into_phase * 100/29)
    if light > 100:
        light = abs(light - 200);
    date = "%d%s%d" % (day, months[month-1], year)
    return date, status, light

def slew_angle(ra,dec,ral,decl):
    ral, decl, idl = np.array(ral), np.array(decl), []
    while len(idl)<len(ral):
        dist = np.sqrt(((ral-ra)*15*np.cos(dec))**2+\
                       (decl-dec)**2)
        _new,_nn = False,0
        while not _new:
            if _nn < len(dist):
                _ids = np.argsort(dist)[_nn]
            else:
                _new=True
            if _ids in idl:
                _nn+=1
            else:
                _new=True
                idl.append(_ids)
            ra,dec = ral[_ids], decl[_ids]
    return np.array(idl)

def sunset(_timenow,_obs,_lsun):

    flag = True
    _deltat, _nn, _mm = .2, 5, 0
    while flag:
        _time = _timenow + astropy.time.TimeDelta(_deltat*_nn*_mm*3600., format='sec') 
        altazframe = astropy.coordinates.AltAz(obstime=_time, location=_obs)                                               
        sunaltaz = astropy.coordinates.get_sun(_time).transform_to(altazframe)
        if sunaltaz.alt.value <= _lsun:
            _mm =  int((_mm*_nn)/(_nn-1))
            _nn -= 1
        else:
            _mm += 1
        if _nn == 1: flag = False
        if _deltat*_nn*_mm >= 24: 
            print ('!!! Error: check sunset time...')
            flag = False
    return _time

def innight(_timel,_obs,_lsun):
    _timel1 = []
    for _time in _timel:
        altazframe = astropy.coordinates.AltAz(obstime=_time, location=_obs)                                               
        sunaltaz = astropy.coordinates.get_sun(_time).transform_to(altazframe)
        if sunaltaz.alt.value <= _lsun: _timel1.append(_time)
    return _timel1
