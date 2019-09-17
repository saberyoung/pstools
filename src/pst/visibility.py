"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input

def airmass(ra, dec, lat, lon, alt,
            time_input, airmass_min=1, airmass_max=5.8):
    """ Airmass calculation at a given time in a particular site for an input sky position (ra, dec) in degrees.
    The airmass is calculated in the range [airmass_min, airmass_max]."""
    
    import astropy
    from astropy import units as u
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.time import TimeDelta
    from astropy.time import Time

    observatory = astropy.coordinates.EarthLocation(lat = lat*u.deg,
                                                    lon = lon*u.deg,
                                                    height = alt*u.m)
    sky_coord = SkyCoord(ra = ra*u.deg,
                         dec=dec*u.deg, frame='icrs')     
    time = Time(time_input) 
    altaz = sky_coord.transform_to(AltAz(obstime=time,
                                         location=observatory))                                                                             
    airmass_value = altaz.secz                  
        
    if airmass_value < airmass_min or airmass_value > airmass_max:
        airmass_value =  "nan"
        return airmass_value
    else:     
        airmass_value = round(float(airmass_value), 2)
        return airmass_value

def prob_observable(m, header, _lat, _lon, _height, timelist):
    """
    Determine the integrated probability contained in a gravitational-wave
    sky map that is observable from a particular ground-based site at a
    particular time.

    Bonus: make a plot of probability versus UTC time!
    """

    # Determine resolution of sky map
    npix = len(m)
    nside = hp.npix2nside(npix)

    # Geodetic coordinates of observatory
    observatory = astropy.coordinates.EarthLocation(
        lat=_lat*u.deg, lon=_lon*u.deg, height=_height*u.m)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    theta, phi = hp.pix2ang(nside, np.arange(npix))

    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
        ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)

    prob = []
    for time in timelist:
        # Alt/az reference frame at observatory, now
        frame = astropy.coordinates.AltAz(obstime=time, location=observatory)    

        # Transform grid to alt/az coordinates at observatory, now
        altaz = radecs.transform_to(frame)

        # Where is the sun, now?
        sun_altaz = astropy.coordinates.get_sun(time).transform_to(altaz)

        # How likely is it that the (true, unknown) location of the source
        # is within the area that is visible, now? Demand that sun is at
        # least 18 degrees below the horizon and that the airmass
        # (secant of zenith angle approximation) is at most 2.5.
        prob.append(m[(sun_altaz.alt <= -18*u.deg) & (altaz.secz <= 2.5)].sum())

    # Done!
    return prob

def prob_obs(mapligo, header,lat1, lon1, height1,lat2, lon2, height2,_tgw,_tnow, _tlast):     
    _ttlist,_xticks = [],[]
    for ii in range(_tlast):        
        _tt = _tnow + astropy.time.TimeDelta(3600*ii, format='sec')
        _ttlist.append(_tt)
        _xticks.append('\n'.join(_tt.utc.iso.split()))   
    prob1 = prob_observable(mapligo, header,lat1, lon1, height1, _ttlist)
    prob2 = prob_observable(mapligo, header,lat2, lon2, height2, _ttlist)   
    ax = pl.axes([.1,.3,.8,.6])  
    ax.plot(range(_tlast),prob1,'r-',label='PROMPT5')   
    ax.plot(range(_tlast),prob2,'g-',label='MO')     
    ax.set_xticks(range(_tlast))
    ax.set_xticklabels(_xticks, rotation=60,fontsize=10)   
    ax.set_title('GW at %s'%_tgw)   
    plt.ylabel('Probability')
    plt.xlabel('time')
    ax.legend()
    return prob1[0],prob2[0]

def airmass_step(self, ra, dec):
    """Airmass calculation in step of one hour for an input sky position (ra, dec) in degrees.
    10 steps are performed: step <10; dt = 1h."""

    obs_time = Time(self.obs_time)
        
    while self.step < self.end_step:
        time_input = obs_time + self.step*self.dt
        val = self.airmass(ra, dec, self.altitude, self.longitude, self.altitude,
                           time_input)           
        self.airmass_list.append(val)
        self.time_list.append(str(time_input))
        self.step+=1
            
    return self.airmass_list, self.time_list
