"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input

def slew_angle(alt1,az1,alt2,az2):
    return mt.acos(mt.sin(alt1)*mt.sin(alt2)+mt.cos(alt1)*mt.cos(alt2)*mt.cos(abs(az1-az2)))

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

def prob_obs_hpmap(m, header, _lat, _lon, _height, timelist):
   
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

def prob_obs_galaxies(_id,ral,decl,score,_order,nob,limnob,_obslist,_tflist,timenow,timelast,_limsun,_limmoon,_solorobj,_limsolor,_limalt):

    start_time = time.time()
    from collections import defaultdict

    # how much longer for cgecking before night
    _deltat0 = 30*60 # every 30 mins

    # how many OBs can be taken during this time
    _nn = nob

    # read telescopes    
    _deltat = defaultdict(list)
    _duringnight,fdone1,fdone2,_nnight,_telon,_gobs = {},{},{},{},{},{}

    for ntel in _obslist:

        # average time per pointings
        tf = _tflist[ntel]

        # if during night
        _duringnight[ntel]=False

        # log galaxy list
        fdone1[ntel]={}

        # split fdone for different telescopes
        fdone2[ntel]=[]

        # number of tf has been observed
        _nnight[ntel]=0

        # if telescope can be on/available, [On/Off , time]
        _telon[ntel]=[True,0]

        # create delta t list for sampling every some time
        for jj in np.arange(0,timelast*3600,_tflist[ntel]*_nn): _deltat[jj].append(ntel)

    _deltatlist = np.sort(_deltat.keys())
    _deltatlist = np.diff(_deltatlist)

    ''' get OBs'''
    # Convert to RA, Dec
    radecs = astropy.coordinates.SkyCoord(ra=ral*u.deg, dec=decl*u.deg)
    
    # record done fields
    fdone,ii,_iddeltat = [],0,0

    while True:
       
        # stop till: 
        # - 1. end of night;
        # - 2. used up all available pointings, since for some telescope, number of pointings is limited

        iih = float('%.2f'%(ii/3600.))
        _tt = timenow + astropy.time.TimeDelta(ii, format='sec') 
        print('%s hours later (%s secs)'%(iih,ii))
       
        if iih > 3*24:
            # outside 3 days
            logging.info('Done!')
            print('#Done...')
            try:return fdone1,fdone2,_score
            except:return fdone1,fdone2,0

        ### 0 - check telescope available or not
        for ntel in _obslist:
            # if sometime later turn off telon
            if not _telon[ntel][0] and abs(ii - _telon[ntel][1]) >= _tflist[ntel]:_telon[ntel] = [True,ii]

        for ntel in _obslist:
            # check if available fields visible for any telescope

            # if telescope has enough pointings, if yes, go on         
            _tnobsum = 0
            if len(fdone1[ntel].keys())>0:
                for _tnob in fdone1[ntel]:
                    for _tnobnum in range(len(fdone1[ntel][_tnob])):
                        _tnobsum+=1                    
            if _tnobsum >= limnob[ntel]:continue           

            # read observatory
            observatory = _obslist[ntel]
            _show = (' - for telescope %s'%ntel)

            # Alt/az reference frame at observatory, now
            frame = astropy.coordinates.AltAz(obstime=_tt, location=observatory)           

            # Transform grid to alt/az coordinates at observatory, now
            altaz = radecs.transform_to(frame)

            # Where is the sun, now?
            sun_altaz = astropy.coordinates.get_sun(timenow).transform_to(altaz)

            # sun<-12, alt<=2.5        
            _cond = np.logical_and(np.logical_and(altaz.secz <= _limalt, altaz.secz >= 1),sun_altaz.alt <= _limsun*u.deg)
            gobs = _id[_cond]
            gradecs = radecs[_cond]

            if len(gobs) > 0:
                _show += ('\t%s fields visible for telescope %s'%(len(gobs),ntel))
                _duringnight[ntel]=True
                _nnight[ntel]+=1
            else: _show += ('\tNo fields visible for telescope %s'%ntel)            

            # check if there're available telescopes            
            if _telon[ntel][0]:_show += (', which is avalable now')
            else:_show += (', which is not avalable now')
            print(_show)              

            #### 1- go on or not
            if len(gobs) > 0 and _telon[ntel][0]:pass
            else:continue

            # record telon not available from ii to ii+tf
            _telon[ntel] = [False,ii]

            # Where is moon, now?
            # astropy.coordinates.solar_system_ephemeris.bodies
            if _limmoon == 'auto':
                if 'T' in str(_tt):_year, _month, _day = str(_tt.value).split('T')[0].split('-')
                else:_year, _month, _day = str(_tt.value).split()[0].split('-')
                _date, _status, _light = moon_phase(int(_month),int(_day),int(_year))               

                ''' limitation for the moon / TBD '''
                if _light>=80: _limmoon = 30
                elif _light<80 and _light>=40: _limmoon = 20
                elif _light<40 and _light>=10: _limmoon = 10
                else:_limmoon = 1

            _loc = astropy.coordinates.get_body('moon', _tt, observatory).transform_to(altaz)
            sep = gradecs.separation(_loc)              
            _length = len(gobs[np.where(sep.arcsecond<_limmoon)])
            if _length>0:
                print('\t-remove %s sources due to moon'%(_length))
                gobs = np.delete(gobs,np.where(sep.arcsecond<_limsolor))              
                gradecs = radecs[np.isin(_id, gobs)]

            # Where is the other sources inside solor system          
            for _source in _solorobj:

                if _source == 'moon':continue

                # astropy.coordinates.solar_system_ephemeris.bodies
                _loc = astropy.coordinates.get_body(_source, _tt, observatory).transform_to(altaz)
                sep = gradecs.separation(_loc)              
                _length = len(gobs[np.where(sep.arcsecond<_limsolor)])
                if _length>0:
                    print('\t-remove %s sources due to %s'%(_length,_source))
                    gobs = np.delete(gobs,np.where(sep.arcsecond<_limsolor))              
                    gradecs = radecs[np.isin(_id, gobs)]

            #### 2 - go on or not
            if len(gobs) == 0:               
                print('\t!!!for telescope %s No fields visible anymore after removing conditions!'%ntel)
                continue
            else:_gobs[ntel] = gobs

            # print infos
            _ss = '\t%i fields visible for %s, when %s'%(len(gobs),ntel,_tt)
            if iih==0:_ss+=' (Now!)'
            else:_ss+=' (%.2f h)'%(iih)
            print(_ss)                

            #### 3 - ranking fields
            if _order == 1:
                # rank with prob
                pass

            elif _order == 2:
                # from west to east
                _idrank = np.argsort(gradecs.ra)
                gobs = gobs[_idrank]
                gradecs = gradecs[_idrank]

            elif _order in [3,4]:                              

                # from west to east with slewing angle                
                ra1,dec1,_num=np.array(gradecs.ra),\
                               np.array(gradecs.dec),1                                               
                if _order==3: _ids = np.argmin(ra1)
                else:_ids = 0                
                ras,decs = ra1[_ids],dec1[_ids]               
                _idrank = [_ids]
                while _num<len(ra1):   
                    _num+=1
                    dist = np.sqrt(((ra1-ras)*15*mt.cos(decs))**2+(dec1-decs)**2)
                    _OK,_nn=False,1
                    while not _OK:
                        _ids = np.argsort(dist)[_nn]
                        if _ids in _idrank:_nn+=1
                        else:_OK=True        
                    ras,decs = ra1[_ids], dec1[_ids]      
                    _idrank.append(_ids)                    
                _idrank = np.array(_idrank)
                gobs = gobs[_idrank]
                gradecs = gradecs[_idrank]

            else:
                logging.info('Error: wrong order %s!!!'%_order)
                sys.exit()

            # alt az
            galtaz = gradecs.transform_to(frame)

            #
            fdone1[ntel][iih]=[]

            # available fields           
            fdone0 = []
            for _ii in gobs:
                if _ii in fdone:continue
                else:fdone0.append(_ii)
          
            '''
            # adjust n per time
            if len(fdone0)>0 and len(fdone0)<_nn:_nn = len(fdone0)
            '''

            # if too less fields, wait a little bit
#            if _order == 1 and len(fdone0) < len(score)*0.5:
#                fdone2[ntel].append(False)
#                fdone1[ntel][iih].append(False)

            if _nn <= len(fdone0):
                # can only taken _nn frames, the rest go for the next step 
                for _ii in fdone0[:_nn]:
                    fdone.append(_ii)
                    fdone2[ntel].append(_ii)
                    fdone1[ntel][iih].append([_ii,gradecs[np.where(gobs==_ii)],galtaz[np.where(gobs==_ii)]])
            else:
                fdone2[ntel].append(False)
                fdone1[ntel][iih].append(False)

                '''
                # all taken then can go for repetation
                if len(fdone0)==0:
                    # till the point, all galaxies have beed observed            
                    print('!!!Warning: all visiable fields now have already taken!!!\nrepeating previous fields!!!')
                    fdone0=gobs
                else:
                    print('!!!Warning: not enough visiable fields could be taken!!!\nfinishing and repeating previous fields!!!')
                    fdone0 = fdone0 + list(gobs)
                _ok,_mm = True,0
                while _ok:                  
                    for _ii in fdone0:                       
                        if _mm > _nn:
                            _ok=False                        
                            continue
                        _mm+=1
                        fdone.append(_ii)  
                        fdone2[ntel].append(_ii)  
                        fdone1[ntel][iih].append([_ii,gradecs[np.where(gobs==_ii)],galtaz[np.where(gobs==_ii)]])
                '''

            if len(fdone)>0:
                _score = sum(score[np.where(_id==jj)] for jj in np.array(fdone)[np.where(fdone)])[0]
                print('\t\tcoverage:%.5e'%_score)           

                print(fdone)

        ### 4 - add time for the next step
        # which telescope has visible fields
        _telf = []
        for _ee in _gobs:
            if len(_gobs[_ee])>0:_telf.append(_ee)

        if _iddeltat >= len(_deltatlist):
            # till the end of time u set
            logging.info('Done!')
            print('#Done...')
            try:return fdone1,fdone2,_score
            except:return fdone1,fdone2,0

        '''
        if len(_telf) == 0 and sum(_duringnight.values()) == len(_obslist.keys()):
            # no fields visible anymore, and all telescopes have already during night
            # End!
            logging.info('End of night!')
            print('#End of night...')
            try:return fdone1,fdone2,_score
            except:return fdone1,fdone2,0
        '''

        # too much longer, sth wrong!!!
        if  time.time() - start_time > 20*60:
            logging.info('too long waiting!')
            try:return fdone1,fdone2,_score
            except:return fdone1,fdone2,0

        if len(_telf) == 0:                   
            # no fields visible for no telescopes
            # go on...
            print('#go on %s secs'%_deltat0)
            ii += _deltat0 

        else:
            # visible fields for at least one telescope                        
            _deltat1 = _deltatlist[_iddeltat]
            print('#go on %s secs'%_deltat1)
            ii += _deltat1
            _iddeltat+=1
