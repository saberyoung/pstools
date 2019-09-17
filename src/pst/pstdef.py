"""############################################################################ 
2019/1/30 Start
define some functions in assiting/supporting
""" ############################################################################
from __future__ import print_function
from builtins import input

import time
import numpy as np
import healpy as hp
import math as mt
import pylab as pl
import matplotlib.pyplot as plt
import random
import os,sys,glob,shutil
import logging
from wxpy import Bot
import sqlconn

from astropy.io import fits
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
from astroquery.vizier import Vizier
import scipy.stats
from slackclient import SlackClient

from scp import SCPClient
from pst import scheme,priorization,pstplot,scheduler,configure,link
import pst
_pstpath = pst.__path__[0]

#################################################
def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))
    
def RadecToThetaphi(ra,dec):
    return mt.pi/2.-np.radians(dec),np.radians(ra)

def ThataphiToRadec(theta,phi):
    return np.degrees(phi),np.degrees(mt.pi/2.-theta)

def query_ebv(ra,dec,size=2,thresh=.25,verbose=False):
    """
    Take ra,dec and return the E(B-V) number
    URL query, see https://irsa.ipac.caltech.edu/applications/DUST
    """
    import urllib2
    import xmltodict
    url = "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?"
    url += "locstr=%.2f+%.2f+equ+j2000"%(ra,dec)
    if size<2:size=2
    if size>37.5:size=37.5
    url += '&regSize=%.2f'%size # size between 2.0 and 37.5    

    _file = urllib2.urlopen(url)
    data = _file.read()
    _file.close()

    _dict = xmltodict.parse(data)
    _ebv = _dict['results']['result'][0]['statistics']['meanValueSFD']

    ebvalue = float(_ebv.split('(mag)')[0])
    if ebvalue<thresh:
        if verbose:
            print("ra=%.2f dec=%.2f:\tebv=%.2f\tOK"%\
                  (ra,dec,float(_ebv.split('(mag)')[0])))
        return ebvalue
    else:
        if verbose:
            print("ra=%.2f dec=%.2f:\tebv=%.2f\tNo"%\
                  (ra,dec,float(_ebv.split('(mag)')[0])))
        return

def query_vizier(_catname,_size,_interactive,ramin,ramax,decmin,decmax,absmag,distmin,distmax,verbose=False):
    """
    Download galaxy catalog from vizier
    return calalog and its name
    """    

    # find catfile from catalog name
    catalog_list = Vizier.find_catalogs(_catname)    
    if len(catalog_list)==0:sys.exit('Error: unavalable catalogs')
        
    # specify
    if _catname=='GLADE':
        _raname,_decname,_magname,_distname,_columns,_cat = 'RAJ2000', 'DEJ2000', 'BMAG', 'Dist',\
                                                            ['PGC', 'GWGC', 'HyperLEDA', '2MASS', 'RAJ2000', 'DEJ2000', 'BMAG', 'Dist'],\
                                                            'VII/281'
    elif _catname=='GWGC':
        _name,_raname,_decname,_magname,_distname,_columns,_cat = 'Name','RAJ2000', 'DEJ2000','BMAG','Dist',\
                                                            ['Name','RAJ2000', 'DEJ2000','BMAG','Dist'],\
                                                            'VII/267'
    elif _catname=='NED':
        _name,_raname,_decname,_magname,_distname,_columns,_cat = 'MGC','RAJ2000', 'DEJ2000','Bmag','z',\
                                                            ['MGC','RAJ2000', 'DEJ2000','Bmag','z'],\
                                                            'VII/240'
    else:sys.exit('Error: outside galaxy catalogs')
 
    if _magname == 'Bmag':
        # for apparent mag, to abs mag
        #if readkeys.keys(catname)[2] == 'Bmag':mag00 = mag00-5*np.log10(dist00)-25   
        sys.exit('Need TBD for cat with apparent magnitude, now choose another cat with absolute magnitude!!!')

    # download catalog with vizier
    '''
    quite strange vizier:
    if one two conditions for dec selection,
    only the first one works!!!
    so, I set only one condition on dec here and use numpy to set for the second
    '''

    v = Vizier(columns=_columns,
               column_filters={_magname:'<'+str(absmag),\
                               _distname:'%s..%s'%(str(distmin),str(distmax)),\
                               _raname:'%s..%s'%(str(ramin),str(ramax)),\
                               _decname:'%s..%s'%(str(decmin),str(decmax))})
    v.ROW_LIMIT = _size
    catalogs = v.get_catalogs(_cat)[0]    

    if verbose: 
        print("%i galaxies selected from %s in %i sec"%\
              (len(catalogs),_catname,int(time.time()-start_time)))

    # return infos    
    if _catname=='GLADE':
        _name = []
        for ii in range(len(catalogs)):_name.append('%s:%s:%s:%s'%(catalogs['PGC'][ii], \
                                                                   catalogs['GWGC'][ii], \
                                                                   catalogs['HyperLEDA'][ii], \
                                                                   catalogs['_2MASS'][ii]))        
    else:_name = catalogs[_name]
    return _name,np.array(catalogs[_raname]),np.array(catalogs[_decname]),np.array(catalogs[_magname]),np.array(catalogs[_distname])

def rotate_map(hmap, rot_theta, rot_phi):
    """
    Take hmap (a healpix map array) and return another healpix map array 
    which is ordered such that it has been rotated in (theta, phi) by the 
    amounts given.
    """
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))) #theta, phi

    # Define a rotator
    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot)
    
    return rot_map

def vertices(ra,dec,fovw,fovh):
    """Finding the vertices of a FoV given the central location (ra[deg], dec[deg])
    and the FoV size (fovw [deg], fovh [deg]).
    """
    fovw,fovh = fovw/2.,fovh/2.
    vert_ra,vert_dec=[],[]
    ra_rad,dec_rad,fovw_rad,fovh_rad = np.deg2rad(ra), np.deg2rad(dec),\
                                       np.deg2rad(fovw), np.deg2rad(fovh)
    for i,j in zip([-fovw_rad, fovw_rad, fovw_rad, -fovw_rad],\
                     [fovh_rad, fovh_rad, -fovh_rad, -fovh_rad]):
        arg = -i/(np.cos(dec_rad)-j*np.sin(dec_rad))
        v_ra = np.rad2deg(ra_rad+np.arctan(arg))       
        v_dec = np.rad2deg(np.arcsin((np.sin(dec_rad)+j*np.cos(dec_rad))/(1+i**2+j**2)**0.5))
        vert_ra.append(v_ra)
        vert_dec.append(v_dec)
    return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], \
        vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]

def ipix_in_box(ra,dec,width,height,nside):

    width,height = width/2.,height/2.
    v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = vertices(ra, dec, width, height)
    ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                 [v1_dec, v2_dec, v4_dec, v3_dec])

    theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
    phi = np.deg2rad(ra_vertices)
    xyz = hp.ang2vec(theta, phi)
    ipix_fov_box = hp.query_polygon(nside, xyz)
    return ipix_fov_box

def contour(hpx,_contour1=.5,_contour2=.68,_contour3=.9,_contour4=.99,_donorm=False):

    #  READ healpix map    
    hpx1=hpx[:]

    # normlization?
    if _donorm:hpx = hpx/sum(hpx)

    # sorted pixel by probability
    sorted = hpx[hpx.argsort()][::-1]    
    index = hpx.argsort()[::-1]
    cumulative=np.cumsum(sorted)       

    # limit probability
    ilist={}
    for _cc in [_contour1,_contour2,_contour3,_contour4]:       
        if len(sorted[cumulative<_cc])>0:
            limit1 = sorted[cumulative<_cc][-1]             
            index1 = [i for i in index if hpx1[i] >= limit1]# count the area
            ilist[_cc]=index1
        else:ilist[_cc]=[]
    return ilist,hpx 

def gwdist(prob, distmu, distsigma, distnorm,ra,dec):
     
    distmin,distmax,delta,frac = 0,10000,1000,.1
    npix = len(prob)                  
    nside = hp.npix2nside(npix)

    rminlist,rmaxlist = [],[]
    for _ra,_dec in zip(ra,dec):
        theta = 0.5 * np.pi - np.deg2rad(_dec)
        phi = np.deg2rad(_ra)
        ipix = hp.ang2pix(nside, theta, phi)
        r = np.linspace(distmin,distmax,num=delta)   
        dp_dr = r**2 * distnorm[ipix] * scipy.stats.norm(\
                    distmu[ipix], distsigma[ipix]).pdf(r)
        rl = r[np.where(dp_dr>frac*max(dp_dr))]        
        if len(rl)>0: 
            rminlist.append(min(rl))
            rmaxlist.append(max(rl))
        else: 
            rminlist.append(None)
            rmaxlist.append(None)
    return rminlist, rmaxlist

def mcmc(skymap,num,limra=False,limdec=False,fovh=3.,fovw=3.,fig=[1,2,3],radius=False,\
         rot_theta=0.,rot_phi=0.,verbose='2',intractive=False,ordering=False,\
         coord='C',norm = 'hist',color1='grey',color2='k',contour1=.99,contour2=.68):
    
    index1,index2,hpx = contour(skymap,.50,.99)  
    nside=hp.get_nside(hpx)
    dec,ra = IndexToDeclRa(nside,index2)
    if not limra:limra = [min(ra),max(ra)]
    if not limdec:limdec = [min(dec),max(dec)]

    # centering
#    print(np.mean(limra),np.mean(limdec))
#    rot_theta,rot_phi=0,0#np.mean(limra)#np.mean(limra),np.mean(limdec)

    if False:
        # contourview skymap show
        pparams = {'skymap':skymap,'contour1':contour1,'contour2':contour2,\
                   'rot_phi':rot_phi,'rot_theta':rot_theta,\
                   'color1':color1,'color2':color2,'label':'contour','coord':coord}
        optparams = ['rot_theta','rot_phi']
        if intractive:pstplot.intractive_show(pstplot.contourview,pparams,optparams)
        else:pstplot.contourview(pparams)  


    # monte carlo for tiling
    _log=[0.]
    _nloop=1
    for nn in [5,10,15,20]:
        print('searching in fovh/%i fovw/%i'%(nn,nn))
        shifth=fovh/nn
        shiftw=fovw/nn             
        answ = False
        showplot=True
        nloop=0
        while not answ:             
            answ = True
            good_answ=False            
            theta = random.uniform(0,2*np.pi)
            _print='\t with angle: %.2f'%theta
            _shifth,_shiftw = np.sqrt(shifth**2+shiftw**2)*mt.sin(theta),\
                              np.sqrt(shifth**2+shiftw**2)*mt.cos(theta)  
            while not good_answ: 
                good_answ = True                   
                _verbose = (verbose and showplot)                
                if _verbose:
                    # mollview skymap show
                    pparams = {'hpmap':skymap,'title':'skymap','rot_phi':rot_phi,\
                               'rot_theta':rot_theta,'fignum':fig[0],'ordering':ordering,\
                               'coord':coord,'norm':norm}
                    optparams = ['rot_theta','rot_phi']
                    pstplot.mollview(pparams)                     
                
                # generate pointings
                _ral,_decl=scheme.pointings(limdec=limdec,limra=limra,fovh=fovh,fovw=fovw,fig=fig[0],\
                                            shifth=_shifth,shiftw=_shiftw,rot_theta=rot_theta,\
                                            rot_phi=rot_phi,verbose=False,intractive=intractive)
                # cal prob for tiling list
                r2,d2,t = priorization.calprob_tile(skymap,_ral,_decl,fovh,fovw)    

                if _verbose:
                    # highlight selected
                    pparams = {'ra':r2[:num],'dec':d2[:num],'rot_phi':rot_phi,\
                               'rot_theta':rot_theta,'color':'r','fovw':fovw,'fovh':fovh}
                    optparams = ['rot_theta','rot_phi']
                    if intractive:pstplot.intractive_show(pstplot.verticeview,pparams,optparams)
                    else:pstplot.verticeview(pparams)

                if sum(t)>_log[-1]:                    
                    _log.append(sum(t))
                    good_answ=False
                    showplot=True                   

                    print(_print+'\tOK')
#                    plt.figure(fig[0])
                    plt.clf()
                    limra=[limra[0]+_shiftw,limra[1]+_shiftw]
                    limdec=[limdec[0]+_shifth,limdec[1]+_shifth]
                else:
                    nloop+=1
                    showplot=False
                    if nloop<_nloop:
                        answ=False                                                                    
                        print(_print+'\tno')
    
    _oo = open('radec.list','w')
    _ra,_dec,_score = [],[],[]
    print('Priorized pointings:')
    for ii,jj,kk in zip(r2[:num],d2[:num],t[:num]):
        print(ii,jj,kk)
        _ra.append(ii)
        _dec.append(jj)
        _score.append(kk)
        _oo.write('%.6f %.6f %.3f\n'%(ii,jj,kk))
    _oo.close()
    return _ra,_dec,_score

def get_skymap(skymap_url,graceid,_dir):
    """
    Look up URL of sky map in VOEvent XML document,
    download sky map, and parse FITS file.
    """  
    import requests,tempfile,shutil

    if False:
        # Send HTTP request for sky map   
        response = requests.get(skymap_url, stream=True)   

        # Raise an exception unless the download succeeded (HTTP 200 OK)
        response.raise_for_status()

        # Create a temporary file to store the downloaded FITS file
        with tempfile.NamedTemporaryFile() as tmpfile:
            # Save the FITS file to the temporary file
            shutil.copyfileobj(response.raw, tmpfile)
            tmpfile.flush()
      
            # Uncomment to save FITS payload to file       
            shutil.copy(tmpfile.name, _dir + graceid+'_'+os.path.basename(skymap_url))       

    else:
        # wget
        os.system(' '.join(['wget',skymap_url,'-O',_dir + graceid+'_'+os.path.basename(skymap_url)]))

    # Done!
    return _dir + graceid +'_'+os.path.basename(skymap_url)

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

def trigger_validation(skymap, header):

    # validate a trigger by using the header of fits
    ''' judge interests of the trigger '''

    # read and print some values from the FITS header.
    ''' need to be noticed: sometimes failed'''
    header = dict(header)

    # read distance
    try:
        Dmean = header['DISTMEAN']
        Dvar = header['DISTSTD']
        _ld = '%.2f+/-%.2f'%(Dmean,Dvar)
        _str = '#\tDist:\t%s Mpc\n'%_ld
    except:
        _ld = None
        _str = '#\tDist:\tNot available\n'

    # show params of the fits header
    _strremain = ['DATE-OBS', 'CREATOR', 'MJD-OBS']
    for ii in header:
        if ii in _strremain: _str+='#\t%s:\t%s\n'%(ii,header[ii])

    # check the area of skymap
    ilist,hpx = contour(skymap)

    # get area in single pixel
    _areasingle = (hp.nside2resol(hp.get_nside(hpx), arcmin=True)/60.)**2

    # for each contour
    _str += '#\tSky localization:\t'
    _area={}
    for _cc in ilist:
        index1=ilist[_cc]
        if len(index1)>0:
            _str += ' %.2f sq. deg [%i%%]'%(len(index1)*_areasingle,_cc*100)
            _area[_cc]=len(index1)*_areasingle
        else:
            _str += ' NULL [%i%%]'%(_cc*100)
            _area[_cc]=None
    _str+='\n'
    return ilist,_area,_ld,_str,header['DATE-OBS'],header['MJD-OBS']

def read_filelist(flist):

    if flist[0]=='@':             # read file list iraf format
        ff = open(flist[1:])
        _flist = [_ff.replace('\n','') for _ff in ff.readlines() if _ff[0] != '#']       
    elif '*' in flist:
        _flist = sorted(glob.glob(flist))
    else: _flist = flist.split(',')
    _flist = [_ff for _ff in _flist]
    for ff0 in _flist:
        if not os.path.exists(ff0):
            print('!!! ERROR: file '+ff0+' not found !!!')
            sys.exit()
    return _flist

def file_cut(_idc,_ral,_decl,_fskip,_vskip):

    # cut file/file list (,)

    # get center ra,dec of OB
    ra,dec = [],[]
    for _ra,_dec in zip(_ral,_decl):
        ra.append(np.mean(_ra))
        dec.append(np.mean(_dec))    
    _ral0,_decl0,_idc = np.array(ra),np.array(dec),np.array(_idc)

    try: float(_vskip)
    except:sys.exit('Error: vskip!!!')

    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(ra=_ral0*u.deg, dec=_decl0*u.deg)

    for _fskip1 in _fskip.split(','):
        if os.path.exists(_fskip1):                
            for hh in open(_fskip1).readlines():
                if hh[0]=='#':continue
                try:
                    raa,decc = float(hh.split()[0]),float(hh.split()[1])
                    skipradec = astropy.coordinates.SkyCoord(raa*u.deg, decc*u.deg)
                    sep = radecs.separation(skipradec)                        
                    _length = len(radecs[np.where(sep.deg<_vskip)])
                    if _length>0:
                        print('\t-remove %s sources due to %s'%(_length,_fskip1))   
                        _idc = np.delete(_idc,np.where(sep.deg<_vskip))                                     
                        _ral = np.delete(_ral,np.where(sep.deg<_vskip))
                        _decl = np.delete(_decl,np.where(sep.deg<_vskip))
                except:
                    print('!!!Error: format wrong...')
                    input(hh)
        else:print('%s not exists, skipping...'%_fskip1)
    return _idc,_ral,_decl

def ebv_cut(_idc,_ral,_decl,ebvt,verbose,_mode):
    # extinction cut

    # get center ra,dec of OB
    ra,dec = [],[]
    for _ra,_dec in zip(_ral,_decl):
        ra.append(np.mean(_ra))
        dec.append(np.mean(_dec))
    _ral0,_decl0,_idc = np.array(ra),np.array(dec),np.array(_idc)

    if _mode == 'man':
        if len(glob.glob('/'.join(_pstpath.split('/')[:-2])+'/data/ebv_hpmap/*'))>0:
            ebvmapl = []
            for nii,ii in enumerate(glob.glob('/'.join(_pstpath.split('/')[:-2])+'/data/ebv_hpmap/*')):
                print('%i - %s'%(nii,ii))
                ebvmapl.append(ii)
            answ = input('Wanna ebv map instead? which one?')
            try:
                ebvmap = ebvmapl[int(answ)]      
                ''' ebv map '''
                hpx, header = hp.read_map(ebvmap, h=True, verbose=False)    
                _index = DeclRaToIndex(_decl0,_ral0, hp.get_nside(hpx))           
                _id = np.logical_and(_idc,hpx[_index]<ebvt)
                _ral,_decl,_idc = _ral[_id],_decl[_id],_idc[_id]
                return _idc,_ral,_decl             
            except:print('Error: wrong answer, do url approach')

    # url approach, slow, query ebv one by one
    ra0,dec0,id0=[],[],[]
    for _id00,_ra00,_dec00,_ral00,_decl00 in zip(_idc,_ral0,_decl0,_ral,_decl):
        _ebv = query_ebv(_ra00,_dec00,thresh=ebvt,verbose=verbose)
        if _ebv:
            ra0.append(_ral00)
            dec0.append(_decl00)
            id0.append(_id00)
    return np.array(id0),np.array(ra0),np.array(dec0)

def main(mapname,optlist,_mode):

    '''
    for alerts with healpix map: LVC alerts currently
    '''
    notestlist = ['enrico.cappellaro@inaf.it','aniello.grado@inaf.it','fedor.getman@inaf.it','enzo.brocato@inaf.it','paolo.davanzo@inaf.it','massimo.turatto@inaf.it','lina.tomasella@inaf.it','fiore.deluise@inaf.it','stefano.benetti@inaf.it']

    # Main procedure start
    start_time = time.time()

    # try cleaning all the plots
    plt.close('all')

    # define plots setting
    _colorlist,ncolor = ['b','g','k','y','c','m'],0

    # Get JD and time   
    if optlist['arg']['observe']['obstime'] == 'now':
        timenow = astropy.time.Time.now()
        _schsuf = str(timenow).split('T')[0].replace('-','')
    else:
        try: timenow = astropy.time.Time(optlist['arg']['observe']['obstime'], scale='utc')
        except: timenow = astropy.time.Time.now() + astropy.time.TimeDelta(float(optlist['arg']['observe']['obstime'])*60, format='sec')
        _schsuf = str(timenow).split()[0].replace('-','')
    _jdnow = timenow.jd

    ## decide savefig, plot, or no fig
    # pmet: 1.if send email: savefig for specific plots
    #       2.normal verbose: input for all plots
    #       3.inter verbose: inter_plot for all plots
    #       4.no fig
    if _mode=='auto':
        if eval(optlist['arg']['plot']["verbose"]):pmet=1
        else:pmet=4
    if _mode=='man':
        if eval(optlist['arg']['plot']["verbose"]):
            if eval(optlist['arg']['plot']["interactive"]):pmet=3
            else:pmet=2
        else:pmet=4
    if eval(optlist['arg']['email']["sendemail"]) or \
       eval(optlist['arg']['wechat']["activate"]) or \
       eval(optlist['arg']['phone']['activate']):pmet=1
    if pmet in [2,3]:pl.ion()
    if len(optlist['arg']['plot']["showmap"])>0:_showmap = optlist['arg']['plot']["showmap"].split(',')
    else:_showmap = []
        
    ''' 1- Priorization algorithm 
    For all telescopes, strategy, etc'''

    # define a blanket Priorization map
    nside = int(optlist['arg']['priorization']["nside"])
    _pmap = np.zeros(hp.nside2npix(nside))+1

    ''' 1.1 - trigger Priorization '''
    if eval(optlist['arg']['priorization']['trigger']):
        
        # read if exists trigger healpix map
        if not mapname is None:

            # read trigger healpix map            
            try:
                (maptrigger, distmu, distsigma, distnorm), header = \
                            hp.read_map(mapname, field=[0, 1, 2, 3], h=True, \
                            verbose=eval(optlist['arg']['plot']["verbose"]))
            except:
                if eval(optlist['arg']['slack']['activate']):
                    # send message via phone immediately
                    slack_client = SlackClient(optlist['arg']['slack']['slack_bot_token'])
                    for _channel in optlist['arg']['slack']['channel'].split(','):
                        slack_client.api_call("chat.postMessage", channel=_channel,
                                text='new alert, however, wrong healpix fits format!!', \
                                as_user=True)  
                logging.info('wrong fits format for healpy!!') 
                return

            # read and print some values from the FITS header.
            header = dict(header)

            # validate skymap
            _ilist,_alist,_ld,_str,_date_obs,_mjd_obs = trigger_validation(maptrigger, header)
            optlist['arg']['email']['emailcontent'] += _str

            # read areas: 50 68 90 99
            index1,index2,index3,index4=_ilist.values()
            _a50,_a68,_a90,_a99=_alist.values()                      

            # voevent params
            if 'voevent' in optlist['arg'] and header['CREATOR']!='PSTOOLS':
                # for LVC
                params = optlist['arg']['voevent']

                if True:
                    if optlist['arg']['email']['role'] == 'test':_trigger = False
                    else:_trigger = True
                else:_trigger = True
                    
                # for LVC: params of VOEvent to show in email
                _strremain = ['Terrestrial', 'Group', 'EventPage', 'GraceID', \
                              'HasRemnant', 'skymap_png', 'skymap_fits', 'AlertType', \
                              'HasNS', 'BBH', 'BNS', 'Instruments', 'NSBH']
                for ii in params:
                    if ii in _strremain:optlist['arg']['email']['emailcontent']+='#\t%s:\t%s\n'%(ii,params[ii])

                # FAR
                FAR = params['FAR']
                FAR = float(FAR)*360*24*3600 # Hz to per year
                optlist['arg']['email']['emailcontent']+='#\tFAR (per year):\t%.4e\n'%FAR

                # params for alert validation
                grace_id = params['GraceID']
                mapurl = params['skymap_fits']
                BBH = float(params['BBH'])
                BNS = float(params['BNS'])
                HasNS = float(params['HasNS'])
                HasRemnant = float(params['HasRemnant'])
                NSBH = float(params['NSBH'])
                Terrestrial = float(params['Terrestrial'])
                Group = params['Group']
                EventPage = params['EventPage']
                Instruments = params['Instruments']

                #The first selection criteria are about the likelyhood that the event is real and astrophysical. So we will take into consideration only events with a False Alarm Rate (FAR) < 1 event per year. Then we will exclude all the events for which the probability of being a terrestrial event is equal or higher than 90%. 
                
                if FAR < 1 and Terrestrial<0.9:_trigger *= True
                else:_trigger *= False
                
                # !!! WHEN TO TRIGGER ???
                ''' For VST:
                PROB_NS > 0.1
                90% skymap area < 100 deg2
                or
                - PROB_NS < 0.1
                - 90% skymap area < 50 deg2
                - Distance < 100 Mpc
                '''
                if float(HasNS)>.1 and float(HasRemnant)>.1:
                    # interesting for all BNS NS-BH
                    _trigger*=True     
                elif float(HasNS)<.1 and _a90<200:
                    if not _ld is None:
                        Dmean = float(_ld.split('+/-')[0])
                        if Dmean<500:_trigger*=True                
                        else:_trigger*=False
                    else:_trigger*=True    
                else:_trigger*=False

                # mysql insert
                if eval(optlist['arg']['database']['activate']):
                    if _ld is None:Dmean,Dvar='None','None'
                    else:Dmean,Dvar=_ld.split('+/-')[0],_ld.split('+/-')[1]
                    if 'Preliminary' in optlist['arg']['email']['emailcontent']:_stage='Preliminary'
                    elif 'Initial' in optlist['arg']['email']['emailcontent']:_stage='Initial'
                    elif 'Update' in optlist['arg']['email']['emailcontent']:_stage='Update'
                    else:_stage='Null'

                    sqlconn.insert_values('ligoevents',\
                                          {'GraceID':grace_id,\
                                           'JD':_mjd_obs,\
                                           'stage':_stage,\
                                           'type':Group,\
                                           'Dmean':Dmean,\
                                           'Dvar':Dvar,\
                                           'FAR':FAR,\
                                           'BNS':BNS,\
                                           'BBH':BBH,\
                                           'NSBH':NSBH,\
                                           'HasNS':HasNS,\
                                           'HasRemnant':HasRemnant,\
                                           'Terrestrial':Terrestrial,\
                                           'time':str(_date_obs),\
                                           'loc50':str(_a50),\
                                           'loc68':str(_a68),\
                                           'loc90':str(_a90),\
                                           'loc99':str(_a99),\
                                           'detector':Instruments,\
                                           'url':mapurl
                                       })
                  
                # email subject
                if _trigger:  _comments = '[interesting]'
                else:  _comments = '[not interesting]'
                optlist['arg']['email']['emailsub'] += _comments

                # slack
                if eval(optlist['arg']['slack']['activate']):
                    # send message via phone immediately   
                    optlist['arg']['phone']['phonecontent'] += '`%s`\n a %s event detected at %sUT\n *FAR=%.2e per year, HasNS=%s, HasRemnant=%s, localization=%.2f-%.2f-%.2f-%.2f (50-68-90-99) sq deg, Distance=%s Mpc, Terrestrial=%.5e*\nGraceDB: %s\ntelescope scheduler working...'%(grace_id,Group,_date_obs,FAR,HasNS,HasRemnant,_a50,_a68,_a90,_a99,_ld,Terrestrial,EventPage)
                    slack_client = SlackClient(optlist['arg']['slack']['slack_bot_token'])
                    for _channel in optlist['arg']['slack']['channel'].split(','):
                        slack_client.api_call("chat.postMessage", channel=_channel,
                                        text=optlist['arg']['phone']['phonecontent'], as_user=True)

                # send phone message
                if eval(optlist['arg']['phone']['activate']) and \
                   optlist['arg']['email']['role'] == 'observation':                       
                    optlist['arg']['phone']['phonecontent'] += '%s, a %s event detected at %sUT %s: FAR=%.2e per year. HasNS=%s. HasRemnant=%s. localization (90 per)=%.2f sq degree. Distance=%s Mpc. Terrestrial=%.5e. %s'%(grace_id,Group,_date_obs,_comments,FAR,HasNS,HasRemnant,_a90,_ld,Terrestrial,EventPage)
                    for _to in optlist['arg']['phone']['to'].split(','):
                        if eval(optlist['arg']['phone']['activate']):pst.link.phone(optlist['arg']['phone']['account'],optlist['arg']['phone']['token'],optlist['arg']['phone']['from'],_to,optlist['arg']['phone']['phonecontent'])

            elif 'voevent' in optlist['arg'] and header['CREATOR']=='PSTOOLS':

                # for other alerts
                params = optlist['arg']['voevent']
                for ii in params:optlist['arg']['email']['emailcontent']+='#\t%s:\t%s\n'%(ii,params[ii])

            else:
                # update graceid via header: MS190402o -> M328682
                # choose whatever you want
                try:grace_id = header['OBJECT']
                except:grace_id='GWxyz'

            # update nside
            if nside != int(header['NSIDE']):
                nside = int(header['NSIDE'])
                _info = 'Warning: for trigger search, change nside to %i, according to trigger healpix map!'%nside
                if eval(optlist['arg']['plot']["verbose"]):print(_info)
                logging.info(_info)

            # read time of GW
            _jdgw = astropy.time.Time(header['MJD-OBS'], format='mjd') + 2400000.5
            timegw = _jdgw.utc.datetime

            # read distance
            # for offline search, there's no ld available
            try:
                Dmean = header['DISTMEAN']
                Dvar = header['DISTSTD']
                _ld = str(Dmean) + ',' + str(Dvar)    
                _info = 'Distance = %s+/-%s'%(Dmean,Dvar)
            except:
                _ld = None
                _info = 'unmodelled trigger without distance'

            if eval(optlist['arg']['plot']["verbose"]):print(_info)           
            logging.info(_info)

            # for Plotting
            # healpix show
            ''' fignum 1'''
            pparams = {'hpmap':maptrigger,'title':'trigger sky map',\
                       'rot_phi':float(optlist['arg']['plot']["rot_phi"]),\
                       'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'fignum':0,\
                       'ordering':eval(optlist['arg']['plot']["ordering"]),\
                       'coord':optlist['arg']['plot']["coord"],\
                       'norm':str(optlist['arg']['plot']["norm"])}
            optparams = ['rot_theta','rot_phi']

            if pmet==3 and 'trigger' in _showmap:pstplot.interactive_show(pstplot.mollview,pparams,optparams)
            if pmet in [1,2] and 'trigger' in _showmap:
                if all(i for i in _alist.values()):fig_trigger = pstplot.contourview(pparams)
                else:fig_trigger = pstplot.mollview(pparams)
            '''
            if pmet==2 and 'trigger' in _showmap:                
                pparams = {'ra':[255.58],'dec':[-12.4856],'fignum':0,'color':'k','rot_phi':float(optlist['arg']['plot']["rot_phi"]),'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'coord':str(optlist['arg']['plot']["coord"]),'ms':2,'label':'UVOT candidate'}
                pstplot.pointview(pparams)       
                input('trigger map')
            '''

            ## cross map with trigger
            # length = 12*nside**2            
            _pmap = maptrigger
            _info = "%i sec to finish trigger priorization"%int(time.time()-start_time)
            if eval(optlist['arg']['plot']["verbose"]):print(_info)
            logging.info(_info)
    else:grace_id = 'no_astro'    

    ''' 
    1.2 - galaxy Priorization: mass/distance/number/... '''

    _galaxymap,_galaxyobs,_tilingobs,_mcmcobs=0,0,0,0
    # check arglist, decide if need to query galaxy    
    for _mm in [eval(optlist['arg']["priorization"]["mass"]),\
                eval(optlist['arg']["priorization"]["dist"]),\
                eval(optlist['arg']["priorization"]["number"])]:
        _galaxymap+=_mm

    # check optlist, decide if need to query galaxy
    for _tt in optlist:
        if _tt == 'arg':continue        
        if optlist[_tt]["pointings"]["scheduler"] =='G':_galaxyobs+=1
        if optlist[_tt]["pointings"]["scheduler"] == 'T':_tilingobs+=1
        if optlist[_tt]["pointings"]["scheduler"] == 'M':_mcmcobs+=1
    logging.info('galaxymap value=%i; galaxyobs value=%i; tilingobs value=%i; mcmcobs value=%i'%(_galaxymap,_galaxyobs,_tilingobs,_mcmcobs))

    # galaxy mass or number
    if _galaxymap>0 or _galaxyobs>0:
        # read galaxies
        #fignum 2,3,4,5
        _fset=[2,3,4,5]
        limrag,limdecg,limmag = [float(zz) for zz in optlist['arg']['galaxies']["limra"].split(',')],[float(zz) for zz in optlist['arg']['galaxies']["limdec"].split(',')],float(optlist['arg']['galaxies']["limmag"])

        # No distance cut
        mindist,fulldist = 0,2000
        '''
        if _ld is None:
            try:
                mindist,fulldist = float(optlist['arg']['galaxies']["dist"].split(',')[0]),\
                                   float(optlist['arg']['galaxies']["dist"].split(',')[1])
            except:mindist,fulldist = 0,1000
        else:
            _deltat = 1
            _dthre = .1
            fulldist = 1000
            _ddist = np.arange(0,fulldist,_deltat)[np.logical_and(np.array([scipy.stats.norm.cdf(dd,float(Dmean), float(Dvar)) for dd in np.arange(0,fulldist,_deltat)])>_dthre,np.array([scipy.stats.norm.cdf(dd,float(Dmean), float(Dvar)) for dd in np.arange(0,fulldist,_deltat)])<1-_dthre)]
            mindist,fulldist = min(_ddist),max(_ddist)

        logging.info('auto set galaxy distance range: %s-%s'%(mindist,fulldist))
        print('auto set galaxy distance range: %s-%s'%(mindist,fulldist))
        '''

        if len(optlist['arg']['galaxies']["cachefile"])>0:cachefile = '%s.npz'%optlist['arg']['galaxies']["cachefile"]
        else:cachefile=False        

        # query galaxies
        # read full GLADE
        _id0,_gname0,_ra0,_dec0,_mag0,_dist0,fig_g1,fig_g3,fig_g4 = \
                    scheme.galaxies(catname=optlist['arg']['galaxies']["catalog"],\
                    limra=limrag,limdec=limdecg,interactive=eval(optlist['arg']['plot']["interactive"]),\
                    distmin=mindist,distmax=fulldist,absmag=limmag,outfits=False,\
                    size=int(optlist['arg']['galaxies']["size"]),mcolor='k',mcolor2='r',mcolor3='grey',\
                    gcolor='grey',gcolor2='g-',gcolor3='r-',_dir=optlist['arg']['data']["dir"],\
                    rot_theta=float(optlist['arg']['plot']["rot_theta"]),\
                    rot_phi=float(optlist['arg']['plot']["rot_phi"]),\
                    nside=nside,ordering=eval(optlist['arg']['plot']["ordering"]),\
                    coord=optlist['arg']['plot']["coord"],_showmap=_showmap,\
                    norm = str(optlist['arg']['plot']["norm"]),\
                    verbose = eval(optlist['arg']['plot']["verbose"]),
                    pmet=pmet,cachefile=cachefile,ypos1=2,ypos2=-2,\
                    nbindist=fulldist,nbinlums=1,figset=_fset)        

        ### different distance at different directions
        # ra,dec cut
        _indexg = DeclRaToIndex(_dec0,_ra0,nside)
        _indexgin = np.in1d(_indexg, _ilist[.68])
        _id0,_gname0,_ra0,_dec0,_mag0,_dist0 = _id0[_indexgin],_gname0[_indexgin],\
                                               _ra0[_indexgin],_dec0[_indexgin],\
                                               _mag0[_indexgin],_dist0[_indexgin]
        # dist cut
        rminlist,rmaxlist = gwdist(maptrigger, distmu, distsigma, distnorm,_ra0,_dec0)
        _iddc = []
        for _ii in range(len(_id0)):
            _dmin,_dmax,_dgal = rminlist[_ii],rmaxlist[_ii],_dist0[_ii]           
            if _dmin is None or _dmax is None:continue
            if np.logical_and(_dgal>_dmin,_dgal<_dmax):_iddc.append(_ii)
        _iddc = np.array(_iddc)        

        try:
#        if not _iddc is None:           
            _id0,_gname0,_ra0,_dec0,_mag0,_dist0 = _id0[_iddc],_gname0[_iddc],_ra0[_iddc],_dec0[_iddc],_mag0[_iddc],_dist0[_iddc]
        except:
            if eval(optlist['arg']['slack']['activate']):               
                slack_client = SlackClient(optlist['arg']['slack']['slack_bot_token'])
                for _channel in optlist['arg']['slack']['channel'].split(','):
                    slack_client.api_call("chat.postMessage", channel=_channel,
                            text='No Glade galaxies available after dist,ra,dec,mag cut', \
                            as_user=True)  
            _id0,_gname0,_ra0,_dec0,_mag0,_dist0 = None, None, None, None, None, None       

    # define radius for galaxies smoothing
    gradius = float(optlist['arg']['galaxies']["radius"])*2*mt.pi/60/360 # from arcsec to radians

    # convolve galaxy info into score
    if _galaxymap>0:

        _ggc=False
        # galaxy luminosity
        if eval(optlist['arg']['priorization']["mass"]):
            _lums0,_ggc=10**((-1)*(_mag0/2.5)),True
            logging.info('convolve galaxy mass to priorization')

        # galaxy counts
        elif eval(optlist['arg']['priorization']["number"]):
            _lums0,_ggc=np.zeros(len(_ra0))+1,True
            logging.info('convolve galaxy number to priorization')

        if _ggc and not _ra0 is None:
            maplums,fig_g5 = priorization.make_hpfitsmap(_ra0,_dec0,_lums0,_dist0,nside,\
                            optlist['arg']['plot']["coord"],eval(optlist['arg']['plot']["ordering"]),\
                            str(optlist['arg']['plot']["norm"]),eval(optlist['arg']['plot']["verbose"]),\
                            eval(optlist['arg']['plot']["interactive"]),float(optlist['arg']['plot']["rot_phi"]),\
                            float(optlist['arg']['plot']["rot_theta"]),gradius,pmet='4',fig=6)
            ## crossmap
            _pmap *= maplums

        # galaxy distance
        if eval(optlist['arg']["priorization"]["dist"]) and not _dist0 is None:
            #fig num 7
            _dist = priorization.dist_galaxydist(_dist0,_ld)
            mapdist,fig_g6 = priorization.make_hpfitsmap(_ra0,_dec0,_dist,_dist0,nside,\
                            optlist['arg']['plot']["coord"],eval(optlist['arg']['plot']["ordering"]),\
                            str(optlist['arg']['plot']["norm"]),eval(optlist['arg']['plot']["verbose"]),\
                            eval(optlist['arg']['plot']["interactive"]),float(optlist['arg']['plot']["rot_phi"]),\
                            float(optlist['arg']['plot']["rot_theta"]),gradius,pmet='4',fig=7)
            ## cross map
            _pmap *= mapdist
            logging.info('convolve galaxy distance to priorization')

        # done
        print("%i sec finish galaxy priorization"%int(time.time()-start_time))

    # frame plots: sun, moon, horizon, etc
    pst.plot_sky(optlist,timenow,_colorlist)
    #

    ''' 2 - generate pointings       
    For different strategies: one for tiling (if any) and one for galaxy (if any)
    !!!! important:
    for tiling, valid only for multiple telescopes with same or similar FoV
    if not, try to use OB mode to make them have the same FoV
    Otherwise, I stopped'''

    # define net log file
    _netlog,_net,_TOB = {'G':[],'T':[],'M':[]},{},{'w':[],'h':[]}
    for _tt in optlist:        
        if _tt == 'arg':continue
        _net[optlist[_tt]['telescope']['name']]={}

    # define scheduler for each telescope
    for _tt in optlist:        
        if _tt == 'arg':continue

        if optlist[_tt]['pointings']["scheduler"]=='T':
            ''' 2.1 - for tiling mode'''             
            _netlog['T'].append(_tt)
            _TOB['w'].append(int(optlist[_tt]['pointings']['ob'].split('*')[0])*float(optlist[_tt]['telescope']['fovw']))
            _TOB['h'].append(int(optlist[_tt]['pointings']['ob'].split('*')[1])*float(optlist[_tt]['telescope']['fovh']))       

            # if store pointings to database
            if len(optlist[_tt]['pointings']["dbfile"])>0:dbfile=optlist[_tt]['pointings']["dbfile"]+'.npz'
            else:dbfile=False

            ''' fignum 8'''  
            limrap,limdecp = [float(zz) for zz in optlist[_tt]['pointings']["limra"].split(',')],\
                             [float(zz) for zz in optlist[_tt]['pointings']["limdec"].split(',')]                          
          
            _info = 'tel(%s) tiling search with: ra:%s dec:%s ob:%s fov:%s*%s'%(_tt,limrap,limdecp,\
                            optlist[_tt]['pointings']['ob'],optlist[_tt]['telescope']['fovw'],optlist[_tt]['telescope']['fovh'])
            print(_info)
            logging.info(_info)  

            # define procedures
            _iftodo,_iftoread,_iftostore=True,True,True
            if dbfile:
                if os.path.exists(dbfile):
                    if _info == str(np.load(dbfile)['info']):_iftodo,_iftostore=False,False
                    else:
                        os.remove(dbfile)
                        _iftoread = False
                else:_iftoread=False
            else:_iftoread,_iftostore=False,False

            if _iftodo:
                # store to database
                logging.info('Creating pointings')

                # read pointings if any
                if len(optlist[_tt]['pointings']["radecfile"])>0:
                    _ralist,_declist = [],[]
                    for _radecf in open(optlist[_tt]['pointings']["radecfile"]).readlines():
                        _ralist.append(float(_radecf.split()[0]))
                        _declist.append(float(_radecf.split()[1]))
                else:_ralist,_declist = None,None

                # generate pointings            
                _idc,_ral,_decl,fig_np = scheme.pointings(ralist=_ralist,declist=_declist,limdec=limdecp,limra=limrap,fovh=float(optlist[_tt]['telescope']['fovh']),fovw=float(optlist[_tt]['telescope']['fovw']),obx=int(optlist[_tt]['pointings']['ob'].split('*')[0]),oby=int(optlist[_tt]['pointings']['ob'].split('*')[1]),fig=0,shifth=0.,shiftw=0.,rot_theta=float(optlist['arg']['plot']["rot_theta"]),rot_phi=float(optlist['arg']['plot']["rot_phi"]),pmet=4,_obmuilt=False) 
                '''obmuilt: if True, plot OB into different pointings, time consuming;
                Otherwise, plot the full OB''' 

                # skipping pointings due to file
                if len(optlist[_tt]['files']["fskip"])>0: _idc,_ral,_decl = file_cut(_idc,_ral,_decl,optlist[_tt]['files']["fskip"],float(optlist[_tt]['files']["vskip"]))

                # skipping pointings due to ebv
                if eval(optlist['arg']['galaxies']["limebv"]):_idc,_ral,_decl = ebv_cut(_idc,_ral,_decl,float(optlist['arg']['galaxies']["limebv"]),eval(optlist['arg']['plot']["verbose"]),_mode)

                # ranking pointings
                _rac,_decc,_score = priorization.calprob_tile(_pmap,_ral,_decl,float(optlist[_tt]['telescope']['fovw']),float(optlist[_tt]['telescope']['fovh']),int(optlist[_tt]['pointings']['ob'].split('*')[0]),int(optlist[_tt]['pointings']['ob'].split('*')[1]))

            if _iftoread:
                logging.info('Reading pointings')                     
                _idc = np.load(dbfile)['id']
                _rac = np.load(dbfile)['ra']
                _decc = np.load(dbfile)['dec']
                _score = np.load(dbfile)['score']          

            if _iftostore:                
                logging.info('Stored pointings')
                np.savez(dbfile,id=_idc,ra=_rac,dec=_decc,score=_score,info=_info)

        if optlist[_tt]['pointings']["scheduler"]=='G' and not _ra0 is None:
            ''' 2.2 - for galaxy mode'''
            _info = 'tel:%s galaxy search'%_tt
            print(_info)
            logging.info(_info)  
            _netlog['G'].append(_tt)             
            _idc,_rac,_decc,_score = priorization.calprob_gal(_pmap,_ra0,_dec0,_id0,radius=gradius)

        if optlist[_tt]['pointings']["scheduler"]=='M':
            _info = 'tel:%s mcmc tiling with ob:%s fov:%s*%s'%(_tt,optlist[_tt]['pointings']['ob'],optlist[_tt]['telescope']['fovw'],optlist[_tt]['telescope']['fovh'])            
            print(_info)
            logging.info(_info)  
            _netlog['M'].append(_tt)
            ''' 2.2 - for mcmc mode'''
            # not include distance and mass yet        
            if False:
                _rac,_decc,_score = mcmc(_pmap)                        
            logging.info('mcmc: tbd!')

        # store       
        for _key,_value in zip(['ra','dec','score','id'],[_rac,_decc,_score,_idc]):
            try:_net[_tt][_key] = _value
            except:_net[_tt][_key] = None
        print("%i sec finish galaxy priorization"%int(time.time()-start_time))

    if len(_netlog['T'])>0:

        # check size of tiling OB
        for _nOB in _TOB:
            if len(np.unique(_TOB[_nOB]))!=1:
                if _mode == 'man':sys.exit('Error: sth wrong with %s'%_nOB)    
                elif _mode == 'auto':
                    logging.info('Error: sth wrong with nOB')
                    return 'Error: sth wrong with nOB'
            else:
                # define net for tiling
                _TOB[_nOB]=np.unique(_TOB[_nOB])[0]

        _info = 'tiling net: w-%s deg; h-%sdeg'%(_TOB['w'],_TOB['h'])
        logging.info(_info)

    # 
    _loglist = {'T':{},'M':{},'G':{}}
    _scorelist,_numlist = {},{}
    for _approach in _netlog:
        if len(_netlog[_approach])==0:continue
        
        ''' for one scheduler, define one general param file 
        with the first telescope of such scheduler'''

        _tt = _netlog[_approach][0]
        _info = '%s with net of %s'%(_approach,_tt)
        logging.info(_info)

        # pointing list
        _rac,_decc,_idc,_score = _net[_tt]['ra'],\
                                 _net[_tt]['dec'],\
                                 _net[_tt]['id'],\
                                 _net[_tt]['score']
        if _rac is None or _decc is None or _idc is None or _score is None:continue

        # get center ra,dec of OB
        ratl,dectl = [],[]
        for _ra,_dec in zip(_rac,_decc):
            ratl.append(np.mean(_ra))
            dectl.append(np.mean(_dec))
        _rac,_decc = np.array(ratl),np.array(dectl)

        # score normalization
        _score = _score/sum(_score)    

        # for cumshow
        _scorelist[_approach] = _score

        # sort       
        idx=np.argsort(np.asarray(_score))[::-1]
        _rac,_decc,_idc,_score=_rac[idx],_decc[idx],_idc[idx],_score[idx]                

        # auto threshold
        if eval(optlist['arg']['priorization']['trigger']) and not mapname is None:

            if _approach == 'T':
                for _cc in _alist:
                    if _alist[_cc] is None:continue                
                    _np = _alist[_cc]/_TOB['w']/_TOB['h']
                    _index = _ilist[_cc]
                    _info = 'c %i%%: %.2f sq deg ; %.2f pointings'%(_cc*100,_alist[_cc],_np)
                    print(_info)
                    logging.info(_info)

            if True:
                # coverage TBD:
                # if area very small, cover 90%
                # otherwise, cover only 68%          
                if _alist[.5]>30: _sindex = _ilist[.5]
                elif _alist[.68]>30: _sindex = _ilist[.68]
                else: _sindex = _ilist[.9]
            else:
                _sindex = _ilist[.5]

            try:_sindex
            except:
                _covr,_r,_covok = [.68,.9,.99],0,False
                while not _covok:                    
                    try:
                        _sindex = _ilist[_covr[_r]]
                        _covok=True
                    except:
                        _r+=1
                        _sindex = None
                    if _r > 3:_covok,_sindex=True,None

        if not _sindex is None:
            
            if _approach == 'T':#check index inside pointing
                idx=[]
                for _sr in range(len(_rac)):              
                    _pindex = ipix_in_box(_rac[_sr],_decc[_sr],_TOB['w'],_TOB['h'],nside)
                    _cindex = list(set(_pindex)&set(_sindex))                
                    if len(_cindex)>0:idx.append(_sr)                    
            elif _approach == 'G':# check cum prob distribution
                _pindex = DeclRaToIndex(_decc,_rac,nside)                
                _cindex = list(set(_pindex)&set(_sindex))
                idx = []
                for _sr in _cindex:idx.append(np.argwhere(_pindex==_sr)[0][0])

            idx = np.array(idx)        
            if len(idx)>0:
                # cut to threshold  
                _rac,_decc,_idc,_score=_rac[idx],_decc[idx],_idc[idx],_score[idx]
            else:
                print('Warning: Too few pointing inside localization? Very strange, check!!!')
                logging.info('Warning: Too few pointing inside localization? Very strange, check!!!')

                _autoset = False
                for ii in np.arange(0,1,.1):
                    if _autoset:continue
                    if len(_score[np.where(_score>ii)])<10:
                        idx = np.where(_score>ii)                  
                        _autoset = True
                _rac,_decc,_idc,_score = _rac[idx],_decc[idx],_idc[idx],_score[idx]

            if len(_rac)==0:
                print('Error: change arange bin')
                logging.info('Error: change arange bin')
                return
        else:
            print('1 pointing valid? check!!!')
            logging.info('1 pointing valid? anyhow, remain the first one')
            _rac,_decc,_idc,_score = _rac[:4],_decc[:4],_idc[:4],_score[:4]

        print('%i fileds remained'%len(_score))
        
        if False:
            # remain fields
            if _approach == 'T':
                pparams = {'ra':_rac,'dec':_decc,'fignum':0,'color':'k','rot_phi':float(optlist['arg']['plot']["rot_phi"]),'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'fovw':_TOB['w'],'fovh':_TOB['h'],'label':'remianing fields'}
                pstplot.verticeview(pparams)
            elif _approach == 'G':
                pparams = {'ra':_rac,'dec':_decc,'fignum':0,'color':'k','rot_phi':float(optlist['arg']['plot']["rot_phi"]),'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'coord':str(optlist['arg']['plot']["coord"]),'ms':2,'label':'remianing galaxies'}
                pstplot.pointview(pparams)                   
            input('highlight selected')

        #
        _info = '%i fileds remained'%len(_score)
        _numlist[_approach] = len(_score)
        print(_info)
        logging.info(_info)  

        # if no fields
        if len(_score)==0:
            if eval(optlist['arg']['email']['sendemail']):
                for _toaddress in optlist['arg']['email']['emailto'].split(','):
                    
                    # test alerts sent only to me
                    if optlist['arg']['email']['role'] == 'test':
                        if _toaddress in notestlist:continue

                    # send msg
                    _sent = link.sendemail_1(optlist['arg']['email']['email'],\
                                             optlist['arg']['email']['emailpass'],\
                                             optlist['arg']['email']['emailsmtp'],\
                                             optlist['arg']['email']['emailsub'],\
                                             optlist['arg']['email']['email'],\
                                             _toaddress,'Warning: no fields remianed')
                    if _sent:print('>>> Email sent to %s'%_toaddress)
                    else:print('>>> Email sent to %s failed!'%_toaddress)
            else:sys.exit('Error: no fields remianed, change threshold value')
        elif len(_score)<10:
            _info = 'Warning: less than 10 fields remianed!!!'
            print(_info)
            logging.info(_info)  

        ''' 3 - pointing arrangement
        scheuler these pointings for different telescopes, for same pointing approach'''

        # define disctionary for tiling and galaxy
        obslist, tflist, _limnoblist = {},{},{}

        # read params
        for _tt0 in _netlog[_approach]:           

            # define observatory and time
            # Geodetic coordinates of observatory    
            observatory = astropy.coordinates.EarthLocation(lat=float(optlist[_tt0]['telescope']['lat'])*u.deg, \
                                                            lon=float(optlist[_tt0]['telescope']['lon'])*u.deg, \
                                                            height=float(optlist[_tt0]['telescope']['alt'])*u.m)
        
            # time per OB
            ''' !!! if dithering or not; if dither, specify repeat number>0 !!! '''
            _repeat = int(optlist[_tt0]['scheduler']['repeat']) + 1
            tf = (float(optlist[_tt0]['telescope']['exptime'])+float(optlist[_tt0]['telescope']['rottime']))*int(optlist[_tt0]['pointings']['ob'].split('*')[0])*int(optlist[_tt0]['pointings']['ob'].split('*')[1])*_repeat

            # append list
            obslist[_tt0] = observatory
            tflist[_tt0] = tf
            _limnoblist[_tt0] = int(optlist[_tt0]['pointings']['limnob'])

        # arange pointings by using prioritization and visibility
        if optlist['arg']['observe']['limmoon'] == 'auto':limmoon='auto'
        else:limmoon=float(optlist['arg']['observe']['limmoon'])
        fdone1,fdone2,_frac = prob_obs_galaxies(_idc,_rac,_decc,_score,int(optlist['arg']['observe']['order']),int(optlist[_tt]['pointings']['nob']),_limnoblist,obslist,tflist,timenow,float(optlist['arg']['observe']['timelast']),float(optlist['arg']['observe']['limsun']),limmoon,optlist['arg']['observe']['limsolorobject'].split(','),float(optlist['arg']['observe']['limsolorobjectnum']),float(optlist['arg']['observe']['limalt']))

        _loglist[_approach][1] = fdone1
        _loglist[_approach][2] = fdone2
        _loglist[_approach][3] = _frac

        # plot
        # get ra,dec from id
        for _tt0 in _netlog[_approach]:
            _indexf = [np.where(_idc==ii)[0][0] for ii in np.array(fdone2[_tt0])[np.where(fdone2[_tt0])]] # remove False
            _raf = _rac[_indexf]
            _decf = _decc[_indexf]                 
            _scoref = _score[_indexf]
            _idf = _idc[_indexf]    

            if _approach in ['T','M']:
                # fignum = 0 or 1
                pparams = {'ra':_raf,'dec':_decf,'fignum':0,'color':_colorlist[ncolor],'rot_phi':float(optlist['arg']['plot']["rot_phi"]),'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'fovw':_TOB['w'],'fovh':_TOB['h'],'label':'%s tilings'%_tt0}
                optparams = ['rot_theta','rot_phi']
                if pmet==3 and 'trigger' in _showmap:pstplot.interactive_show(pstplot.verticeview,pparams,optparams)
                elif pmet in [1,2] and 'trigger' in _showmap:
                    _fig_trigger = pstplot.verticeview(pparams)
                    if _fig_trigger is None:pass
                    else:fig_trigger=_fig_trigger
                if pmet==2 and 'trigger' in _showmap:input('highlight selected')

            elif _approach == 'G':
                pparams = {'ra':_raf,'dec':_decf,'fignum':0,'color':_colorlist[ncolor],'rot_phi':float(optlist['arg']['plot']["rot_phi"]),'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'coord':str(optlist['arg']['plot']["coord"]),'ms':2,'label':'%s galaxies'%_tt0}
                optparams = ['rot_theta','rot_phi']
                if pmet==3 and 'trigger' in _showmap:pstplot.interactive_show(pstplot.pointview,pparams,optparams)
                elif pmet in [1,2] and 'trigger' in _showmap:
                    _fig_trigger = pstplot.pointview(pparams)
                    if _fig_trigger is None:pass
                    else:fig_trigger=_fig_trigger
                if pmet==2 and 'trigger' in _showmap:input('highlight selected')

            ncolor+=1
            print('trigger show')
            logging.info('trigger show')

    # 2 - cumshow
    for _approach in _netlog:
        if len(_netlog[_approach]) == 0:continue       
        pparams = {'score':_scorelist[_approach],'fignum':99,'number':_numlist[_approach],'highlight':_loglist[_approach][2],'color':_colorlist}
        optparams = ['number']
        if pmet==3 and 'cum' in _showmap:pstplot.interactive_show(pstplot.cumshow,pparams,optparams)
        elif pmet in [1,2] and 'cum' in _showmap:fig_cum = pstplot.cumshow1(pparams) 
    if pmet==2 and 'cum' in _showmap:input('cunshow')
    print('cum show')
    logging.info('cum show')

    # output scheduler name for different telescopes
    # and write email contents as str for them
    fnamelist,_strlist = {},{}

    # for a general email
    _strlist['all'] = ''
    fnamelist['all'] = optlist['arg']['data']['dir']+'_'.join([grace_id,_schsuf,optlist['arg']['data']['schfile']])
    for _approach in _netlog:
        for _tt0 in _netlog[_approach]:
            fnamelist[_tt0] = optlist['arg']['data']['dir']+'_'.join([grace_id,_schsuf,optlist[_tt0]['files']['schfile']])
            _strlist[_tt0]=''

    # for trigger search, write infos about trigger first, into str['all']
    if grace_id != 'nGW':
        # lvc infos
        _strlist['all'] += optlist['arg']['email']['emailcontent']
        _strlist['all'] += '#\ttimenow (UT):\t%s\n'%timenow

    # then, telescope infos  
    for _approach in _netlog:     
        for _tt0 in _netlog[_approach]:
            _strlist[_tt0]+='#telescope used:\n'        
            _tel,_lon,_lat,_height,_fovw,_fovh,_ob,_frac,_exptime,_rottime = optlist[_tt0]['telescope']["name"],\
                                                                             optlist[_tt0]['telescope']["lon"],\
                                                                             optlist[_tt0]['telescope']["lat"],\
                                                                             optlist[_tt0]['telescope']["alt"],\
                                                                             optlist[_tt0]['telescope']["fovw"],\
                                                                             optlist[_tt0]['telescope']["fovh"],\
                                                                             optlist[_tt0]['pointings']["ob"],\
                                                                             _loglist[optlist[_tt0]['pointings']["scheduler"]][3],\
                                                                             optlist[_tt0]['telescope']["exptime"],\
                                                                             optlist[_tt0]['telescope']["rottime"]
            _strlist[_tt0] += '#\ttel: %s\tlon:%s\tlat:%s\theight:%s\tfovw:%s*%s\tOB size:%s\tstrategy:%s(coverage:%.2f)\texposure:%s\tpointing reassignmant:%s\n'%(_tel,_lon,_lat,_height,_fovw,_fovh,_ob,_approach,_frac,_exptime,_rottime)
            # OB info
            _strlist[_tt0]+='#OB infos:\n'
            if _approach=='G':_strlist[_tt0]+='#Name(PGC:GWGC:HyperLEDA:2MASS) RA(J2000) DEC(J2000) Dist(Mpc) absBmag Probability(%) airmass alt az time telescope nOB\n'
            else:_strlist[_tt0]+='#RA(J2000) DEC(J2000) Probability(%) airmass alt az time telescope nOB\n'

    # then, pointing criteria
    _strlist['all'] += '\n#pointing criteria:\n'
    _trigger,_mass,_dist,_number,_limsun,_limmoon,_limsolorobject,_limsolorobjectnum,_limalt,_cov = eval(optlist['arg']['priorization']["trigger"]),eval(optlist['arg']['priorization']["mass"]),eval(optlist['arg']['priorization']["dist"]),eval(optlist['arg']['priorization']["number"]),optlist['arg']['observe']["limsun"],optlist['arg']['observe']["limmoon"],optlist['arg']['observe']["limsolorobject"],optlist['arg']['observe']["limsolorobjectnum"],optlist['arg']['observe']["limalt"],optlist['arg']['observe']["cov"]

    # limmoon for auto setting
    if _limmoon == 'auto':
        if 'T' in str(timenow):_year, _month, _day = str(timenow.value).split('T')[0].split('-')
        else:_year, _month, _day = str(timenow.value).split()[0].split('-')
        _date, _status, _light = moon_phase(int(_month),int(_day),int(_year))               

        ''' limitation for the moon / TBD '''
        if _light>=80: _limmoon = 30
        elif _light<80 and _light>=40: _limmoon = 20
        elif _light<40 and _light>=10: _limmoon = 10
        else:_limmoon = 1

    _strlist['all'] += '#\tsun<%s(deg)\tmoon>%s(deg, %s)\tother solor:%s(>%s deg)\tairmass<%s\tcov threshould:%s\n#\tGW prob:%s\tmass:%s\tdistance:%s\tnumber:%s\n'%(_limsun,_limmoon,_status,_limsolorobject,_limsolorobjectnum,_limalt,_cov,_trigger,_mass,_dist,_number)

    # galaxy infos
    _strlist['all']+='\n#Galaxy selection:\n'
    _cat,_filter,_distrange,_limrag,_limdecg,_limmag = optlist['arg']['galaxies']["catalog"],\
                                                       optlist['arg']['galaxies']["filter"],\
                                                       optlist['arg']['galaxies']["dist"],\
                                                       optlist['arg']['galaxies']["limra"],\
                                                       optlist['arg']['galaxies']["limdec"],\
                                                       optlist['arg']['galaxies']["limmag"]
    _strlist['all']+='#\tcatalog:%s\tfilter:%s\tdistance limit:%s\tra limit:%s\tdec limit:%s\tmagnitude limit:%s\n'%(_cat,_filter,_distrange,_limrag,_limdecg,_limmag)

    # write OB
    for _approach in _netlog:
        # for T/G/M
        for ntel in _netlog[_approach]:
            # for each telescope ntel
            _rac,_decc,_score,_idc = _net[ntel]['ra'],\
                                     _net[ntel]['dec'],\
                                     _net[ntel]['score'],\
                                     _net[ntel]['id']

            # norm
            _score = _score/sum(_score)

            ww = 0 # number of OB
            fdone1 = _loglist[_approach][1][ntel]
            for _llk in np.sort(fdone1.keys()):
                # _llk time after 
                nww=0 # detailed time after
                if not fdone1[_llk][0]:
                    # no visible fields currently
                    _strlist[ntel]+='Null\t\tNull\t\tNull\t\tNull\t\tNull\t\tNull\t\t%s\t%s\tNull\n' % \
                        (timenow+astropy.time.TimeDelta(float(_llk)*3600, format='sec'),ntel)
                    continue
                for _lls in fdone1[_llk]:                    
                    ww+=1
                    _lls1,_lls2,_lls3 = _lls # _lls1: id; _lls2: radecs; _lls3: altaz            
                    _ll = np.where(_idc == _lls1)                    
                    for _ra,_dec in zip(_rac[_ll][0],_decc[_ll][0]):
                        # do dithering and repeatation
                        _dither = float(optlist[ntel]['scheduler']['dither'])/3600.
                        for _nrepeat,_repeat in enumerate(range(int(optlist[ntel]['scheduler']['repeat'])+1)):
                            _timeafter = float(_llk)*3600 + \
                                         nww*(float(optlist[ntel]['telescope']['exptime'])+\
                                         float(optlist[ntel]['telescope']['rottime']))

                            if _approach == 'G':
                                # for galaxy search, add name, mag, dist
                                _gname = np.array(_gname0)[np.logical_and(np.array(_ra0)==_ra,np.array(_dec0)==_dec)][0]
                                _mag = np.array(_mag0)[np.logical_and(np.array(_ra0)==_ra,np.array(_dec0)==_dec)]
                                _dist = np.array(_dist0)[np.logical_and(np.array(_ra0)==_ra,np.array(_dec0)==_dec)]

                                _ntaball,_ntab = 6,int(len(_gname)/8)
                                if _ntab>6:_ntab = 1
                                else:_ntab = _ntaball-_ntab

                                _strlist[ntel] += _gname
                                _strlist[ntel] += '\t'*_ntab
                                _strlist[ntel] += '%.5f\t%.5f\t%.2f\t%.2f\t%.2e\t%.2f\t%.2f\t%.2f\t%.19s\t%s\t%i\n'%\
                                                  (_ra,_dec,_dist,_mag,_score[_ll],_lls3.altaz.secz,_lls3.alt.deg,\
                                                   _lls3.az.deg,timenow+astropy.time.TimeDelta(_timeafter, format='sec'),ntel,ww)

                            elif _nrepeat == 0:
                                _strlist[ntel]+='%.5f\t%.5f\t%.2e\t%.2f\t%.2f\t%.2f\t%.19s\t%s\t%i\n' % \
                                    (_ra,_dec,_score[_ll],_lls3.altaz.secz,_lls3.alt.deg,_lls3.az.deg,\
                                     timenow+astropy.time.TimeDelta(_timeafter, format='sec'),ntel,ww)

                            else:
                                if ntel in ['VST','schmidt']:# for schmidt/VST, there's auto dither mode
                                    pass
                                else:# for telescope without auto dither mode
                                    # dither when repeat > 1
                                    _radt,_decdt = _ra+_dither*np.cos(random.uniform(0,360)),_dec+_dither*np.sin(random.uniform(0,360))
                                    _strlist[ntel]+='%.5f\t%.5f\t%.2e\t%.2f\t%.2f\t%.2f\t%.19s\t%s\t%i\n' % \
                                        (_radt,_decdt,_score[_ll],_lls3.altaz.secz,_lls3.alt.deg,_lls3.az.deg,\
                                         timenow+astropy.time.TimeDelta(_timeafter, format='sec'),ntel,ww)
                            nww+=1

    # save figures
    # figlist
    figlist={}
    try:figlist['tmap']=fig_trigger
    except:pass
    try:figlist['g1']=fig_g1
    except:pass
    try:figlist['g3']=fig_g3        
    except:pass
    try:figlist['g4']=fig_g4        
    except:pass    
    try:figlist['cum']=fig_cum        
    except:pass

    if pmet==1:
        for _fig,_figname in zip(figlist.values(),figlist.keys()):           
            if _fig:               
                _fig.savefig(optlist['arg']['data']["dir"]+grace_id+'_'+_figname+'.png')
                optlist['arg']['email']['images'].append(optlist['arg']['data']["dir"]+grace_id+'_'+_figname+'.png')

    if 'schmidt' in optlist.keys():
        scheduler.schmidt(_strlist['schmidt'],fnamelist['schmidt'].replace('.txt','.OB.txt'),optlist['schmidt']['scheduler']['filter'],int(optlist['schmidt']['scheduler']['repeat']),float(optlist['schmidt']['telescope']['exptime']),float(optlist[ntel]['scheduler']['dither']),grace_id)
        optlist['arg']['email']['files'].append(fnamelist['schmidt'].replace('.txt','.OB.txt'))
        print('schmidt scheduler done')
        logging.info('schmidt scheduler done')
#        _strlist['all']+='#check schmidt OB in %s\n'%fnamelist['schmidt'].replace('.txt','.OB.txt')

    # write schedule
    for _approach in _netlog:
        for _tt0 in _netlog[_approach]:
            fname = fnamelist[_tt0]
            if os.path.exists(fname):os.remove(fname)
            _oo = open(fname,'w')           
            _oo.write(_strlist[_tt0])        
            _oo.close()
            print('>>> Output a scheduler file:%s'%fname)
            logging.info('>>> Output a scheduler file:%s'%fname)
            optlist['arg']['email']['files'].append(fname)

    # for all
    fname = fnamelist['all']
    if os.path.exists(fname):os.remove(fname)
    _oo = open(fname,'w')
    _oo.write(_strlist['all'])
    _oo.write('\n#Attachments:\n')
    #for _approach in _netlog:
    #    for _tt0 in _netlog[_approach]:_oo.write(_strlist[_tt0]) 
    _oo.close()
    print('>>> Output a scheduler file:%s'%fname)    
    logging.info('>>> Output a scheduler file:%s'%fname)
    optlist['arg']['email']['files'].append(fname)

    if eval(optlist['arg']['email']['sendemail']):
        ####### send email
        with open(fnamelist['all'], 'r') as myfile:
            note=myfile.read()                               
            if eval(optlist['arg']['plot']['verbose']):print(note)
            for _toaddress in optlist['arg']['email']['emailto'].split(','):

                # test alerts sent only to me
                if optlist['arg']['email']['role'] == 'test':
                    if _toaddress in notestlist:continue

                _sent = link.sendemail_2(optlist['arg']['email']['email'],optlist['arg']['email']['emailpass'],optlist['arg']['email']['emailsmtp'],optlist['arg']['email']['emailsub'],optlist['arg']['email']['email'],_toaddress,note,optlist['arg']['email']['images'],optlist['arg']['email']['files'])
                if _sent:_emailinfo = '>>> Email sent to %s'%_toaddress                    
                else:_emailinfo = '>>> Email sent to %s failed!'%_toaddress
                print(_emailinfo)
                logging.info(_emailinfo)

    # auto scheduler generator for specific telescopes
    if 'VST' in optlist.keys():                      
        if eval(optlist['VST']['scheduler']['api']):  #and \
#           optlist['arg']['email']['role'] == 'observation':

            _weburl = scheduler.VST(_strlist['VST'],grace_id,optlist['VST']['scheduler']['environment'],optlist['VST']['scheduler']['username'],optlist['VST']['scheduler']['password'],int(optlist['VST']['scheduler']['containerid']),int(optlist['VST']['scheduler']['userpriority']),float(optlist['VST']['scheduler']['airmass']),optlist['VST']['scheduler']['skytransparency'],float(optlist['VST']['scheduler']['fli']),float(optlist['VST']['scheduler']['seeing']),optlist['VST']['scheduler']['filter'],float(optlist['VST']['telescope']['exptime']),float(optlist[ntel]['scheduler']['dither']),int(optlist['VST']['scheduler']['repeat']))
#            _strlist['VST'] += '# !!! check p2 web link:%s\n'%_weburl
            print('VST scheduler done')
            logging.info('VST scheduler done')
#            _strlist['all']+='#check VST OB in: https://www.eso.org/p2demo/home/container/2297733\n'

    if eval(optlist['arg']['wechat']['activate']) and os.path.exists('wxpy.pkl'):

        # wechat robot
        bot = Bot(console_qr=True, cache_path=True)

        # find friend
        if len(optlist['arg']['wechat']['friends'].split(','))>0:
            for _toaddress in optlist['arg']['wechat']['friends'].split(','):          
                if len(bot.friends().search(_toaddress))==0:print('%s: wechat friend found'%_toaddress)
                else:
                    my_friend = bot.friends().search(_toaddress)[0]
                    print('>>> Wechat sent to %s'%_toaddress)

                    # send msg
                    my_friend.send(optlist['arg']['phone']['phonecontent'])

                    # send img
                    for _img in optlist['arg']['email']['images']:my_friend.send_image(_img)

        # find group
        if len(optlist['arg']['wechat']['groups'].split(','))>0:
            for _toaddress in optlist['arg']['wechat']['groups'].split(','):          
                if len(bot.groups().search(_toaddress))==0:print('%s: wechat group found'%_toaddress)
                else:
                    my_group = bot.groups().search(_toaddress)[0]
                    print('>>> Wechat sent to %s'%_toaddress)

                    # send msg
                    my_group.send(optlist['arg']['phone']['phonecontent'])

                    # send img
                    for _img in optlist['arg']['email']['images']:my_group.send_image(_img)

    if eval(optlist['arg']['web']['activate']):       
        _ftransferlist = np.append(optlist['arg']['email']['images'],\
                                   optlist['arg']['email']['files'])
        if len(optlist['arg']['web']['ssh'])>0:            
            _ssh1,_ssh2 = optlist['arg']['web']['ssh'].split('@')
            _user,_pswd = _ssh1.split(':')
            _server,_port = _ssh2.split(':')
            link.sftp_put_dir(_server, int(_port), _user, _pswd, _ftransferlist, optlist['arg']['web']['dir'] + '/' + grace_id)
            print('web done')
            logging.info('web done')

        # slack
        if eval(optlist['arg']['slack']['activate']):
            # send message via phone immediately               
            slack_client = SlackClient(optlist['arg']['slack']['slack_bot_token'])
            for _channel in optlist['arg']['slack']['channel'].split(','):
                slack_client.api_call("chat.postMessage", channel=_channel,
                                      text='GRAWITA cheduler: http://sngroup.oapd.inaf.it/gw/pic/%s/'%grace_id, as_user=True)
