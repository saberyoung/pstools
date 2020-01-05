#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/PstGetGalaxies.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import numpy as np
import healpy as hp
import os
import random
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
from astroquery.vizier import Vizier
import logging

__all__ = ('PstGetGalaxies')

class PstGetGalaxies():
        """
        * implement the following functions:
        -> __init__(self)
        -> 
        """        

        # Static version info
        version = 1.0
   
        def __init__(self, name=None, ra=None, dec=None,
                     distance=None, gname=None, mag=None,
                     index=None, defconf=None, logger=None):
                """
                generate tiling network

                Parameters:
                -----------
                ra          `list`           tiling ra list
                dec         `list`           tiling dec list               
                defconf     `dict`           default configure
                logger      `class`          logger
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set galaxy properties ----- #
                self.data   =   {'ra':ra, 'dec':dec, 'n':index,
                                 'dist':distance, 'name':gname,
                                 'mag':mag}

                # ----- define default parameters ----- #
                self.run_config(defconf)

                # ----- set name for telescope ----- #
                self.set_name(name)
                
        def set_name(self, name):
                
                if name is None:
                        self.name = '%.2f-%.2f-%i' % (self.conf['lon'],
                                self.conf['lat'], self.conf['alt'])
                else:
                        self.name = name

        def run_config(self, defconf):
                     
                self.conf = {
                        'catalog':    'GLADE',
                        'filter':     'B',
                        'size':       -1,
                        'limra':      [0,360.],        # `range`   [0,360.]
                        'limdec':     [-89,89],        # `range`   [-89.,89.]
                        'limdist':    [0,1000],        # `range`   [0,1000]
                        'limmag':     [-18,-99],       # `range`   [-10,-99]                        
                        'nside':      512,             # `int`     default healpix resolution
                        'nest':       False,     # `bool`  healpix map defualt ordering: nest or ring
                        'wdir':       '/tmp/',         # `str`     working directory
                        'filetype':   'npz',           # `str`     options: `txt`, `npz`
                        'filename':   'pst_galaxies',  # `str`     tilings file name
                        'obstime':    None,            #
                        'lat':        24.625,          # la parma
                        'lon':        70.403,          # la parma
                        'alt':        2635             # la parma
                }
        
                if defconf is None: return        
                for _k in self.conf.keys():
                        if _k in defconf.keys():
                                self.conf[_k] = defconf[_k]
                        else:
                                self.logger.info ('### Warning: use default value for %s'%_k)
                
        def generate(self,limra=None,limdec=None,
                     catalog=None,filtro=None,size=None,
                     limdist=None,limmag=None):
                """
                create galaxies by querying Vizier
                """
                
                _hp = self.checkdata()
                if _hp: return
                
                if limra is None:    limra   =  self.conf['limra']
                if limdec is None:   limdec  =  self.conf['limdec']
                if limdist is None:  limdist =  self.conf['limdist']
                if limmag is None:   limmag  =  self.conf['limmag']
                if catalog is None:  catalog =  self.conf['catalog']
                if filtro is None:   filtro  =  self.conf['filter']
                if size is None:     size    =  self.conf['size']                

                # specify columns
                if catalog == 'GLADE':
                        if not filtro in ['B', 'K']:
                                self.logger.info ('### Error: wrong filters for GLADE')
                                return
                        catid, columns = 'VII/281', \
                                ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                                 'Dist', 'PGC', 'GWGC', 'HyperLEDA', '2MASS']
                elif catalog == 'GWGC':
                        if not filtro in ['B']:
                                self.logger.info ('### Error: wrong filters for GWGC')
                                return
                        catid, columns = 'VII/267', \
                                ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                                 'Dist', 'Name']
                else:
                        self.logger.info ('### Error: wrong galaxy catalogs')
                        return

                # download catalog with vizier                
                v = Vizier(columns=columns, \
                           column_filters={columns[0]:'%s..%s'%(str(limra[0]),str(limra[1])),\
                                           columns[1]:'%s..%s'%(str(limdec[0]),str(limdec[1])),\
                                           columns[2]:'%s..%s'%(str(limmag[0]),str(limmag[1])),\
                                           columns[3]:'%s..%s'%(str(limdist[0]),str(limdist[1]))
                           })
                v.ROW_LIMIT = size
                catalogs = v.get_catalogs(catid)[0]

                self.logger.info ("%i galaxies selected from %s"%(len(catalogs),catid))

                # return infos    
                if catalog == 'GLADE':
                        _name = []
                        for ii in range(len(catalogs)):
                                _name.append('%s:%s:%s:%s'%(catalogs[columns[4]][ii], \
                                                            catalogs[columns[5]][ii], \
                                                            catalogs[columns[6]][ii], \
                                                            catalogs['_%s'%columns[7]][ii]))        
                else:
                        _name = catalogs[columns[4]]

                self.data = {'n':       np.arange(len(catalogs[columns[0]])),
                             'name':    np.array(_name),
                             'ra':      np.array(catalogs[columns[0]]),
                             'dec':     np.array(catalogs[columns[1]]),
                             'mag':     np.array(catalogs[columns[2]]),
                             'dist':    np.array(catalogs[columns[3]])}                

        def checkdata(self):
                if not self.data['ra'] is None and \
                   not self.data['dec'] is None and \
                   not self.data['dist'] is None and \
                   not self.data['mag'] is None:
                        self.logger.info ('galaxies has already been parsed')
                        return True
                else:
                        return False

        def hpmap(self, nside=None, nest=None):
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if nside is None: nside = self.conf['nside']
                if nest is None: nest = self.conf['nest']
                
                hpx = np.zeros(hp.nside2npix(nside))
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=nest)
                hpx[ipix] += 10**((-1)*(self.data['mag']/2.5))                
                hpx = hpx/sum(hpx)
                return hpx
        
        def read(self, filename=None, filetype=None, wdir=None):

                _res = self.readfile(filename=filename, filetype=filetype, wdir=wdir)
                if not _res is None:                        
                        self.data = {'ra':_res['ra'], 'dec':_res['dec'],'n':_res['n'],
                                'dist':_res['dist'], 'mag':_res['mag'], 'name':_res['name']}
                
        def readfile(self, filename=None, filetype=None, wdir=None):

                if filename is None: filename = self.conf['filename']   
                if filetype is None: filetype = self.conf['filetype']                
                if not wdir is None: self.conf['wdir'] = wdir

                if filename is None:
                        self.logger.info ('### Warning: filename not defined')
                        return

                _data = {}                
                if filetype is None:     return
                
                elif filetype == 'npz':
                        cachefile = '%s/%s.npz' % (self.conf['wdir'], filename)
                        if os.path.exists(cachefile):
                                for _k in ['ra','dec','dist','mag','n','name']:
                                        try:
                                                _data[_k] = np.load(cachefile)[_k]
                                        except:
                                                self.logger.info ('### Warning: missing keyword %s'%_k)
                                                return                               
                        else:
                                self.logger.info ('### Warning: %s not found'%cachefile)
                                return                        
                
                elif filetype == 'txt':
                        cachefile = '%s/%s.txt' % (self.conf['wdir'], filename)
                        if os.path.exists(cachefile):
                                _datat, _keys = {}, {}
                                _skipl = 1
                                for ll in open(cachefile).readlines():
                                        if ll[0] == '#':
                                                # read header
                                                for _nk,_k in enumerate(ll.replace('#','').split()):
                                                        _keys[_nk] = _k
                                                        _datat[_k] = []
                                        else:
                                                # read data
                                                for _nk,_k in enumerate(ll.split()):
                                                        if _nk in _keys:
                                                                try:
                                                                        _datat[_keys[_nk]].append(float(_k))
                                                                except:
                                                                        self.logger.info ('### Warning: skip line %i'%_skipl)
                                                                        _skipl += 1
                                for _k in ['ra','dec','dist','mag','n','name']:
                                        if _k in _datat.keys():
                                                _data[_k] = _datat[_k]
                                        else:
                                                self.logger.info ('### Warning: %s not found'%_k)
                                                
                        else:
                                self.logger.info ('### Warning: %s not found'%cachefile)
                                return                        
                else:
                        self.logger.info ('### Warning: filetype %s unknown'%filetype)
                        return                

                return _data        
                        
        def astrotab(self):
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                
                return Table([self.data['n'],
                              self.data['name'],
                              self.data['ra'],
                              self.data['dec'],
                              self.data['dist'],
                              self.data['mag']],
                             names=('n', 'name', 'ra', 'dec', 'distance', 'mag'))
        
        def save(self, filetype=None, filename=None, wdir=None):

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return

                if not filetype is None: self.conf['filetype'] = filetype
                if not filename is None: self.conf['filename'] = filename
                if not wdir is None: self.conf['wdir'] = wdir
                
                if self.conf['filetype'] is None:
                        return
                
                elif self.conf['filetype'] == 'npz':
                        cachefile = '%s/%s.npz' % (self.conf['wdir'], self.conf['filename'])
                        if os.path.exists(cachefile): os.remove(cachefile)
                        
                        # store to cachefile                        
                        np.savez(cachefile,ra=self.data['ra'],dec=self.data['dec'],
                                 dist=self.data['dist'],mag=self.data['mag'],
                                 name=self.data['name'],n=self.data['n'])
                                
                elif self.conf['filetype'] == 'txt':
                        cachefile = '%s/%s.txt' % (self.conf['wdir'], self.conf['filename'])
                        if os.path.exists(cachefile): os.remove(cachefile)

                        # store to cachefile
                        ww = open(cachefile,'w')                       
                        ww.write('# n ra dec distance mag name \n')
                        for _ra,_dec,_dist,_mag,_n,_name in zip(self.data['ra'], self.data['dec'],
                                                self.data['dist'],self.data['mag'],
                                                self.data['n'],self.data['name']):
                                ww.write('%i %.5f %.5f %.2f %.2f %s \n'%(_n,_ra,_dec,_dist,_mag,_name))
                        ww.close()
                else:
                        self.logger.info ('### Warning: filetype %s unknown'%self.conf['filetype'])

        def calc_prob_loc(self, triggerobj, nest=None):
                '''
                triggerobj   PstParseTriggers object
                '''
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if not nest is None: self.conf['nest'] = nest
                
                from pst.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                elif is_seq(triggerobj.data['hpmap']):
                        hpx = triggerobj.data['hpmap']
                else: return
                
                probs, ns = [], []
                nside = hp.get_nside(hpx)
                pixarea = hp.nside2pixarea(nside, degrees=True)
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=self.conf['nest'])
                probs = hpx[ipix]/pixarea                        
                return Table([np.arange(len(probs)), probs], names=('n', 'prob'))  
        
        def calc_prob_dis(self, triggerobj, nest=None):
                '''
                triggerobj   PstParseTriggers object
                '''
                                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if not nest is None: self.conf['nest'] = nest
                
                from pst.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                else:
                        return

                from scipy.stats import norm        
                nside = hp.get_nside(hpx)
                pixarea = hp.nside2pixarea(nside, degrees=True)
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=self.conf['nest']) 
                dmu, dsigma, dnorm = hpd1[ipix], hpd2[ipix], hpd3[ipix]
                probs = dnorm * norm(dmu, dsigma).pdf(self.data['dist'])/pixarea
                return Table([self.data['n'], probs], names=('n', 'prob'))

        def altaz(self, obstime=None, lat=None, lon=None, alt=None):

                # observiting time
                if obstime is None: obstime = self.conf['obstime']                
                obstime = self.obstime(obstime)

                # observatory
                if lat is None: lat = self.conf['lat']
                if lon is None: lon = self.conf['lon']
                if alt is None: alt = self.conf['alt']     
                observatory = astropy.coordinates.EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=alt*u.m)
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg, dec=self.data['dec']*u.deg)               

                # Alt/az reference frame at observatory, now
                frame = astropy.coordinates.AltAz(obstime=obstime, location=observatory)

                # Transform grid to alt/az coordinates at observatory, now
                altaz = radecs.transform_to(frame)
                return altaz, obstime, observatory
                
        def calc_airmass(self, obstime=None, lat=None, lon=None, alt=None):

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                return Table([self.data['n'], altaz.secz], names=('n', 'airmass'))  
        
        def calc_sun(self, obstime=None, lat=None, lon=None, alt=None):

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the sun, now?
                sun_altaz = astropy.coordinates.get_sun(obstime).transform_to(altaz)
                return sun_altaz.alt

        def calc_solar(self, sobj, obstime=None, lat=None, lon=None, alt=None):
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no galaxies found')
                        return
                if not sobj in ['moon','mercury','venus','mars','jupiter','saturn','uranus','neptune']:
                        self.logger.info ('sobj should be selected from: '+
                                          'mercury, venus, mars, jupiter, saturn, uranus, neptune')
                        return
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg, dec=self.data['dec']*u.deg)
                
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the solar object, now
                moonloc = astropy.coordinates.get_body(sobj, obstime, observatory).transform_to(altaz)
                return Table([self.data['n'], radecs.separation(moonloc)], names=('n', 'dist'))

        @staticmethod
        def obstime(t=None):
                """ define observing time
                t: None, or float, unit in sec, or astropy.time
                """                
                if t is None:
                        obstime = astropy.time.Time.now()                 
                elif type(t) in [int,float]:
                        obstime = astropy.time.Time.now() + \
                                astropy.time.TimeDelta(t, format='sec')            
                elif type(t) == str:
                        try:
                                obstime = astropy.time.Time(t, scale='utc')
                        except:
                                obstime = astropy.time.Time.now()                
                return obstime
