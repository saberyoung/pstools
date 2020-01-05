#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/tilings.py
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
import logging

__all__ = ('PstGetTilings')

class PstGetTilings():
        """
        * implement the following functions:
        -> __init__(self, ra=None, dec=None,
                    fovra=None, fovdec=None,
                    defconf=None, logger = None)
        -> 
        """        

        # Static version info
        version = 1.0
   
        def __init__(self, name=None, ra=None, dec=None,
                     fovra=None, fovdec=None,
                     index=None, defconf=None,
                     logger=None):
                """
                generate tiling network

                Parameters:
                -----------
                ra          `list`           tiling ra list
                dec         `list`           tiling dec list
                fovra       `list`           tiling field of view list in ra direction
                fovdec      `list`           tiling field of view list in dec direction
                index       `list`           tiling index list
                defconf     `dict`           default configure
                logger      `class`          logger
                """
                        
                # ----- define logger ----- #
                if logger is None:
                        logging.basicConfig(level = logging.INFO)
                        self.logger = logging.getLogger(__name__)
                else:
                        self.logger = logger

                # ----- set tiling properties ----- #
                self.data   =   {'ra':ra, 'dec':dec,'n':index,
                                 'fovra':fovra, 'fovdec':fovdec}

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
                        'limra':      [0,360.],  # `range`   [0,360.]
                        'limdec':     [-89,89],  # `range`   [-89.,89.]
                        'fovra':      1.,        # `float`   > 0
                        'fovdec':     1.,        # `float`   > 0
                        'shiftra':    0.,        # `float`   abs(shiftra) < fovra
                        'shiftdec':   0.,        # `float`   abs(shiftdec) < fovdec
                        'obra':       1,         # `int`
                        'obdec':      1,         # `int`
                        'skipfile':   None,      # `str`     options: `txt`, `npz`
                        'skipfrac':   0.,        # `float`   [0, 1]
                        'nside':      512,       # `int`     default healpix resolution
                        'nest':       False,     # `bool`  healpix map defualt ordering: nest or ring
                        'wdir':       '/tmp/',   # `str`     working directory
                        'filetype':   'npz',     # `str`     options: `txt`, `npz`
                        'filename':   'pst_tilings', # `str`     tilings file name
                                                     #           relative path, use wdir set directory
                                                     #           without suffix
                        'obstime':    None,          #
                        'lat':        24.625,        # la parma
                        'lon':        70.403,        # la parma
                        'alt':        2635           # la parma
                }
        
                if defconf is None: return        
                for _k in self.conf.keys():
                        if _k in defconf.keys():
                                self.conf[_k] = defconf[_k]
                        else:
                                self.logger.info ('### Warning: use default value for %s'%_k)
                
        def generate(self,limra=None,limdec=None,
                     fovra=None,fovdec=None,
                     shiftra=None,shiftdec=None):
                """
                create pointings by tiling sky
                """
                
                _hp = self.checkdata()
                if _hp: return
                
                if not limra is None: self.conf['limra'] = limra
                if not limdec is None: self.conf['limdec'] = limdec
                if not fovra is None: self.conf['fovra'] = fovra
                if not fovdec is None: self.conf['fovdec'] = fovdec
                if not shiftra is None: self.conf['shiftra'] = shiftra
                if not shiftdec is None: self.conf['shiftdec'] = shiftdec                
        
                # limit shift to proper range
                if abs(self.conf['shiftra']) >= self.conf['fovra']:
                        self.logger.info ('# Warning: abs(shiftra) < fovra')
                        return
                
                if abs(self.conf['shiftdec']) >= self.conf['fovdec']:
                        self.logger.info ('# Warning: abs(shiftdec) < fovdec')
                        return
                
                # cut ra,dec
                ramin, ramax = min(self.conf['limra'])+self.conf['shiftra'], \
                        max(self.conf['limra'])+self.conf['shiftra']
                decmin, decmax = min(self.conf['limdec'])+self.conf['shiftdec'], \
                        max(self.conf['limdec'])+self.conf['shiftdec']
                ramin = max(ramin, 0)
                ramax = min(ramax, 360)
                decmin = max(decmin, -89)
                decmax = min(decmax, 89)
                
                # dec range
                decrange= np.arange(decmin,decmax,self.conf['fovdec'])

                # get network, i.e. ra,dec list
                ralist,declist=[],[]
                fovralist,fovdeclist=[],[]
                for _dec in decrange:
                        npoint = 360*np.cos(_dec*np.pi/180)/self.conf['fovra']
                        for nn in np.arange(0,npoint,1):  
                                _ra = 360./npoint*nn+self.conf['shiftra']
                                if _ra < ramax and _ra > ramin:
                                        ralist.append(float('%.5f'%_ra))
                                        declist.append(float('%.5f'%_dec))
                                        fovralist.append(float('%.5f'%self.conf['fovra']))
                                        fovdeclist.append(float('%.5f'%self.conf['fovdec']))
                                        
                self.data = {'n':       np.arange(len(ralist)),
                             'ra':      np.array(ralist),
                             'dec':     np.array(declist),                             
                             'fovra':   np.array(fovralist),
                             'fovdec':  np.array(fovdeclist)}

        def generate_mc(self, skymap, num,
                        limra=None,limdec=None,
                        fovra=None,fovdec=None):
                """
                create pointings by tiling sky
                monte carlo approach to maxmize trigger probability with num pointings
                """
                _hp = self.checkdata()
                if _hp: return
                
                if not limra is None: self.conf['limra'] = limra
                if not limdec is None: self.conf['limdec'] = limdec
                if not fovra is None: self.conf['fovra'] = fovra
                if not fovdec is None: self.conf['fovdec'] = fovdec                

                # monte carlo for tiling
                _log, _nloop = [0.], 3
                shifthi, shiftwi = 0, 0
                for nn in [5.,10.,20.]:
                        if verbose: 
                                print(' - searching in fovh/%i fovw/%i'%(nn,nn))
                        shifth, shiftw=fovh/nn, fovw/nn

                        nloop = 0
                        answ1 = False
                        while not answ1:  # if angle OK: loop 100 times
                                angle = random.uniform(0,2*np.pi)
                                if verbose: print('\t %i with angle: %.2f'%(nloop,angle))
                                _shifth,_shiftw = np.sqrt(shifth**2+shiftw**2)*np.sin(angle),\
                                        np.sqrt(shifth**2+shiftw**2)*np.cos(angle) 
            
                                answ2 = False
                                while not answ2:  # if OK, go on, if no, change
                                        shifthi += _shifth
                                        shiftwi += _shiftw
                                        
                                        # generate pointings
                                        _ral,_decl = pst.gen_pointings(limdec=limdec,limra=limra,\
                                                fovh=fovh,fovw=fovw,shifth=shifthi,shiftw=shiftwi)
                                        # cal prob for tiling list
                                        t = pst.calprob_tile(skymap,_ral,_decl,fovh,fovw)  

                                        # judge converge or not
                                        if sum(t)>_log[-1]:  # accept, direction correct
                                                _log.append(sum(t))
                                                print ('\t\tcovered %.5e probs'%_log[-1])
                                        else:  # reject, change direction
                                                nloop+=1
                                                answ2 = True
                                                if nloop>=_nloop: answ1=True
                                        _ral,_decl = pst.gen_pointings(limdec=limdec,limra=limra,\
                                                fovh=fovh,fovw=fovw,shifth=shifthi,shiftw=shiftwi)
                                        _ral,_decl = pst.remove_fields(skipfile,_ral,_decl,\
                                                                np.zeros(len(limra))+fovw,\
                                                                np.zeros(len(limra))+fovh,verbose)
                return _ral,_decl


        def checkdata(self):
                if not self.data['ra'] is None and \
                   not self.data['dec'] is None and \
                   not self.data['fovra'] is None and \
                   not self.data['fovdec'] is None:
                        self.logger.info ('tilings has already been parsed')
                        return True
                else:
                        return False
                
        def read(self, filename=None, filetype=None, wdir=None):

                _res = self.readfile(filename=filename, filetype=filetype, wdir=wdir)
                if not _res is None:                        
                        self.data = {'ra':_res['ra'], 'dec':_res['dec'],'n':_res['n'],
                                     'fovra':_res['fovra'], 'fovdec':_res['fovdec']}
                
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
                                for _k in ['ra','dec','fovra','fovdec','n']:
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
                                for _k in ['ra','dec','fovra','fovdec','n']:
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
                        self.logger.info ('### Warning: no pointings found')
                        return
                
                return Table([self.data['n'],
                              self.data['ra'],
                              self.data['dec'],
                              self.data['fovra'],
                              self.data['fovdec']],
                             names=('n', 'ra', 'dec', 'fovra', 'fovdec'))                 
               

        def remove_fileds_coo(self, ra, dec, fovra, fovdec,
                              nside=None, nest=None, skipfrac=None):
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                
                if not skipfrac is None: self.conf['skipfrac'] = skipfrac
                if not nside is None: self.conf['nside'] = nside
                if not nest is None: self.conf['nest'] = nest
               
                _frac = self.overlapregion(self.data['ra'],self.data['dec'],
                                           self.data['fovra'],self.data['fovdec'],
                                           ra,dec,fovra,fovdec,self.conf['nside'],
                                           self.conf['nest'])               
                _idx = np.where(np.array(_frac) <= self.conf['skipfrac'])
                self.data['ra']     =   self.data['ra'][_idx]
                self.data['dec']    =   self.data['dec'][_idx]
                self.data['fovra']  =   self.data['fovra'][_idx]
                self.data['fovdec'] =   self.data['fovdec'][_idx]
                self.data['n']      =   self.data['n'][_idx]
        
        def remove_fields_file(self, filename=None, filetype=None, wdir=None,
                               nside=None, skipfrac=None):

                if filename is None: filename = self.conf['filename']  
                if filetype is None: filetype = self.conf['filetype']                
                if wdir is None: wdir = self.conf['wdir']
                if nside is None: nside = self.conf['nside']  
                if skipfrac is None: skipfrac = self.conf['skipfrac']
                
                _res = self.readfile(filename=filename, filetype=filetype, wdir=wdir)
                if _res is None: return

                # tilings should be skipped
                ra,dec,fovra,fovdec = _res['ra'], _res['dec'], _res['fovra'], _res['fovdec']

                # remove coolist
                self.remove_fileds_coo(ra, dec, fovra, fovdec,
                                       nside=nside, skipfrac=skipfrac)
        
        def save(self, filetype=None, filename=None, wdir=None):

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
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
                                 fovra=self.data['fovra'],fovdec=self.data['fovdec'],
                                 n=self.data['n'])
                                
                elif self.conf['filetype'] == 'txt':
                        cachefile = '%s/%s.txt' % (self.conf['wdir'], self.conf['filename'])
                        if os.path.exists(cachefile): os.remove(cachefile)

                        # store to cachefile
                        ww = open(cachefile,'w')                       
                        ww.write('# ra dec fovra fovdec n \n')
                        for _ra,_dec,_fra,_fdec,_n in zip(self.data['ra'], self.data['dec'],
                                        self.data['fovra'],self.data['fovdec'],self.data['n']):
                                ww.write('%i %.5f %.5f %.2f %.2f \n'%(_n,_ra,_dec,_fra,_fdec))
                        ww.close()
                else:
                        self.logger.info ('### Warning: filetype %s unknown'%self.conf['filetype'])

        def calc_prob_loc(self, triggerobj, nest=None):
                '''
                triggerobj   PstParseTriggers object
                '''
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if nest is None: nest = self.conf['nest']
                
                from pst.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                elif is_seq(triggerobj.data['hpmap']):
                        hpx = triggerobj.data['hpmap']
                else: return
                
                probs, ns = [], []
                nside = hp.get_nside(hpx)
                for _ra,_dec,_fovw,_fovh,_n in zip(self.data['ra'], self.data['dec'],
                                self.data['fovra'], self.data['fovdec'], self.data['n']):
                        ipix_poly=(self.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,nest))
                        _probs = hpx[ipix_poly].sum()
                        probs.append(_probs)
                        ns.append(_n)
                return Table([ns, probs], names=('n', 'prob'))  
        
        def calc_prob_dis(self, triggerobj, nest=None, limdist=400):
                '''
                triggerobj   PstParseTriggers object
                '''
                                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if nest is None: nest = self.conf['nest']
                
                from pst.cookbook import is_seq, is_seq_of_seq
                if is_seq_of_seq(triggerobj.data['hpmap']):
                        (hpx, hpd1, hpd2, hpd3) = triggerobj.data['hpmap']
                else:
                        return

                from scipy.stats import norm
                # decide limiting distance                
                r = np.linspace(0, limdist)

                nside = hp.get_nside(hpx)
                pixarea = hp.nside2pixarea(nside, degrees=True)
                theta, phi = np.pi/2.-np.radians(self.data['dec']),np.radians(self.data['ra'])
                ipix = hp.ang2pix(nside,theta,phi,nest=nest)                
                dmu, dsigma, dnorm = hpd1[ipix], hpd2[ipix], hpd3[ipix]
                probl = [dnorm * norm(dmu, dsigma).pdf(rr)/pixarea for rr in r]
                probs = [max(ii) for ii in list(map(list, zip(*probl)))]               
                return Table([self.data['n'], probs], names=('n', 'prob'))  

        def calc_prob_mass(self, galaxyobj, nside=None, nest=None):
                '''
                galaxyobj   PstGetGalaxies object
                '''
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if nside is None: nside = self.conf['nside']
                if nest is None: nest = self.conf['nest']
                
                hpx = galaxyobj.hpmap(nside=nside)
                probs, ns = [], []               
                for _ra,_dec,_fovw,_fovh,_n in zip(self.data['ra'], self.data['dec'],
                                self.data['fovra'], self.data['fovdec'], self.data['n']):
                        ipix_poly=(self.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,nest))
                        _probs = hpx[ipix_poly].sum()
                        probs.append(_probs)
                        ns.append(_n)
                return Table([ns, probs], names=('n', 'prob'))        
                
        def lightcurve(self, lcfile=None):
                return

        def OACAPI(self):
                
                return

        def altaz(self, obstime=None, lat=None, lon=None, alt=None):

                # observiting time
                if obstime is None: obstime = self.conf['obstime']                
                obstime = self.obstime(obstime)

                # observatory
                if not lat is None: self.conf['lat'] = lat
                if not lon is None: self.conf['lon'] = lon
                if not alt is None: self.conf['alt'] = alt
                observatory = astropy.coordinates.EarthLocation(lat=self.conf['lat']*u.deg,
                                lon=self.conf['lon']*u.deg, height=self.conf['alt']*u.m)
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg,
                                                      dec=self.data['dec']*u.deg)               

                # Alt/az reference frame at observatory, now
                frame = astropy.coordinates.AltAz(obstime=obstime, location=observatory)

                # Transform grid to alt/az coordinates at observatory, now
                altaz = radecs.transform_to(frame)
                return altaz, obstime, observatory
                
        def calc_airmass(self, obstime=None, lat=None, lon=None, alt=None):

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                return Table([self.data['n'], altaz.secz], names=('n', 'airmass'))  
        
        def calc_sun(self, obstime=None, lat=None, lon=None, alt=None):

                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the sun, now?
                sun_altaz = astropy.coordinates.get_sun(obstime).transform_to(altaz)
                return sun_altaz.alt

        def calc_solar(self, sobj, obstime=None, lat=None, lon=None, alt=None):
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                if not sobj in ['moon','mercury','venus','mars','jupiter','saturn','uranus','neptune']:
                        self.logger.info ('sobj should be selected from: '+
                                          'mercury, venus, mars, jupiter, saturn, uranus, neptune')
                        return
                
                # ra dec of all fields
                radecs = astropy.coordinates.SkyCoord(ra=self.data['ra']*u.deg,
                                                      dec=self.data['dec']*u.deg)
                
                altaz, obstime, observatory = self.altaz(obstime=obstime, lat=lat, lon=lon, alt=alt)
                
                # Where is the solar object, now
                moonloc = astropy.coordinates.get_body(sobj, obstime, observatory).transform_to(altaz)
                return Table([self.data['n'], radecs.separation(moonloc)], names=('n', 'dist'))        
        
        def divide_OB(self, nobw=None, nobh=None):
                
                # nobw, number of pointings in OB, in ra direction
                # nobh, number of pointings in OB, in dec direction
                
                _hp = self.checkdata()
                if not _hp:
                        self.logger.info ('### Warning: no pointings found')
                        return
                
                if not nobw is None: self.conf['obra'] = nobw
                if not nobh is None: self.conf['obdec'] = nobh                        

                ralist, declist, fovralist, fovdeclist, oblist = [], [], [], [], []
                _nn = 0
                for _ra,_dec,_fovw,_fovh in zip(self.data['ra'], self.data['dec'],
                                                self.data['fovra'], self.data['fovdec']):
                        _ra0,_dec0,_fovw0,_fovh0 = self.divide_OB_one(_ra, _dec, _fovw, _fovh,
                                                        self.conf['obra'], self.conf['obdec'])
                        _nn+=1
                        for _ra00,_dec00,_fovw00,_fovh00 in zip(_ra0,_dec0,_fovw0,_fovh0):
                                ralist.append(float('%.5f'%_ra00))
                                declist.append(float('%.5f'%_dec00))
                                fovralist.append(float('%.5f'%_fovw00))
                                fovdeclist.append(float('%.5f'%_fovh00))
                                oblist.append(_nn)
                        
                self.data = {'ra':     np.array(ralist),
                             'dec':    np.array(declist),
                             'n':      np.array(oblist),
                             'fovra':  np.array(fovralist),
                             'fovdec': np.array(fovdeclist)}                

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

        @staticmethod
        def divide_OB_one(rac, decc, fovw, fovh, nobw, nobh):
                """divide a pointing to a list of sub-pointings
                """
                _ndec = np.arange(nobh)-(nobh-1)/2               
                _decdiff = fovh
                _decspread=decc+_ndec*_decdiff

                ralist,declist, fovralist,fovdeclist = [],[],[],[]
                for _dec in _decspread:
                        npoint = 360*np.cos(_dec*np.pi/180)/fovw
                        _radiff = 360/npoint
                        _nra = np.arange(nobw)-(nobw-1)/2
                        _raspread=rac+_nra*_radiff
                
                        for _ra in _raspread:      
                                ralist.append(_ra)
                                declist.append(_dec)
                                fovralist.append(fovw/nobw)
                                fovdeclist.append(fovh/nobh)
                return  ralist,declist, fovralist,fovdeclist 
        
        @staticmethod
        def ipix_in_box(ra,dec,width,height,nside,nest):
                """Finding the healpix indices of a given box
                """
                v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = \
                        PstGetTilings.vertices(ra, dec, width, height)
                ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                             [v1_dec, v2_dec, v4_dec, v3_dec])                
                theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
                phi = np.deg2rad(ra_vertices)
                xyz = hp.ang2vec(theta, phi)
                ipix_fov_box = hp.query_polygon(nside, xyz, nest=nest)
                return ipix_fov_box

        @staticmethod
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
                        v_dec = np.rad2deg(np.arcsin((np.sin(dec_rad)+\
                                j*np.cos(dec_rad))/(1+i**2+j**2)**0.5))
                        vert_ra.append(v_ra)
                        vert_dec.append(v_dec)
                return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], \
                        vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]

        @staticmethod
        def overlapregion(ra1,dec1,fovw1,fovh1,ra2,dec2,fovw2,fovh2,nside,nest):
                # ra1,dec1,fovw1,fovh1: input tilings
                # ra2,dec2,fovw2,fovh2: remove tilings
                # nside: for calculating area, unit in sq. deg
                
                from pst.cookbook import is_seq
                
                if ra1 is None or ra2 is None:   return None                
                if not is_seq(ra1) or not is_seq(ra2): return None
                               
                areasingle =  hp.nside2pixarea(nside, degrees=True)                
                index1, index2 = [],[]                
                for _index,_ral,_decl,_fovwl,_fovhl in zip([index1,index2],
                                [ra1,ra2],[dec1,dec2],[fovw1,fovw2],[fovh1,fovh2]):
                        for _ra,_dec,_fovw,_fovh in zip(_ral,_decl,_fovwl,_fovhl):
                                _idx = PstGetTilings.ipix_in_box(_ra,_dec,_fovw,_fovh,nside,nest)
                                _index.append(_idx)

                # all slices should be skipped
                index2 = [j for i in index2 for j in i]

                return [len(set(index1[ii]) & set(index2))*areasingle/fovw1[ii]/fovh1[ii]
                        for ii in range(len(index1))] 
