#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/tilings.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import numpy as np
import healpy as hp
import random
import astropy.units as u
from astropy.table import Table
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
   
        def __init__(self, ra=None, dec=None,
                     fovra=None, fovdec=None,
                     defconf=None, logger=None):
                """
                generate tiling network

                Parameters:
                -----------
                ra          `list`           tiling ra list
                dec         `list`           tiling dec list
                fovra       `list`           tiling field of view list in ra direction
                fovdec      `list`           tiling field of view list in dec direction
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
                self.data   =   {'ra':ra, 'dec':dec,
                                 'fovra':fovra, 'fovdec':fovdec}

                # ----- define default parameters ----- #
                self.run_config(defconf)

        def run_config(self, defconf):
                     
                self.conf = {
                        'limra':      [0,360.],  # `range`   [0,360.]
                        'limdec':     [-90,90],  # `range`   [-90.,90.]
                        'fovra':      1.,        # `float`   > 0
                        'fovdec':     1.,        # `float`   > 0
                        'shiftra':    0.,        # `float`   abs(shiftra) < fovra
                        'shiftdec':   0.,        # `float`   abs(shiftdec) < fovdec
                        'skipfile':   None,      # `str`     options: `txt`, `npz`
                        'skipfrac':   0.,        # `float`   [0, 1]
                        'wdir':       '/tmp/',   # `str`     working directory
                        'savetype':   'npz'      # `str`     options: `txt`, `npz`
                }
        
                if defconf is None: return        
                for _k in self.conf.keys():
                        if _k in defconf.keys():
                                self.conf[_k] = defconf[_k]
                
        def gen_pointings(self,limra=None,limdec=None,
                          fovra=None,fovdec=None,
                          shiftra=None,shiftdec=None,                          
                          skipfile=None, skiparea=None,
                          savetype=None, wdir=None):
                """
                create pointings by tiling sky
                """
                
                if limra is None: limra = self.conf['limra']
                else: limra = limra
                if limdec is None: limdec = self.conf['limdec']
                else: limdec = limdec
                if wdir is None: wdir = self.conf['wdir']
                else: wdir = wdir
                if wdir is None: wdir = self.conf['wdir']
                else: wdir = wdir
                if wdir is None: wdir = self.conf['wdir']
                else: wdir = wdir
                if wdir is None: wdir = self.conf['wdir']
                else: wdir = wdir
                if wdir is None: wdir = self.conf['wdir']
                else: wdir = wdir
                if wdir is None: wdir = self.conf['wdir']
                else: wdir = wdir
                
                # limit shift to proper range
                if abs(shiftra) >= fovra:
                        self.logger.info ('# Warning: abs(shiftra) < fovra')
                        return
                
                if abs(shiftdec) >= fovdec:
                        self.logger.info ('# Warning: abs(shiftdec) < fovdec')
                        return
                
                # cut ra,dec
                ramin, ramax = min(limra)+shiftra, max(limra)+shiftra
                decmin, decmax = min(limdec)+shiftdec, max(limdec)+shiftdec
                ramin = max(ramin, 0)
                ramax = min(ramax, 360)
                decmin = max(decmin, -90)
                decmax = min(decmax, 90)
                
                # dec range
                decrange= np.arange(decmin,decmax,fovdec)

                # get network, i.e. ra,dec list
                ralist,declist=[],[]
                fovralist,fovdeclist=[],[]
                for _dec in decrange:
                        npoint = 360*np.cos(_dec*np.pi/180)/fovra
                        for nn in np.arange(0,npoint,1):  
                                _ra = 360./npoint*nn+shiftra
                                if _ra < ramax and _ra > ramin:
                                        ralist.append(_ra)
                                        declist.append(_dec)
                                        fovralist.append(fovra)
                                        fovdeclist.append(fovdec)  
                self.data = {'ra':ralist, 'dec':declist,
                             'fovra':fovralist, 'fovdec':fovdeclist}

        def read_pointings(self, cachefile):
                self.file2list(cachefile)
                return
        
        def file2list(self, fname):
                return

        def astrotab(self, dic):
                return Table([self.data['ra'],
                              self.data['dec'],
                              self.data['fovra'],
                              self.data['fovdec']],
                             names=('ra', 'dec', 'fovra', 'fovdec'))
        
        def nparray(self):                
                return {'ra':     np.array(self.data['ra']),
                        'dec':    np.array(self.data['dec']),
                        'fovra':  np.array(self.data['fovra']),
                        'fovdec': np.array(self.data['fovdec'])}

        def savez(self, dic):               
                if savetype == 'npz':
                        cachefile = '%s.npz' % save
                        if os.path.exists(cachefile):
                                os.remove(cachefile)

                        # store to cachefile
                        np.savez(cachefile,ra=ralist,dec=declist)
                else:
                        cachefile = '%s.txt' % save
                        if os.path.exists(cachefile):
                                os.remove(cachefile)
                                
                        # store to cachefile
                        np.savez(cachefile,ra=ralist,dec=declist)
        
        def calc_prob(self, hpmap, lcfile=None, limmag=19.):
                return

        def calc_vis(self):
                return
        
#        def priorizatin(self):
#                return
        
        def divide_OB(rac,decc,fovw,fovh,nobw,nobh):
    
                '''
                npoint = 360*np.cos(decc*np.pi/180)/fovw
                _radiff,_decdiff = 360/npoint,fovh
                _nra,_ndec = np.arange(nobw)-(nobw-1)/2,np.arange(nobh)-(nobh-1)/2
                _raspread=rac+_nra*_radiff
                _decspread=decc+_ndec*_decdiff

                ralist,declist=[],[]
                for _ra in _raspread:
                for _dec in _decspread:
                ralist.append(_ra)
                declist.append(_dec)  
                '''

                _ndec = np.arange(nobh)-(nobh-1)/2
                _decdiff = fovh
                _decspread=decc+_ndec*_decdiff

                ralist,declist=[],[]
                for _dec in _decspread:
                        npoint = 360*np.cos(_dec*np.pi/180)/fovw
                        _radiff = 360/npoint
                        _nra = np.arange(nobw)-(nobw-1)/2
                        _raspread=rac+_nra*_radiff

                for _ra in _raspread:      
                        ralist.append(_ra)
                        declist.append(_dec)          
                return ralist,declist

        

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

                v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = vertices(ra, dec, width, height)
                ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra],\
                                             [v1_dec, v2_dec, v4_dec, v3_dec])

                theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
                phi = np.deg2rad(ra_vertices)
                xyz = hp.ang2vec(theta, phi)
                ipix_fov_box = hp.query_polygon(nside, xyz)
                return ipix_fov_box

def compute_contours(proportions,samples):
    r''' Plot containment contour around desired level.
    E.g 90% containment of a PDF on a healpix map
    '''

    levels = []
    sorted_samples = list(reversed(list(sorted(samples))))
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > \
                       proportion).tolist().index(True)
        level = (sorted_samples[level_index] + \
                 (sorted_samples[level_index+1] \
                  if level_index+1 < len(samples) else 0)) / 2.0
        levels.append(level)        

    try: import meander
    except: sys.exit('### Error: install meander via pip...')
    contours_by_level = meander.spherical_contours(sample_points, samples, levels)

    theta_list = {}; phi_list={}
    for cc,contours in enumerate(contours_by_level):
        _cnt = proportions[cc]
        try:theta_list[_cnt]; phi_list[_cnt]
        except: theta_list[_cnt] = []; phi_list[_cnt] = []
        for contour in contours:            
            theta, phi = contour.T
            phi[phi<0] += 2.0*np.pi
            theta_list[_cnt].append(theta)
            phi_list[_cnt].append(phi)
    return theta_list, phi_list


def radec2skycell(rac,decc,fovw,fovh,idlist=None,ralist=None,declist=None,obx=1,oby=1,shifth=0.,shiftw=0.):

    if idlist is None:
        # generate skycells
        cell_list,ralist,declist,fig = pointings(limdec=[-89.,89.],limra=[0,360],fovh=fovh,fovw=fovw,obx=obx,oby=oby,shifth=shifth,shiftw=shiftw)
    else: cell_list,ralist,declist = idlist,ralist,declist

    # which cell
    _netlist = {'ra':[],'dec':[],'cell':[]}
    for _skycell,_ra,_dec in zip(cell_list,ralist,declist):        
        for _subcell in range(len(_ra)):
            _netlist['ra'].append(_ra[_subcell])
            _netlist['dec'].append(_dec[_subcell])
            _netlist['cell'].append('%i-%i'%(_skycell,_subcell))      

    _ralist,_declist = np.array(_netlist['ra']),np.array(_netlist['dec'])
    _dist = np.sqrt((rac-_ralist)**2 + (decc-_declist)**2)
    _index = np.argmin(_dist)
    _skycell,_pra,_pdec = np.array(_netlist['cell'])[_index],_ralist[_index],_declist[_index]
    return _skycell,_pra,_pdec

def skycell2radec(skycell,subcell,fovw,fovh,ralist=None,declist=None,obx=1,oby=1,shifth=0.,shiftw=0.):

    # generate skycells
    cell_list,ralist,declist,fig = pointings(ralist=ralist,declist=declist,limdec=[-89.,89.],limra=[0,360],fovh=fovh,fovw=fovw,obx=obx,oby=oby,shifth=shifth,shiftw=shiftw)

    # which cell
    _netlist = {'ra':[],'dec':[],'cell':[]}
    for _skycell,_ra,_dec in zip(cell_list,ralist,declist):        
        for _subcell in range(len(_ra)):
            _netlist['ra'].append(_ra[_subcell])
            _netlist['dec'].append(_dec[_subcell])
            _netlist['cell'].append('%i-%i'%(_skycell,_subcell))      

    _index = np.where(np.array(_netlist['cell']) == '%i-%i'%(skycell,subcell))
    return np.array(_netlist['ra'])[_index],np.array(_netlist['dec'])[_index]

def remove_fields(skipfile,tra,tdec,fovw,fovh,verbose):
    # remove fields
    if os.path.exists(skipfile):# if exists
        _indict = np.load(skipfile)
        _cs = [ii in _indict for ii in ['ra','dec','fovw','fovh']]
        if sum(_cs) != len(_cs):
            # missing fields in cachefile
            if verbose: print('\tmissing fields in skipfile')
        else:
            _ra,_dec,_fovw,_fovh = \
                _indict['ra'],_indict['dec'],\
                _indict['fovw'],_indict['fovh']

        _tra,_tdec,nn = [],[],0
        for _rap,_decp,_fovwp,_fovhp in zip(tra,tdec,_fovw,_fovh):
            if verbose:
                nn+=1
                _ii = int(float(nn)/len(tra)*100)                
                sys.stdout.write(('='*_ii)+(''*(100-_ii))+("\r [ %d"%_ii+"% ] "))
                sys.stdout.flush()                     
            _cov = pst.overlapregion(_rap,_decp,fovw,fovh,_ra,_dec,_fovwp,_fovhp)
            _frac = _cov/fovwp/fovhp
            if _frac>uarea:
                if verbose:
                    print ('\tskip pointing: ra=%.5f dec=%.5f %s was covered'%\
                           (_rap,_decp,_frac))
                _tra.append(_rap)
                _tdec.append(_decp)
        tra,tdec = np.array(_tra),np.array(_tdec)
    if verbose: print("\t%i pointings remained after removing skip files\n"%(len(tra)))
    return tra,tdec

def pointngsshift(skymap,num,limra=False,limdec=False,fovh=3.,fovw=3.,\
               verbose=False,skipfile=''):

    if not limra and not limdec:
        theta,phi = pst.compute_contours([.99],skymap)
        ra,dec = pst.ThataphiToRadec(theta,phi)
        limra = [min(ra),max(ra)]
        limdec = [min(dec),max(dec)]

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

def overlapregion(ra1,dec1,fovw1,fovh1,ra2,dec2,fovw2,fovh2,nside=1024):
    areasingle = (hp.nside2resol(nside,arcmin=True)/60.)**2      
    index1 = pst.ipix_in_box(ra1,dec1,fovw1,fovh1,nside)
    _ilist = []
    if len(ra2) == 0: return None
    for _ra,_dec,_fovw,_fovh in zip(ra2,dec2,fovw2,fovh2):
        _none = False
        for _ele in [_ra,_dec,_fovw,_fovh]:
            if _ele is None: _none = True
        if _none: continue
        index2 = pst.ipix_in_box(_ra,_dec,_fovw,_fovh,nside)
        for _ii in index2: _ilist.append(_ii)
    _ilist = np.unique(_ilist)
    _cpix = len(list(set(index1)&set(_ilist)))
    return _cpix*areasingle

def overlapregioncut(ra1,dec1,ra2,dec2,radius=100):
    ra2, dec2 = np.array(ra2), np.array(dec2)
    dist = np.sqrt(((ra2-ra1)*15*np.cos(dec1))**2+\
                   (dec2-dec1)**2)
    idx = np.where(dist<=radius)
    return idx
