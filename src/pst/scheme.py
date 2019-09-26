"""############################################################################ 
2019/1/30 Start
1. points - Divide full sky into pointings/tiles
2. galaxies - Select galaxies
3. outat - Generate a galaxy list with number of selected galaxies as weight
""" ############################################################################

# python 2 & 3
from __future__ import print_function
from builtins import input

# Python standard library imports
import sys
import os
import time

# Third-party imports
import numpy as np
import healpy as hp
import math as mt
import astropy.units as u
from astropy.table import Table
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import logging

# self import
import pst
    
def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))
    
def RadecToThetaphi(ra,dec):
    return mt.pi/2.-np.radians(dec),np.radians(ra)

def ThataphiToRadec(theta,phi):
    return np.degrees(phi),np.degrees(mt.pi/2.-theta)

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

def compute_contours(proportions,samples):
    r''' Plot containment contour around desired level.
    E.g 90% containment of a PDF on a healpix map
    '''

    levels, indexlist = [], {}
    # sorted_samples = list(reversed(list(sorted(samples))))    
    index = samples.argsort()[::-1]
    sorted_samples = samples[index]
    cumulative_samples = np.cumsum(sorted_samples)
    nside = hp.pixelfunc.get_nside(samples)
    sample_points = np.array(hp.pix2ang(nside,np.arange(len(samples)))).T
    for proportion in proportions:
        level_index = (np.cumsum(sorted_samples) > \
                       proportion).tolist().index(True)
        level = (sorted_samples[level_index] + \
                 (sorted_samples[level_index+1] \
                  if level_index+1 < len(samples) else 0)) / 2.0
        levels.append(level)        
        _limit = sorted_samples[cumulative_samples<proportion][-1]
        indexlist[proportion] = [i for i in index if samples[i] >= _limit]        

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

    return theta_list, phi_list, indexlist

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

def outcat(ralist,declist,_cat,_fovw=1.,_fovh=1.,_outcat='tmp_fields.asc'):
    _ra,_dec = _cat['RAJ2000'],_cat['DEJ2000']
    _ocat = open(_outcat,'w')
    _ocat.write('#index\tRa\tDec\tNgal\n')
    for nn,(ra,dec) in enumerate(zip(ralist,declist)):
        dx,dy = abs(_ra-ra),abs(_dec-dec)
        _ra,_dec = _ra[np.logical_and(dx<_fovw/2, dy<_fovh/2)],\
                   _dec[np.logical_and(dx<_fovw/2, dy<_fovh/2)]                        
        _ocat.write('%i\t%.2f\t%.2f\t%i\n'%(nn,ra,dec,len(_ra)))
    _ocat.close()   

##################################################
##################################################
def pointings(ralist=None,declist=None,limdec=[-89.,89.],\
              limra=[0,360],fovh=1.,fovw=1.,obx=1,oby=1,shifth=0.,shiftw=0.):
    """
    create pointings by tiling sky with specific FoV and:
    - ra,dec cut
    - shift mode
    - rotation mode, TBD
    """

    start_time = time.time()

    if ralist is None and declist is None:
        # create pointings
        if shifth>fovh*oby or shiftw>fovw*obx:
            print('--shifting outside fov--skippng')
            return

        # cut ra,dec
        if limdec:decmin,decmax = min(limdec)+shifth,max(limdec)+shifth
        else:decmin,decmax = -89.,89.
        decrange= np.arange(decmin,decmax,fovh*oby)

        if limra:ramin,ramax = min(limra),max(limra)
        else:ramin,ramax = 0.,360.

        if ramax>360-fovw*obx:ramax=360-fovw*obx
        if ramin<fovw*obx:ramin=fovw*obx
        if decmax>90-fovh*oby:decmax=90-fovh*oby
        if decmin<-90+fovh*oby:decmin=-90+fovh*oby

        #
        ralist,declist=[],[]
        for dec in decrange:
            npoint = 360*np.cos(dec*np.pi/180)/fovw/obx           
            for nn in np.arange(0,npoint,1):  
                _ra = 360/npoint*nn+shiftw
                if _ra < ramax and _ra > ramin:
                    ralist.append(_ra)
                    declist.append(dec)                 
        declist,ralist = np.array(declist),np.array(ralist)    
        idlist = np.arange(len(ralist))

    pralist,pdeclist=[],[]
    for _ra,_dec in zip(ralist,declist):
        # split OB into pointings              
        _tralist,_tdeclist = divide_OB(_ra,_dec,fovw,fovh,obx,oby)       
        pralist.append(_tralist)
        pdeclist.append(_tdeclist)       
    pidlist = np.arange(len(ralist))

    print("%i pointings generated in %i sec"%(len(ralist),int(time.time()-start_time)))
    return np.array(pidlist),np.array(pralist),np.array(pdeclist)

def galaxies(catalog=1,filtro='B',limra=[0,360.],limdec=[-20,90],\
             limdist=[0,1000.],limmag=[-10.,-18.],size=-1,\
             verbose=False,cachemode=3,cachefile='tmp_glade.npz'):
    """
    Galaxy selection with:
    - catalog:GLADE/GWGC/2MASS/etc
    - ra,dec cut
    - distance cut
    - extinction cut
    - absolute magnitude cut
    - size, size of catalog, -1 for full

    the other parameters are for healpix plotting:
    - rot_theta,rot_phi
    - nside: divide sky into 12*nside**2 pieces with equal area
    - ordering T:nest F:string
    - coord: C,G,E 
       C/equational:ra/dec <-> G/galactic:lon/lat theta/phi
    - norm: healpy color normalization: log,hist,None
    
    the other parameters
    - verbose    
    - cache/cachefile
    - nbindist/nbinlums
    """

    start_time = time.time()
    if verbose: 
        print('cat:%s, filter:%s, ra:%s, dec:%s, mag:%s, dist:%s'%(catalog,filtro,limra,limdec,limmag,limdist))    
    for _dd in [limra,limdec,limmag,limdist]:
        if len(_dd)==2:pass
        else:sys.exit('check limra,limdec,limmag,limdist...')

    # get galaxies
    if verbose: print('cache mode: %i; cachefile: %s'%(cachemode,cachefile))
    if not '.npz' in cachefile: cachefile = cachefile+'.npz'
    if cachemode in [1,2]:
        if os.path.exists(cachefile):# if exists
            _indict = np.load(cachefile)
            _cs = [ii in _indict for ii in ['limra','limdec',\
                                            'limmag','limdist',\
                                            'cat','filtro','name',\
                                            'ra','dec','mag','dist']]
            if sum(_cs) != len(_cs):
                # missing fields in cachefile
                if verbose: print('missing fields in cachefile')                
                os.remove(cachefile) # clobber

                # download galaxy catalog 
                if verbose: print('query vizier')
                gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                            size,filtro,limra,limdec,\
                                            limmag,limdist,verbose)
            else:
                if verbose: print('cachefile exists, read')
                # read infos
                _limra,_limdec,_limmag,_limdist,_cat,_filter = \
                                _indict['limra'],\
                                _indict['limdec'],\
                                _indict['limmag'],\
                                _indict['limdist'],\
                                _indict['cat'],\
                                _indict['filtro']            

                # check if it's OK
                if catalog==_cat and filtro == _filter and \
               limra[0] >= _limra[0] and limra[1] <= _limra[1] and \
               limdec[0] >= _limdec[0] and limdec[1] <= _limdec[1] and \
               limmag[0] >= _limmag[0] and limmag[1] <= _limmag[1] and \
               limdist[0] >= _limdist[0] and limdist[1] <= _limdist[1]:
                    if verbose: print('cachefile meets galaxy settings')
                    gname,gra,gdec,gmag,gdist = _indict['name'],\
                                                _indict['ra'],\
                                                _indict['dec'],\
                                                _indict['mag'],\
                                                _indict['dist']
                else: # no, generate galaxies
                    if verbose: print('cachefile do not meet galaxy settings, remove')
                    os.remove(cachefile) # clobber

                    # download galaxy catalog 
                    if verbose: print('query vizier')
                    gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                                size,filtro,limra,limdec,\
                                                limmag,limdist,verbose)
        else:           
            if verbose: print('cachefile not exists, query vizier')
            # download galaxy catalog 
            gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                        size,filtro,limra,limdec,\
                                        limmag,limdist,verbose)
        if cachemode == 1 and \
           not os.path.exists(cachefile):
            if verbose: print('store cachefile')
            # store to cachefile
            np.savez(cachefile,cat=catalog,filtro=filtro,\
                     limra=limra,limdec=limdec,\
                     limmag=limmag,limdist=limdist,\
                     name=gname,ra=gra,dec=gdec,\
                     mag=gmag,dist=gdist)

    elif cachemode in [3,4]:
        if verbose: print('query vizier')
        gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                size,filtro,limra,limdec,\
                                limmag,limdist,verbose)
        if cachemode == 3 and os.path.exists(cachefile):            
            os.remove(cachefile) # clobber
            if verbose: print('remove cachefile')
            # store to cachefile
            if verbose: print('store cachefile')
            np.savez(cachefile,cat=catalog,filtro=filtro,\
                     limra=limra,limdec=limdec,\
                     limmag=limmag,limdist=limdist,\
                     name=gname,ra=gra,dec=gdec,\
                     mag=gmag,dist=gdist)

    print("%i galaxies obtained in %i sec"%(len(gra),int(time.time()-start_time)))
    return np.arange(len(gra)),gname,gra,gdec,gmag,gdist
