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
import random
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
    return -np.degrees(theta-np.pi/2.),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))
    
def RadecToThetaphi(ra,dec):
    return np.pi/2.-np.radians(dec),np.radians(ra)

def ThataphiToRadec(theta,phi):
    return np.degrees(phi),np.degrees(np.pi/2.-theta)

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

def compute_contours_1(proportions,samples):

    indexlist = {}
    index = samples.argsort()[::-1]
    sorted_samples = samples[index]
    cumulative_samples = np.cumsum(sorted_samples)
    for proportion in proportions:
        if len(sorted_samples[cumulative_samples<proportion]) > 0:
            _limit = sorted_samples[cumulative_samples<proportion][-1]
            indexlist[proportion] = [i for i in index if samples[i] >= _limit]   
    return indexlist

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


##################################################
##################################################
def pointings(tel='VST',limra=[0,360.],limdec=[-20,90],\
              fovh=1.,fovw=1.,shifth=0.,shiftw=0.,\
              verbose=False,cachemode=3,cachefile='',\
              skipfile='',uarea=0.):
    """
    create pointings by tiling sky with specific FoV and:
    - ra,dec cut
    - shift mode
    - rotation mode, TBD
    """   

    start_time = time.time()
    if verbose: 
        print('\ttel:%s, ra:%s, dec:%s, fov:%s*%s, shift:%s %s'%\
              (tel,limra,limdec,fovw,fovh,shiftw,shifth))
    for _dd in [limra,limdec]:
        if len(_dd)==2:pass
        else:sys.exit('check limra,limdec ...')

    # get pointings
    if verbose: print('\tcache mode: %i; cachefile: %s'%(cachemode,cachefile))
    if len(cachefile)>0 and not '.npz' in cachefile: cachefile += '.npz'
    if len(skipfile)>0 and not '.npz' in skipfile: skipfile += '.npz'

    if cachemode in [1,2]:
        if os.path.exists(cachefile):# if exists
            _indict = np.load(cachefile)
            _cs = [ii in _indict for ii in ['tel','limra','limdec',\
                    'fovh','fovw','shifth','shiftw',\
                    'index','ra','dec']]

            if sum(_cs) != len(_cs):
                # missing fields in cachefile
                if verbose: print('\tmissing fields in cachefile')                
                os.remove(cachefile) # clobber

                # download galaxy catalog 
                if verbose: print('\tgenerate pointings')
                tra,tdec = gen_pointings(limra,limdec,\
                            fovh,fovw,shifth,shiftw)
            else:
                if verbose: print('\tcachefile exists, read')
                # read infos
                _tel,_limra,_limdec,_fovh,_fovw,_shifth,_shiftw = \
                    _indict['tel'],_indict['limra'],_indict['limdec'],\
                    _indict['fovh'],_indict['fovw'],_indict['shifth'],\
                    _indict['shiftw']

                # check if it's OK
                if tel==_tel and _fovh==fovh and _fovw==fovw and \
                   _shifth==shifth and _shiftw==shiftw and \
                   limra[0] >= _limra[0] and limra[1] <= _limra[1] and \
                   limdec[0] >= _limdec[0] and limdec[1] <= _limdec[1]:
                    if verbose: print('\tcachefile meets galaxy settings')
                    tra,tdec = _indict['ra'],_indict['dec']
                else: # no, generate pointings
                    if verbose: print('\tcachefile do not meet pointing settings, remove')
                    os.remove(cachefile) # clobber

                    # download galaxy catalog 
                    if verbose: print('\tcreate pointings')
                    tra,tdec = gen_pointings(limra,limdec,\
                                    fovh,fovw,shifth,shiftw)
        else:           
            if verbose: print('\tcachefile not exists, create pointings')
            # create pointings
            tra,tdec = gen_pointings(limra,limdec,\
                                fovh,fovw,shifth,shiftw)

        if cachemode == 1 and \
           not os.path.exists(cachefile):
            if verbose: print('\tstore cachefile')
            # store to cachefile
            np.savez(cachefile,tel=tel,ra=tra,dec=tdec,\
                     limra=limra,limdec=limdec,\
                     fovw=fovw,fovh=fovh,\
                     shiftw=shiftw,shifth=shifth)

    elif cachemode in [3,4]:
        if verbose: print('\tgenerate pointings')
        tra,tdec = gen_pointings(limra,limdec,\
                            fovh,fovw,shifth,shiftw)

        if cachemode == 3 and os.path.exists(cachefile):            
            os.remove(cachefile) # clobber
            if verbose: print('\tremove cachefile')
            # store to cachefile
            if verbose: print('\tstore cachefile')
            np.savez(cachefile,tel=tel,ra=tra,dec=tdec,\
                     limra=limra,limdec=limdec,\
                     fovw=fovw,fovh=fovh,\
                     shiftw=shiftw,shifth=shifth)
    if verbose: print("\t%i pointings generated in %i sec"%(len(tra),int(time.time()-start_time)))
    tra,tdec = pst.remove_fields(skipfile,tra,tdec,\
                    np.zeros(len(tra))+fovw,\
                    np.zeros(len(tra))+fovh,verbose)
    return tra,tdec

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

def gen_pointings(limra=[0,360.],limdec=[-20,90],\
                  fovh=1.,fovw=1.,shifth=0.,shiftw=0.):

    # constrain shift
    if shifth>fovh:shifth-=fovh
    if shiftw>fovw:shiftw-=fovw
    if shifth<-fovh:shifth+=fovh
    if shiftw<-fovw:shiftw+=fovw

    # cut ra,dec
    if limdec:decmin,decmax = min(limdec)+shifth,max(limdec)+shifth
    else:decmin,decmax = -89.,89.
    decrange= np.arange(decmin,decmax,fovh)

    if limra:ramin,ramax = min(limra),max(limra)
    else:ramin,ramax = 0.,360.

    if ramax>360-fovw:ramax=360-fovw
    if ramin<fovw:ramin=fovw
    if decmax>90-fovh:decmax=90-fovh
    if decmin<-90+fovh:decmin=-90+fovh

    # network: ra,dec list
    ralist,declist=[],[]
    for dec in decrange:
        npoint = 360*np.cos(dec*np.pi/180)/fovw
        for nn in np.arange(0,npoint,1):  
            _ra = 360/npoint*nn+shiftw
            if _ra < ramax and _ra > ramin:
                ralist.append(_ra)
                declist.append(dec)                 
    declist,ralist = np.array(declist),np.array(ralist)
    return ralist,declist

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
        print('\tcat:%s, filter:%s, ra:%s, dec:%s, mag:%s, dist:%s'%(catalog,filtro,limra,limdec,limmag,limdist))    
    for _dd in [limra,limdec,limmag,limdist]:
        if len(_dd)==2:pass
        else:sys.exit('check limra,limdec,limmag,limdist...')

    # get galaxies
    if verbose: print('\tcache mode: %i; cachefile: %s'%(cachemode,cachefile))
    if len(cachefile)>0 and not '.npz' in cachefile: cachefile += '.npz'
    if cachemode in [1,2]:
        if os.path.exists(cachefile):# if exists
            _indict = np.load(cachefile)
            _cs = [ii in _indict for ii in ['limra','limdec',\
                                            'limmag','limdist',\
                                            'cat','filtro','name',\
                                            'ra','dec','mag','dist']]
            if sum(_cs) != len(_cs):
                # missing fields in cachefile
                if verbose: print('\tmissing fields in cachefile')                
                os.remove(cachefile) # clobber

                # download galaxy catalog 
                if verbose: print('\tquery vizier')
                gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                            size,filtro,limra,limdec,\
                                            limmag,limdist,verbose)
            else:
                if verbose: print('\tcachefile exists, read')
                # read infos
                _limra,_limdec,_limmag,_limdist,_cat,_filter = \
                                _indict['limra'],\
                                _indict['limdec'],\
                                _indict['limmag'],\
                                _indict['limdist'],\
                                _indict['cat'],\
                                _indict['filtro']            

                # check if it's OK
                if catalog==_cat and str(filtro) == str(_filter) and \
               limra[0] >= _limra[0] and limra[1] <= _limra[1] and \
               limdec[0] >= _limdec[0] and limdec[1] <= _limdec[1] and \
               limmag[0] >= _limmag[0] and limmag[1] <= _limmag[1] and \
               limdist[0] >= _limdist[0] and limdist[1] <= _limdist[1]:
                    if verbose: print('\tcachefile meets galaxy settings')
                    gname,gra,gdec,gmag,gdist = _indict['name'],\
                                                _indict['ra'],\
                                                _indict['dec'],\
                                                _indict['mag'],\
                                                _indict['dist']
                else: # no, generate galaxies
                    if verbose: print('\tcachefile do not meet galaxy settings, remove')
                    os.remove(cachefile) # clobber

                    # download galaxy catalog 
                    if verbose: print('\tquery vizier')
                    gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                                size,filtro,limra,limdec,\
                                                limmag,limdist,verbose)
        else:           
            if verbose: print('\tcachefile not exists, query vizier')
            # download galaxy catalog 
            gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                        size,filtro,limra,limdec,\
                                        limmag,limdist,verbose)
        if cachemode == 1 and \
           not os.path.exists(cachefile):
            if verbose: print('\tstore cachefile')
            # store to cachefile
            np.savez(cachefile,cat=catalog,filtro=filtro,\
                     limra=limra,limdec=limdec,\
                     limmag=limmag,limdist=limdist,\
                     name=gname,ra=gra,dec=gdec,\
                     mag=gmag,dist=gdist)

    elif cachemode in [3,4]:
        if verbose: print('\tquery vizier')
        gname,gra,gdec,gmag,gdist = pst.query_vizier(catalog,\
                                size,filtro,limra,limdec,\
                                limmag,limdist,verbose)
        if cachemode == 3 and os.path.exists(cachefile):            
            os.remove(cachefile) # clobber
            if verbose: print('\tremove cachefile')
            # store to cachefile
            if verbose: print('\tstore cachefile')
            np.savez(cachefile,cat=catalog,filtro=filtro,\
                     limra=limra,limdec=limdec,\
                     limmag=limmag,limdist=limdist,\
                     name=gname,ra=gra,dec=gdec,\
                     mag=gmag,dist=gdist)

    print("\t%i galaxies obtained in %i sec\n"%(len(gra),int(time.time()-start_time)))
    return gname,gra,gdec,gmag,gdist

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
