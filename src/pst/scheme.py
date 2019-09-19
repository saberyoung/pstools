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
import astropy.units as u
from astropy.table import Table
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import logging

# self import
from pst import pstdef,pstplot
    
def IndexToDeclRa(NSIDE,index):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))
    
def RadecToThetaphi(ra,dec):
    return mt.pi/2.-np.radians(dec),np.radians(ra)

def ThataphiToRadec(theta,phi):
    return np.degrees(phi),np.degrees(mt.pi/2.-theta)

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
def pointings(ralist=None,declist=None,limdec=[-89.,89.],limra=[0,360],fovh=1.,fovw=1.,obx=1,oby=1,fig=4,shifth=0.,shiftw=0.,rot_theta=0.,rot_phi=0.,verbose=False,interactive=False,pmet=False,_obmuilt=False):

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

    # plot
    colorlist,ii = ['r','g','b','k','y'],0
    if _obmuilt:
        for ralist,declist in zip(pralist,pdeclist):
            ii+=1
            pparams = {'ra':ralist,'dec':declist,'fignum':fig,'label':'tiles','rot_phi':rot_phi,\
                       'rot_theta':rot_theta,'color':colorlist[ii%5],'fovw':fovw,'fovh':fovh}
            optparams = ['rot_theta','rot_phi']

            if pmet==3:pstplot.interactive_show(pstplot.verticeview,pparams,optparams)
            elif pmet in [1,2]:fig = pstplot.verticeview(pparams)          
    else:
        ralist2,declist2 = [],[]
        for _ra,_dec in zip(ralist,declist):
            ralist2.append(np.mean(_ra))
            declist2.append(np.mean(_dec))
        pparams = {'ra':ralist2,'dec':declist2,'fignum':fig,'label':'tiles','rot_phi':rot_phi,\
                   'rot_theta':rot_theta,'color':colorlist[ii%5],'fovw':fovw*obx,\
                   'fovh':fovh*oby}
        optparams = ['rot_theta','rot_phi']

        if pmet==3:pstplot.interactive_show(pstplot.verticeview,pparams,optparams)
        elif pmet in [1,2]:fig = pstplot.verticeview(pparams)

    print("%i pointings generated in %i sec"%(len(ralist),int(time.time()-start_time)))
    if pmet==2:input('pointings')

    return np.array(pidlist),np.array(pralist),np.array(pdeclist),fig

def galaxies(catname='GLADE',outfits=False,limra=[0,360.],limdec=[-20,90],\
             distmin=0,distmax=1000.,absmag=-18.,size=-1,mcolor='k',mcolor2='r',mcolor3='grey',\
             gcolor='grey',gcolor2='g-',gcolor3='r-',_dir='',_showmap=[],\
             rot_theta=0.,rot_phi=0.,nside=1024,ordering=False,\
             coord='C',norm = 'hist',verbose=False,pmet=False,interactive=False,\
             cachefile='tmp_glade.npz',ypos1=2,ypos2=-2,nbindist=80,\
             nbinlums=1,figset=[11,12,13,14]):
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
    - interactive
    - cache/cachefile
    - nbindist/nbinlums
    """

    start_time = time.time()

    ramin,ramax = min(limra),max(limra)
    decmin,decmax = min(limdec),max(limdec)
    _ginfo = '%s : [%.2f,%.2f]\t%s : [%.2f,%.2f]\t%s : < %.2f\t%s : [%.2f,%.2f]'%('ra',ramin,ramax,'dec',decmin,decmax,'mag',absmag,'dist',distmin,distmax)
    if verbose:print(_ginfo)

    #
    _iftodo,_iftoread,_iftostore=True,True,True
    if cachefile:
        if os.path.exists(cachefile):
            if _ginfo == str(np.load(cachefile)['info']):_iftodo,_iftostore=False,False
            else:
                #os.remove(cachefile)
                _iftoread = False
        else: _iftoread=False
    else:_iftoread,_iftostore=False,False    

    if _iftodo:
        _info = 'create %s for galaxies'%(cachefile)
        if verbose:print(_info)
        logging.info(_info)

        # download galaxy catalog 
        _gname,ra0,dec0,mag0,dist0 = pstdef.query_vizier(catname,size,interactive,ramin,ramax,decmin,decmax,absmag,distmin,distmax) 

    if _iftoread:
        if verbose:print('read %s for galaxies'%(cachefile))

        # if cache, store a tempfile, in order to save time next time 
        # read: ra,dec,absmag,dist
        _gname,ra0,dec0,mag0,dist0 = np.load(cachefile)['galaxyid'],\
                                     np.load(cachefile)['ra'],\
                                     np.load(cachefile)['dec'],\
                                     np.load(cachefile)['mag'],\
                                     np.load(cachefile)['dist']
        _info = 'read %i galaxies from %s: %s'%(len(ra0),cachefile,_ginfo)
        if verbose:print(_info)
        logging.info(_info)

    if _iftostore:
        _info = 'store %s for galaxies'%(cachefile)
        if verbose:print(_info)
        logging.info(_info)
        
        np.savez(cachefile,catname=catname,info=_ginfo,galaxyid=_gname,ra=ra0,dec=dec0,mag=mag0,dist=dist0)

    '''
    # create healpix map if u want
    galpixels_Range_num= np.zeros(hp.nside2npix(nside))
    pix_num_Range    = (pstdef.DeclRaToIndex(dec0,ra0,nside))
    galpixels_Range_num[pix_num_Range]+=1
    galpixels_Range_num=hp.sphtfunc.smoothing(galpixels_Range_num,sigma=0.1)
    if outfits:hp.write_map(outfits,galpixels_Range_num,nest=ordering,coord=coord,overwrite=True)
    '''

    # show galaxy properties       
    # 1
    fig_g1,fig_g2,fig_g3,fig_g4 = False,False,False,False

    pparams = {'ra':ra0,'dec':dec0,'rot_phi':rot_phi,'fignum':figset[0],\
               'rot_theta':rot_theta,'color':'r','coord':coord,'label':'full','ms':4}
    optparams = ['rot_theta','rot_phi']
    if pmet==3 and 'galaxy' in _showmap:pstplot.interactive_show(pstplot.pointview,pparams,optparams)
    elif pmet in [1,2] and 'galaxy' in _showmap:fig_g1 = pstplot.pointview(pparams) 
    if pmet==2 and 'galaxy' in _showmap:input('galaxy1')

    '''
    # 2
    pparams = {'hpmap':galpixels_Range_num,'title':_ginfo,'rot_phi':rot_phi,\
               'rot_theta':rot_theta,'fignum':figset[1],'ordering':ordering,\
               'coord':coord,'norm':norm}
    optparams = ['rot_theta','rot_phi']
    if pmet==3:pstplot.interactive_show(pstplot.mollview,pparams,optparams)
    elif pmet in [1,2]:fig_g2 = pstplot.mollview(pparams)
    if pmet==2:input('galaxy2')
    '''

    # 3               
    if distmin==np.nan and distmax==np.nan:
        sys.exit('input distmin & distmax')
    else:           
        pparams = {'distmin':distmin,'distmax':distmax,'distfull_mat':dist0,\
                   'dist_mat':dist0,'color1':mcolor,'color2':mcolor2,'color3':mcolor3,\
                   'nbin':nbindist,'label1':'full','label2':_ginfo,'fignum':figset[2],\
                   'ypos1':ypos1,'ypos2':ypos2}
        optparams = ['distmin','distmax','nbin']
        if pmet==3 and 'galaxy' in _showmap:pstplot.interactive_show(pstplot.distview,pparams,optparams)
        elif pmet in [1,2] and 'galaxy' in _showmap:fig_g3 = pstplot.distview(pparams)    
        if pmet==2 and 'galaxy' in _showmap:input('galaxy3')

    # 4              
    if distmin==np.nan and distmax==np.nan:
        sys.exit('input distmin & distmax')
    else:       
        pparams = {'distmin':distmin,'distmax':distmax,'magfull_mat':mag0,\
                   'distfull_mat':dist0,'dist_mat':dist0,'fignum':figset[3],\
                   'mag_mat':mag0,'color1':gcolor,'color2':gcolor2,'color3':gcolor3,\
                   'nbin':nbinlums,'label':_ginfo}
        optparams = ['distmin','distmax','nbin']
        if pmet==3 and 'galaxy' in _showmap:pstplot.interactive_show(pstplot.lumsview,pparams,optparams)
        elif pmet in [1,2] and 'galaxy' in _showmap:fig_g4 = pstplot.lumsview(pparams)
        if pmet==2 and 'galaxy' in _showmap:input('galaxy4')

    print("%i galaxies generated in %i sec"%(len(ra0),int(time.time()-start_time)))
    return np.arange(len(ra0)),_gname,ra0,dec0,mag0,dist0,fig_g1,fig_g3,fig_g4
