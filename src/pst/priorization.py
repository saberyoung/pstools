"""############################################################################ 
2019/1/30 Start
priorzation algorithm: by convolution of galaxt mass map, trigger probability, galaxy distance distribution (so far), and more?
""" ############################################################################

from __future__ import print_function
from builtins import input
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import pylab as pl
import matplotlib.colors as colors
from matplotlib import mlab
import os,sys,re,shutil,glob,json,subprocess
import math as mt
import smtplib,tempfile
import requests
from astropy.io import fits
import astropy.coordinates
import astropy.time
import astropy.units as u
from pst import pstplot,pstdef,scheme

#################################################
def make_hpfitsmap(_ra,_dec,_cat,_dist,nside,coord,ordering,norm,verbose,\
                   interactive,rot_phi,rot_theta,radius,fig=4,pmet=False):

    galpixels_Range_lum= np.zeros(hp.nside2npix(nside))
    pix_num_Range    = (pstdef.DeclRaToIndex(_dec,_ra,nside)) 
    galpixels_Range_lum[pix_num_Range]+=_cat
    
    if radius:lmap = hp.sphtfunc.smoothing(galpixels_Range_lum,fwhm=radius) # mass map
    else:lmap = galpixels_Range_lum       
                  
    pparams = {'hpmap':lmap,'title':'catalog','rot_phi':rot_phi,\
               'rot_theta':rot_theta,'fignum':fig,'ordering':ordering,\
               'coord':coord,'norm':norm}
    
    optparams = ['rot_theta','rot_phi']
    if pmet==3:pstplot.interactive_show(pstplot.mollview,pparams,optparams)
    elif pmet in [1,2]:fig = pstplot.mollview(pparams)       
    if pmet==2:input('galaxy5')
    return lmap,fig

def dist_galaxydist(_dist,tdist,verbose=False):   

    score = []
    if not tdist:return
    dmean,dsigma = tdist.split(',')
    dmean,dsigma = float(dmean),float(dsigma)
    for dist in _dist:
        _score = np.e**(-(dist-dmean)**2/2./dsigma**2)        
        score.append(_score)
    if verbose:pstplot.dist_gauss(dmean,dsigma)
    return np.array(score)

def calprob_gal(skymap,ra,dec,id0,radius=False):   

    if radius>0:
        problist,nside = [],hp.get_nside(skymap)
        for ii,xx in enumerate(ra):
            theta,phi = pstdef.RadecToThetaphi(ra[ii],dec[ii])        
            vec = hp.ang2vec(theta,phi)
            pix = hp.query_disc(nside,vec,radius)
            _theta,_phi = hp.pix2ang(nside,pix)
            _prob = hp.get_interp_val(skymap, _theta, _phi)   
            prob = sum(_prob)                      
            problist.append(prob)             
        problist=np.array(problist)        

    else:
        theta,phi = pstdef.RadecToThetaphi(ra,dec)
        problist = hp.get_interp_val(skymap, theta, phi)   

    # transform to OB format
    ral,decl,idl=[],[],[]
    for _ra,_dec,_id in zip(ra,dec,id0):
        ral.append([_ra])
        decl.append([_dec])
        idl.append(_id)
    return np.array(idl),np.array(ral),np.array(decl),problist

def calprob_tile(skymap,ra,dec,fovh,fovw,obx,oby): 

    if True:
        # cal via pointings
        problist,nside = [],hp.get_nside(skymap)
        for _ral,_decl in zip(ra,dec):
            _probs=0
            for _ra,_dec in zip(_ral,_decl):
                if _ra>360-fovw:_ra=360-fovw
                if _ra<fovw:_ra=fovw
                if _dec>90-fovh:_dec=90-fovh
                if _dec<-90+fovh:_dec=-90+fovh
                ipix_poly=(pstdef.ipix_in_box(_ra,_dec,fovh,fovw,nside))
                _probs+=skymap[ipix_poly].sum()
            problist.append(_probs)
    else:
        # cal via OBs
        _ra2,_dec2 = [],[]
        for _ra,_dec in zip(ra,dec):
            _ra2.append(np.mean(_ra))
            _dec2.append(np.mean(_dec))

        # 
        problist,nside = [],hp.get_nside(skymap)
        for _ra,_dec in zip(_ra2,_dec2):
            if _ra>360-fovw*obx:_ra=360-fovw*obx
            if _ra<fovw*obx:_ra=fovw*obx
            if _dec>90-fovh*oby:_dec=90-fovh*oby
            if _dec<-90+fovh*oby:_dec=-90+fovh*oby
            ipix_poly=(pstdef.ipix_in_box(_ra,_dec,fovh*oby,fovw*obx,nside))
            problist.append(skymap[ipix_poly].sum())
        
    return np.array(ra),np.array(dec),np.array(problist)
