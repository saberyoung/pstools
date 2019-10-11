"""############################################################################ 
2019/1/30 Start
priorzation algorithm: by convolution of galaxt mass map, 
trigger probability, galaxy distance distribution (so far), and more?
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
import pst

#################################################
def make_hpfitsmap(_ra,_dec,_cat,nside):

    galpixels_Range_lum= np.zeros(hp.nside2npix(nside))
    pix_num_Range    = (pst.DeclRaToIndex(_dec,_ra,nside)) 
    galpixels_Range_lum[pix_num_Range]+=_cat    
    lmap = galpixels_Range_lum    
    lmap = lmap/sum(lmap) # normalization
    return lmap

def dist_galaxydist(_dist,tdist,verbose=False):   

    score = []
    if not tdist:return
    dmean,dsigma = tdist.split(',')
    dmean,dsigma = float(dmean),float(dsigma)
    for dist in _dist:
        _score = np.e**(-(dist-dmean)**2/2./dsigma**2)        
        score.append(_score)
    if verbose:pst.dist_gauss(dmean,dsigma)
    return np.array(score)

def calprob_gal(skymap,ra,dec,radius=False):   

    if radius:
        problist,nside = [],hp.get_nside(skymap)
        for ii,xx in enumerate(ra):
            theta,phi = pst.RadecToThetaphi(ra[ii],dec[ii])        
            vec = hp.ang2vec(theta,phi)
            pix = hp.query_disc(nside,vec,radius)
            _theta,_phi = hp.pix2ang(nside,pix)
            _prob = hp.get_interp_val(skymap, _theta, _phi)   
            prob = sum(_prob)                      
            problist.append(prob)                                  
    else:
        theta,phi = pst.RadecToThetaphi(ra,dec)
        problist = hp.get_interp_val(skymap, theta, phi) 
    return problist

def calprob_tile(skymap,ral,decl,fovh,fovw): 

    problist,nside = [],hp.get_nside(skymap)
    for _ra,_dec in zip(ral,decl):
        if _ra>360-fovw:_ra=360-fovw
        if _ra<fovw:_ra=fovw
        if _dec>90-fovh:_dec=90-fovh
        if _dec<-90+fovh:_dec=-90+fovh
        ipix_poly=(pst.ipix_in_box(_ra,_dec,fovh,fovw,nside))
        _probs = skymap[ipix_poly].sum()
        problist.append(_probs)
    return np.array(problist)
