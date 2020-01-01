"""############################################################################ 
2019/1/30 Start
define some functions in assiting/supporting
""" ############################################################################
from __future__ import print_function
from builtins import input

import time
import numpy as np
import healpy as hp
import pylab as pl
import matplotlib.pyplot as plt
import random
import os,sys,glob,shutil

from astropy.io import fits
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
import scipy.stats
import pst
from pst.default import *

_pstpath = pst.__path__[0]
#################################################

class main(object):
    
    def __init__(self,optlist):

        ''' Main procedure start '''
        start_time = time.time()
        
        # try cleaning all the plots if any
        plt.close('all')

        # how to deal with plots        
        self.showmode = int(optlist['arg']['show']["showmode"])        
        if self.showmode in [3,4,5]: pl.ion()

        # which plots to show
        if len(optlist['arg']['show']["showmap"])>0:
            self.showmap = [int(gg) for gg in \
            optlist['arg']['show']["showmap"].split(',')]
        else: self.showmap = []

        # read params
        self.optlist = optlist
        self.str1 = self.optlist['arg']['email']['emailcontent']
        self.str2 = self.optlist['arg']['phone']['phonecontent']
        self.str3 = self.optlist['arg']['slack']['slackcontent']
        self.verbose = eval(self.optlist['arg']['show']['verbose'])

        # define attachments
        self.distdict, self.indict = {}, {}
        self.filelist, self.figlist = self.optlist['tmp']['files'], []
        self.sname = '%s/$telname$_%s.txt'%(\
                    self.optlist['arg']['data']['dir'],\
                    str(astropy.time.Time.now().value)[:-10])

        # define nside
        if len(self.optlist['tmp']['tmap'])>0:
            self.nside = hp.get_nside(self.optlist['tmp']['tmap'])           
        else: self.nside = int(nside)

        ''' 1 - Decide telescope schedule '''
        if self.verbose:
            print (' [1] Read telescope and strategy: \n')
        self.dec_scheduler()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 2 - Check triggers (hp map) '''
        if self.verbose: 
            print (' [2] Check trigger infos: \n')
        self.triggers(contours)
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        # action = 1 will stop here
        if int(self.optlist['arg']['react']['action'])==1:
            return
        else:
            # go on, refine files
            self.filelist, self.figlist = [], []

        ''' 3 - Obtain galaxies '''
        if self.verbose:
            print (' [3] Obtain galaxies: \n')
        self.galaxies()       
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 4 - Priorization algorithm '''
        if self.verbose: 
            print (' [4] Generate 3D priorization map: \n')
        self.priorization()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 5 - generate pointings 
        for tiling search '''
        if self.verbose: 
            print (' [5] Obtain pointings: \n')
        self.pointings()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 6 - rank pointings/galaxies '''
        if self.verbose: 
            print (' [6] Ranking these fields: \n')
        self.ranking()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 7 - remain only important pointings/galaxies '''
        if self.verbose: 
            print (' [7] Remain interesting fields: \n')
        self.cutfields()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 8 - scheduler pointings/galaxies '''
        if self.verbose:
            print (' [8] arrange fields: \n')
        self.scheduler()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 9 - Visualization '''
        if self.verbose: 
            print (' [9] Visualize above processes: \n')
        self.visualization2()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

        ''' 10 - distributing to users/api/... '''
        if self.verbose: 
            print (' [10] Distribute to users/api: \n')
        self.send()
        if self.verbose: 
            print("\t...Finished in %i secs \n"%(int(time.time()-start_time)))

    def dec_scheduler(self):
        # decide which strategy to use

        self.schlist = {}
        self.schlist['T'] = []
        self.schlist['G'] = []
        for _ff in self.optlist:            
            if _ff in ['arg','tmp']:continue
            if self.optlist[_ff]['pointings']['scheduler'] == 'A':
                if float(self.optlist[_ff]['telescope']['fovw']) > .2 and \
                   float(self.optlist[_ff]['telescope']['fovh']) > .2:
                    # greater than .2*.2 sq deg FoV          
                    self.optlist[_ff]['pointings']['scheduler'] = 'T'
                    self.schlist['T'].append(self.optlist[_ff])
                else:
                    self.optlist[_ff]['pointings']['scheduler'] = 'G'
                    self.schlist['G'].append(self.optlist[_ff])
            elif self.optlist[_ff]['pointings']['scheduler'] in ['G','T']:
                self.schlist[self.optlist[_ff]['pointings']['scheduler']].\
                    append(self.optlist[_ff])
            if self.verbose:
                print ('\t%s: %s search\n'%\
                       (self.optlist[_ff]['telescope']['name'], \
                        self.optlist[_ff]['pointings']['scheduler']))

    def trigger_validation(self):
        # validate a trigger by using the header of fits
        ''' judge interests of the trigger '''

        # read more from the xml, complementry to the fits header                
        for _ii in strxml:
            if _ii in self.optlist['tmp']['voevent']:
                self.str1 += '#\t%s:\t%s\n'%(_ii,self.optlist['tmp']['voevent'][_ii])
                self.str2 += '%s:%s, '%(_ii,self.optlist['tmp']['voevent'][_ii])
                self.str3 += '`%s`:%s, \n'%(_ii,self.optlist['tmp']['voevent'][_ii])
                self.indict['`%s`'%_ii] = self.optlist['tmp']['voevent'][_ii]
                if _ii == 'GraceID':
                    self.graceid = self.optlist['tmp']['voevent'][_ii]
                    # use trigger name and date
                    self.sname = '%s/$telname$_%s_%s.txt'%(\
                        self.optlist['arg']['data']['dir'],self.graceid,\
                        str(astropy.time.Time.now().value)[:-10])
            else:
                self.str1 += '#\t%s:\tNot available\n'%_ii
                self.str2 += '%s:None, '%_ii
                self.str3 += '`%s`:None, \n'%_ii

        # show params of the fits header
        for _ii in strhdr:
            if _ii in self.header:            
                self.str1 += '#\t%s:\t%s\n'%(_ii,self.header[_ii])
                self.str2 += '%s:%s, '%(_ii,self.header[_ii])
                self.str3 += '`%s`:%s, \n'%(_ii,self.header[_ii])
                self.indict['`%s`'%_ii] = self.header[_ii]
            else:
                self.str1 += '#\t%s:\tNone\n'%_ii
                self.str2 += '%s:None, '%_ii
                self.str3 += '`%s`:None, \n'%_ii

        # check the area of skymap     
        # get area in single pixel, unit in deg
        self.areasingle = (hp.nside2resol(\
                            self.nside,arcmin=True)/60.)**2

        self.area={}
        for _cc in self.indexlist: # levels: .9,.5,....
            _num = len(self.indexlist[_cc])          
            if _num>0:
                self.str1 += ' %.2fsky: %.f sq deg\n'%\
                             (_cc,_num*self.areasingle)
                self.str2 += '%.2fsky: %.f sq deg, '%\
                             (_cc,_num*self.areasingle)
                self.str3 += '`%.2fsky`: %.f sq deg, \n'%\
                             (_cc,_num*self.areasingle)
                self.distdict[_cc] = ' %.f%% sky:%.f $deg^2$\n'%\
                            (_cc*100,_num*self.areasingle)
                self.indict['`loc%i`'%(_cc*100)] = _num*self.areasingle
                self.area[_cc] = _num*self.areasingle               
            else:
                self.str1 += ' %.2f%sky:None\n'%(_cc)
                self.str2 += '%.2fsky:None, '%(_cc)
                self.str3 += '`%.2fsky`:None, \n'%(_cc)
                self.distdict[_cc] = ' %.f%% sky:None\n'%(_cc*100)
                self.area[_cc] = None       

        self.str1 += '\n'
        self.str2 = self.str2[:-2]
        self.str3 = self.str3[:-2]
        if self.verbose: print (self.str1)

    def triggers(self,contours):
        # check if need trigger

        if eval(self.optlist['arg']['priorization']['trigger']): # read if yes

            # read and print some values from the FITS header.
            self.header = dict(self.optlist['tmp']['header'])

            # obtain containment contour of GW PDF
            # binnned map to lower resolution in order to save time
            _probs = hp.pixelfunc.ud_grade(self.optlist['tmp']['tmap'],nside1)
            _probs = _probs/np.sum(_probs)

            # computer contour lines and sky area at specific levels
            # add contour for tel
            for ff in self.optlist:
                try:
                    contours.append(float(self.optlist[ff]['observe']['remcov']))
                except: pass
            contours = np.unique(contours)

            # compute contours for plot
            self.theta_contour, self.phi_contour = \
                    pst.compute_contours(contours,_probs)   
            
            # compute contours index
            self.indexlist = pst.compute_contours_1(contours, \
                        self.optlist['tmp']['tmap'])            

            # validate skymap
            self.trigger_validation()

            # immediately report before telescopes scheduling
            self.visualization1()

            # fast plot
            try: 
                self.fig_2d.savefig(self.sname.replace('$telname$_','').replace('.txt','.png'))
                self.figlist.append(self.sname.replace('$telname$_','').replace('.txt','.png'))
            except: pass
            
            # - email
            if eval(self.optlist['arg']['email']['activate']):
                for _toaddress in self.optlist['arg']['email']['to'].split(','):
                    _sent = pst.sendemail(self.optlist['arg']['email']['from'],\
                            self.optlist['arg']['email']['pwd'],\
                            self.optlist['arg']['email']['smtp'],\
                            self.optlist['arg']['email']['sub'],\
                            self.optlist['arg']['email']['from'],\
                            _toaddress,self.str1,self.figlist,self.filelist)
                    if self.verbose:
                        if _sent: print ('email sent successful to %s'%_toaddress)
                        else: print ('email sent failed to %s'%_toaddress)             

            # - slack
            if eval(self.optlist['arg']['slack']['activate']):
                for _toaddress in self.optlist['arg']['slack']['channel'].split(','):
                    _sent = pst.slack(self.optlist['arg']['slack']['token'], \
                                      _toaddress, self.str3, self.figlist+self.filelist)
                    if self.verbose:
                        if _sent: print ('slack sent successful to %s'%_toaddress)
                        else: print ('slack sent failed to %s'%_toaddress) 

            # - SMS
            if eval(self.optlist['arg']['phone']['activate']):
                for _toaddress in self.optlist['arg']['phone']['to'].split(','):
                    _sent = pst.phone(self.optlist['arg']['phone']['account'],\
                                     self.optlist['arg']['phone']['token'],\
                                     self.optlist['arg']['phone']['from'],\
                                     _usr,self.str2)
                    if self.verbose:
                        if _sent: print ('SMS sent successful to %s'%_toaddress)
                        else: print ('SMS sent failed to %s\tno twilio'%_toaddress)

            # - excute python codes
            #   for special functions, e.g. insert to DB               
            if len(self.optlist['arg']['data']['py1'])>0:
                import subprocess
                for _py in self.optlist['arg']['data']['py1'].split(','):
                    pid = subprocess.Popen(['python',_py,str(self.indict)],\
                                    stdout=subprocess.PIPE,\
                                    stderr=subprocess.PIPE)
                    output,error = pid.communicate()
                    if self.verbose: 
                        if len(output)>0: print('### Info: %s'%output)
                        if len(error)>0: print('### Error: %s'%error)
            if self.verbose: print('\n')
        else: self.theta_contour, self.phi_contour = None, None

    def galaxies(self):
        # check if need galaxies

        if int(self.optlist['arg']["priorization"]["galaxy"]) in [1,2] or \
           len(self.schlist['G'])>0:

            self.catalog = int(self.optlist['arg']['galaxies']['catalog'])
            filtro = self.optlist['arg']['galaxies']['filter']
            limrag = [float(self.optlist['arg']['galaxies']['limra'].split(',')[0]),
                     float(self.optlist['arg']['galaxies']['limra'].split(',')[1])]
            limdecg = [float(self.optlist['arg']['galaxies']['limdec'].split(',')[0]),
                      float(self.optlist['arg']['galaxies']['limdec'].split(',')[1])]
            limdist = [float(self.optlist['arg']['galaxies']['limdist'].split(',')[0]),
                       float(self.optlist['arg']['galaxies']['limdist'].split(',')[1])]
            limmag = [float(self.optlist['arg']['galaxies']['limmag'].split(',')[0]),
                      float(self.optlist['arg']['galaxies']['limmag'].split(',')[1])]

            # query galaxies via vizier
            self.gname,self.gra,self.gdec,self.gmag,self.gdist = \
                    pst.galaxies(catalog=self.catalog,size=size,\
                    limra=limrag, limdec=limdecg, limdist=limdist, \
                    limmag=limmag, filtro=filtro, cachemode=cachemode,\
                    cachefile=cachefileg,verbose = self.verbose)

    def pointings(self):
        # generate pointings

        if len(self.schlist['T']) > 0:
            for _nt,_tt in enumerate(self.schlist['T']):
                # read OB mode
                _obx, _oby = int(_tt['observe']["ob"].split(',')[0]),\
                             int(_tt['observe']["ob"].split(',')[1])
                _fovw, _fovh = float(_tt['telescope']['fovw'])*_oby,\
                               float(_tt['telescope']['fovh'])*_obx
                # update fov
                self.schlist['T'][_nt]['telescope']['fovw'] = _fovw
                self.schlist['T'][_nt]['telescope']['fovh'] = _fovh

                _shift = int(_tt['pointings']["shift"])                
                if _shift < 0: _shift = 0
                limrat,limdect = \
                    [float(zz) for zz in _tt['pointings']["limra"].split(',')],\
                    [float(zz) for zz in _tt['pointings']["limdec"].split(',')]
                if _shift == 0: # with shfitw=shifth=0
                    if self.verbose:
                        print ('\tgenerate pointing without shift mode')
                    if '%s' in cachefilet:
                        _cachefile = cachefilet%_tt['telescope']['name']
                    else:
                        _cachefile = cachefilet
                    _ral,_decl = pst.pointings(
                        tel=_tt['telescope']['name'],limdec=limdect,\
                        limra=limrat,fovh=_fovh,fovw=_fovw,\
                        shifth=0.,shiftw=0.,verbose=self.verbose,\
                        uarea=float(_tt['pointings']['uarea']),\
                        cachemode=cachemode,cachefile=_cachefile,\
                        skipfile=_tt['pointings']['skipfile'])
                else: # running shfitw, shifth
                    if self.verbose:
                        print ('\tgenerate pointing with shift mode')
                    _ral,_decl = pst.pointngsshift(
                        self.pmap,_shift,verbose=self.verbose,\
                        limdec=limdect,limra=limrat,fovh=_fovh,\
                        fovw=_fovw,skipfile=_tt['pointings']['skipfile'])
                self.schlist['T'][_nt]['ra'] = _ral
                self.schlist['T'][_nt]['dec'] = _decl                

    def priorization(self):
        # priorization with trigger
        if eval(self.optlist['arg']["priorization"]["trigger"]): # trigger
            if self.verbose:
                print ('\t+ trigger prob')
            self.pmap = self.optlist['tmp']['tmap']
        else:
            self.pmap = np.zeros(hp.nside2npix(self.nside))+1

        # consider dist distrbution
        if eval(self.optlist['arg']["priorization"]["dist"]) and \
           hasattr(self, 'gra') and self.optlist['tmp']['distmu'] is \
           not None and self.optlist['tmp']['distsigma'] is not None \
           and self.optlist['tmp']['distnorm'] is not None:
            # if 3D GW loc available
            # query the min to max distance from catalog
            # meanwhile, get dist range for every piece of sky 
            # (use small nside, since it would not change too much)
            # the dist range of each piece will be then included,
            # together with selected galaxies (gaussian), calculated
            # for a score, contributing to the final list                
            if self.verbose:
                print ('\t+ distance distribution')

            self.gdscore = pst.gwdist(self.nside,\
                    self.optlist['tmp']['distmu'], \
                    self.optlist['tmp']['distsigma'], \
                    self.optlist['tmp']['distnorm'], \
                    self.gra,self.gdec,self.gdist)

        # convolve galaxy info into score
        if int(self.optlist['arg']["priorization"]["galaxy"]) in [1,2]:            
            if int(self.optlist['arg']["priorization"]["galaxy"])==1: # lums
                _cat = 10**((-1)*(self.gmag/2.5))
                if self.verbose: print ('\t+ galaxy lums')
            elif int(self.optlist['arg']["priorization"]["galaxy"])==2: # counts
                _cat = np.zeros(len(self.gra))+1
                if self.verbose: print ('\t+ galaxy counts')
            try:_cat *= self.gdscore
            except:pass

            self.gmap = pst.make_hpfitsmap(self.gra,self.gdec,_cat,self.nside)

            # cross map with trigger map
            self.pmap *= self.gmap

        # normalize
        self.pmap = self.pmap/sum(self.pmap)
        if self.verbose: print ('\n')

    def cutfields(self):
        # remove fields, accoroding to the 2D contour and the distance constrains
        # 1- galaxies
        if len(self.schlist['G']) > 0:            
            for _nt,_tt in enumerate(self.schlist['G']):
                try:_cov = float(_tt['observe']['remcov'])
                except:_cov = False
                if _cov<=0 and _cov>=1:_cov=False

                # cut galaxies on score
                if self.verbose:
                    print ('\tG: fields threshold set to %s'%remfields)
                for _key,_val in zip(['ra','dec','score',\
                        'name','dist','mag'],[self.gra,self.gdec,\
                        self.gscore,self.gname,self.gdist,self.gmag]):
                    self.schlist['G'][_nt][_key] = _val[:remfields]
                if self.verbose: 
                    print ('\tremian %i galaxies for %s'%\
                           (len(self.schlist['G'][_nt]['ra']),\
                            _tt['telescope']['name']))

                # cut galaxies on 2d sky
                if self.verbose:
                    print ('\tG: 2d region set to %s'%_cov)
                if hasattr(self, 'indexlist') and _cov:
                    gtheta,gphi = pst.RadecToThetaphi(self.schlist['G'][_nt]['ra'],\
                                                      self.schlist['G'][_nt]['dec'])
                    _aindex = hp.ang2pix(self.nside, gtheta, gphi, \
                            nest=ordering,lonlat=False)
                    _cindex = list(set(_aindex)&set(self.indexlist[_cov]))
                    if len(_cindex)>0: #covered
                        _index = np.argwhere([rr in _cindex for rr in _aindex])
                    else: _index = []

                    for _key in ['ra','dec','score','name','dist','mag']:
                        if len(_index) > 0:
                            self.schlist['G'][_nt][_key] = \
                            self.schlist['G'][_nt][_key][_index]
                        else: self.schlist['G'][_nt][_key] = []

                if self.verbose: 
                    print ('\tremian %i galaxies for %s\n'%\
                           (len(self.schlist['G'][_nt]['ra']),\
                            _tt['telescope']['name']))

        # 2- tilings
        if len(self.schlist['T']) > 0:
            for _nt,_tel in enumerate(self.schlist['T']):
                try:_cov = float(_tel['observe']['remcov'])
                except:_cov = False
                if _cov<=0 and _cov>=1:_cov=False

                # cut pointings on score
                if self.verbose: 
                    print ('\tT: fields threshold set to %s'%remfields)
                for _key in ['ra','dec','score']:
                    self.schlist['T'][_nt][_key] = \
                            self.schlist['T'][_nt][_key][:remfields]
                if self.verbose: 
                    print ('\tremian %i pointings for %s'%\
                           (len(self.schlist['T'][_nt]['ra']),\
                            _tel['telescope']['name']))

                # cut pointings on 2d sky
                if self.verbose:
                    print ('\tT: 2d region set to %s'%_cov)
                if hasattr(self, 'indexlist') and _cov:
                    _xx,_yy = float(_tel['telescope']['fovw']),\
                              float(_tel['telescope']['fovh'])            
                    _indexl = []
                    for _index, (_pra,_pdec) in \
                        enumerate(zip(self.schlist['T'][_nt]['ra'], \
                                      self.schlist['T'][_nt]['dec'])):
                        _aindex = pst.ipix_in_box(_pra,_pdec,_xx,_yy,self.nside)
                        _cindex = list(set(_aindex)&set(self.indexlist[_cov]))
                        if len(_cindex)>0: _indexl.append(_index)
                
                    for _key in ['ra','dec','score']:
                        if len(_indexl) > 0:
                            self.schlist['T'][_nt][_key] = \
                                    self.schlist['T'][_nt][_key][np.array(_indexl)]
                        else:self.schlist['T'][_nt][_key] = []
                if self.verbose: 
                    print ('\tremian %i pointings for %s\n'%\
                        (len(self.schlist['T'][_nt]['ra']),_tel['telescope']['name']))

    def ranking(self):
        # ranking for tiling or galaxies

        if len(self.schlist['G'])>0: # galaxy search
            self.gscore = pst.calprob_gal(self.pmap,self.gra,self.gdec)            

            # sort galaxies
            idx=np.argsort(np.asarray(self.gscore))[::-1]
            self.gscore,self.gname,self.gra,self.gdec,self.gmag,self.gdist=\
                    self.gscore[idx],self.gname[idx],self.gra[idx],\
                    self.gdec[idx],self.gmag[idx],self.gdist[idx]
            if self.verbose: print ('\tranking galaxies done')

        if len(self.schlist['T'])>0: # tiliing search
            for _nt,_tel in enumerate(self.schlist['T']):
                _tra,_tdec,_fovw,_fovh = _tel['ra'],_tel['dec'],\
                    float(_tel['telescope']['fovw']),\
                    float(_tel['telescope']['fovh'])
                _tscore = pst.calprob_tile(self.pmap,_tra,_tdec,_fovw,_fovh)

                # sort pointings
                idx=np.argsort(np.asarray(_tscore))[::-1]
                _tscore,_tra,_tdec=_tscore[idx],_tra[idx],_tdec[idx]
                self.schlist['T'][_nt]['score'] = _tscore
                self.schlist['T'][_nt]['ra'] = _tra
                self.schlist['T'][_nt]['dec'] = _tdec
                if self.verbose: 
                    print ('\tranking %s fields done'%_tel['telescope']['name'])
        if self.verbose:print('\n')

    def scheduler(self):

        # pointing arrangement      
        self.tellist = {}
        self.tellist['G'] = {}
        self.tellist['T'] = {}

        # 1- galaxies
        if len(self.schlist['G']) > 0:  
            wlist = {}
            for _nt,_tt in enumerate(self.schlist['G']):
                wlist[_nt] = float(_tt['telescope']['weight']) 

            _cralist, _cdeclist = [], []
            for weight in sorted(np.unique(list(wlist.values()))): 
                if self.verbose:
                    print ('\n\t ### schedule telescope with weight=%s'%weight) 
                tellist = self.prob_obs(np.array(self.schlist['G'])[\
                        np.where(list(wlist.values())==weight)],\
                        _cralist, _cdeclist, False, False)
                self.tellist['G'][weight] = tellist
                for _nt,_tt in enumerate(tellist):
                    for _ra,_dec in zip(tellist[_nt]['ra'],\
                                        tellist[_nt]['dec']):
                        _cralist.append(_ra)
                        _cdeclist.append(_dec)

        # 2- tilings
        if len(self.schlist['T']) > 0:
            wlist = {}
            for _nt,_tt in enumerate(self.schlist['T']):
                wlist[_nt] = float(_tt['telescope']['weight']) 
            _cralist, _cdeclist, _cfovw, _cfovh = [], [], [], []
            for weight in sorted(np.unique(list(wlist.values()))):   
                if self.verbose:
                    print ('\n\t ### schedule telescope with weight=%s'%weight)
                tellist = self.prob_obs(np.array(self.schlist['T'])[\
                        np.where(list(wlist.values())==weight)],\
                        _cralist, _cdeclist, _cfovw, _cfovh)
                self.tellist['T'][weight] = tellist
                for _nt,_tt in enumerate(tellist):
                    for _ra,_dec,_fovw,_fovh in zip(tellist[_nt]['ra'],\
                                                    tellist[_nt]['dec'],\
                                                    tellist[_nt]['fovw'],\
                                                    tellist[_nt]['fovh']):
                        _cralist.append(_ra)
                        _cdeclist.append(_dec)
                        _cfovw.append(_fovw)
                        _cfovh.append(_fovh)

        # 3- generate scheduler files
        # write files
        self.schfile = {}
        for _sch in self.schlist: # schedule: T/G
            if self.verbose:
                print ('\n\t ### output schedule for %s strategy'%_sch)

            for _nt,_tt in enumerate(self.schlist[_sch]): # telescope
                if self.verbose:
                    print ('\n\t\tfor tel:%s'%_tt['telescope']['name'])
                fovw,fovh = float(_tt['telescope']['fovw']),\
                            float(_tt['telescope']['fovh'])
                _obx,_oby = int(_tt['observe']["ob"].split(',')[0]),\
                            int(_tt['observe']["ob"].split(',')[1])
                _nob = 0
                # get selected ra,dec to tellist
                tellist = None
                for _weight in self.tellist[_sch]:
                    for _ntt in self.tellist[_sch][_weight]:
                        if self.tellist[_sch][_weight][_ntt]['name'] ==\
                           _tt['telescope']['name']:
                            tellist = self.tellist[_sch][_weight][_ntt]

                if not tellist is None:
                    if len(_tt['scheduler']['schfile'])>0:
                        if _tt['scheduler']['schfile'] == 'auto':
                            self.schfile[_tt['telescope']['name']] = \
                                        self.sname.replace('$telname$',\
                                        _tt['telescope']['name'])
                        else:
                            self.schfile[_tt['telescope']['name']] = \
                                        _tt['scheduler']['schfile']
                    else:
                        self.schfile[_tt['telescope']['name']] = None

                    if self.schfile[_tt['telescope']['name']] is not None:
                        # read ra,dec,time
                        _strlist,sumscore = '',0
                        _hdr,_whdr = '#RA(J2000)   DEC(J2000)   Filter      Nexp     Prob(%)'+\
                                     ' '*6,True
                        ra,dec,date = tellist['ra'],\
                                      tellist['dec'],\
                                      tellist['timelist']
                        _obs, _tf, _tfl = self.def_obs(_tt)
                        [_exp,_rot,_nf,_nexp] = _tfl

                        # if dither or filter
                        _dit = float(_tt['scheduler']['dither'])/3600
                        _bands = _tt['scheduler']['filter'].split(',')

                        for _zz in range(len(date)):
                            # at different time
                            _time = date[_zz]
                            _num = 0

                            for _nfi in range(_nf):
                                # number of fields
                                try:
                                    _rac,_decc = ra[_zz*_nf+_nfi],\
                                                 dec[_zz*_nf+_nfi]
                                    _nob += 1
                                except: 
                                    continue
                                s,t = np.where(_tt['ra']==_rac),\
                                      np.where(_tt['dec']==_decc)

                                # find index at specific ra,dec
                                index = np.intersect1d(s, t)[0]
                                score = pst.decomposit(_tt['score'])[index]
                                sumscore+=score
                                
                                # OB mode
                                ral, decl = pst.divide_OB(_rac,_decc,fovw,fovh,_obx,_oby)
                                for _pointing,(_ra,_dec) in enumerate(zip(ral, decl)):
                                    # for 1 pointing of OB
                                    for _nfilt in range(len(_bands)):
                                        # for different filter
                                        _filt = _bands[_nfilt]

                                        # write ra,dec,score,[name,dist,mag],....
                                        _strlist += '%-*.5f %-*.5f %-*s %-*s %-*.2e '%\
                                                    (12,_ra,12,_dec,12,_filt,12,_nexp,12,score)
                                        if _sch == 'G':
                                            _name = pst.decomposit(_tt['name'])[index]
                                            _strlist += '%-*s '%(30,_name)
                                            if _whdr:_hdr += (' '*12+'name'+' '*16)
                                            for _key in ['dist','mag']:
                                                _val = pst.decomposit(_tt[_key])[index]
                                                _strlist += '%-*.2f '%(8,_val)
                                                if _whdr:_hdr += '%s      '%_key
                                        # get time
                                        _timenow = _time + _exp*_num
                                        _num += _nexp

                                        # calculate sun and airmass
                                        radecs = astropy.coordinates.SkyCoord(ra=_ra*u.deg, \
                                                                              dec=_dec*u.deg)
                                        frame = astropy.coordinates.AltAz(obstime=_timenow, \
                                                                          location=_obs)     
                                        altaz = radecs.transform_to(frame)
                                        sun_altaz = astropy.coordinates.\
                                                    get_sun(_timenow).transform_to(altaz)
                                        # write time,sun,airmass
                                        _strlist += '%-*s %-*.2f %-*s %i-%i \n'%\
                                                    (6,str(sun_altaz.alt)[:5],8,\
                                                     altaz.secz,20,str(_timenow)[:19],\
                                                     _nob,_pointing+1)
                                        if _whdr:
                                            _hdr += 'sun  airmass        time'
                                            _hdr += ' '*10
                                            _hdr += 'OB \n'
                                            _whdr = False
                        # write
                        if os.path.exists(self.schfile[_tt['telescope']['name']]):
                            os.remove(self.schfile[_tt['telescope']['name']])
                        rr = open(self.schfile[_tt['telescope']['name']],'w')
                        rr.write('# %i fields covered %.2f%% Prob\n\n'%\
                                 (len(ra),sumscore*100))
                        rr.write(_hdr)
                        rr.write(_strlist)
                        rr.close()
                        if self.verbose:
                            print ('\t%s generated'%self.schfile[_tt['telescope']['name']])
                        self.filelist.append(self.schfile[_tt['telescope']['name']])

    def def_time(self,_tlist): # define time
        _obstime = _tlist['observe']['obstime']
        if _obstime == 'now': 
            timenow = astropy.time.Time.now()  
        elif _obstime == 'night':
            # define next sunset/sunrise time
            _timenow = astropy.time.Time.now()  
            _obs, _tf, _tfl = self.def_obs(_tlist)
            timenow = pst.sunset(_timenow,_obs,\
                float(_tlist['observe']['limsun']))
            _timedif = timenow - _timenow
            if self.verbose:
                print ('\t%.2f hours later sunset for %s'%\
                (_timedif.value*24,_tlist['telescope']['name']))
        elif type(_obstime) in [int,float]:
            timenow = astropy.time.Time.now() + \
                      astropy.time.TimeDelta(float(_obstime)*60, \
                                             format='sec')
        else:
            try: 
                timenow = astropy.time.Time(_obstime, scale='utc')
            except: 
                timenow = astropy.time.Time.now()  
        jdnow = timenow.jd
        return timenow, jdnow

    def def_obs(self,_tlist): # define observatory and time per frame

        # Geodetic coordinates of observatory    
        observatory = astropy.coordinates.EarthLocation(\
                    lat=float(_tlist['telescope']['lat'])*u.deg, \
                    lon=float(_tlist['telescope']['lon'])*u.deg, \
                    height=float(_tlist['telescope']['alt'])*u.m)

        _nexp = int(_tlist['scheduler']['nexp']) + 1
        _obx,_oby = int(_tlist['observe']["ob"].split(',')[0]),\
                    int(_tlist['observe']["ob"].split(',')[1])
        _tf = float(_tlist['telescope']['exptime'])*_nexp*nfields*_obx*_oby +\
              float(_tlist['telescope']['ovhtime'])
        _tf = astropy.time.TimeDelta(_tf, format='sec') 
        return observatory,_tf,\
            [astropy.time.TimeDelta(float(_tlist['telescope']['exptime']),\
            format='sec'),\
            astropy.time.TimeDelta(float(_tlist['telescope']['ovhtime']),\
            format='sec'),\
             nfields,_nexp]

    def prob_obs(self,_glist,_cralist,_cdeclist,_cfovw,_cfovh):

        start_time = time.time()

        _tellist = {}
        for _nt,_tt in enumerate(_glist):

            _tellist[_nt] = {}

            # define observatory and time per frame
            _obs, _tf, [_exp,_rot,_numfield,_numexp] = \
                                self.def_obs(_tt)

            # define time
            _timenow, _jd = self.def_time(_tt)

            # time last
            _deltat = float(_tt['observe']['timelast'])            
            _deltat = astropy.time.TimeDelta(_deltat*3600., format='sec') 

            # time list
            _timel = np.arange(_timenow, _timenow+_deltat, _tf)

            # store them
            # innight remove time, when the sun is high
            _tellist[_nt]['timelist'] = \
                    pst.innight(_timel,_obs,\
                    float(_tt['observe']['limsun']))
            _tellist[_nt]['tf'] = _tf
            _tellist[_nt]['name'] = _tt['telescope']['name']
            _tellist[_nt]['ra'] = []
            _tellist[_nt]['dec'] = []
            if _glist[_nt]['pointings']['scheduler'] == 'T':
                _tellist[_nt]['fovw'] = []
                _tellist[_nt]['fovh'] = []
            _tellist[_nt]['obs'] = _obs

        # define time list
        _timelist = []
        for zz in _tellist:
            _timelist.extend(_tellist[zz]['timelist'])

        _timelist = np.unique(_timelist)
        _timelist = np.sort(_timelist)

        # For slew angle case, start from latest fields
        ralatest, declatest = {}, {}
        for ntel,tel in enumerate(_glist):
            ralatest[ntel], declatest[ntel] = None, None

        # field assignment start
        for _nt,_time in enumerate(_timelist):

            # stop till: 
            # - 1. the end of time;
            # - 2. used up all available pointings

            _timedif = _time - _timelist[0]
            if self.verbose:
                print ('\n\t+ %.2f hs'%(_timedif.value*24))

            for ntel,tel in enumerate(_glist):

                # check telescope available or not
                if _time in _tellist[ntel]['timelist']:
                    if self.verbose:
                        print ('\n\t - tel: %s avalable now'%\
                               tel['telescope']['name'])
                else: # not available
                    if self.verbose:
                        print ('\n\t - tel: %s not avalable now'%\
                               tel['telescope']['name'])
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # if reached pointing limit
                try:
                    nf = float(_glist[ntel]['observe']['limfields'])
                    if sum([bb is not None for bb in \
                            _tellist[ntel]['ra']]) >= nf:
                        if self.verbose:
                            print ('\t - tel: %s reached pointing limit'%\
                                   tel['telescope']['name'])
                        index = np.where(_tellist[ntel]['timelist']==_time)
                        _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                        continue
                except: 
                    pass

                # if used up all pointings
                if sum([bb is not None for bb in \
                    _tellist[ntel]['ra']]) >= \
                    len(_glist[ntel]['ra']):
                    if self.verbose:
                        print ('\t - tel: %s finished all pointings'%\
                               tel['telescope']['name'])
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # remove duplicated fields
                _tralist,_tdeclist = [],[]
                for _ra,_dec in zip(tel['ra'],tel['dec']):
                    
                    _use = True
                    if tel['pointings']['scheduler'] == 'T':
                        # remain pointings for followup
                        # index different for different telescopes
                        # check if 2 pointings have overlap

                        # current pointing
                        fovw, fovh = float(tel['telescope']['fovw']),\
                                     float(tel['telescope']['fovh'])

                        # remove input ra,dec which should be high priority
                        _idx = pst.overlapregioncut(_ra,_dec,_cralist,_cdeclist)
                        _cov = pst.overlapregion(_ra,_dec,fovw,fovh,\
                                                 np.array(_cralist)[_idx],\
                                                 np.array(_cdeclist)[_idx],\
                                                 np.array(_cfovw)[_idx],\
                                                 np.array(_cfovh)[_idx])
                        if not _cov is None:
                            _frac = _cov/fovw/fovh
                            if _frac > float(tel['pointings']['uarea']):
                                _use = False

                        # remove already done for the same priority
                        # different telescopes
                        for _ntel in _tellist:
                            _idx = pst.overlapregioncut(_ra,_dec,\
                                    _tellist[_ntel]['ra'],_tellist[_ntel]['dec'])
                            _cov = pst.overlapregion(_ra,_dec,fovw,fovh,\
                                    np.array(_tellist[_ntel]['ra'])[_idx],\
                                    np.array(_tellist[_ntel]['dec'])[_idx],\
                                    np.array(_tellist[_ntel]['fovw'])[_idx],\
                                    np.array(_tellist[_ntel]['fovh'])[_idx])
                            if not _cov is None:
                                _frac = _cov/fovw/fovh
                                if _frac > float(tel['pointings']['uarea']):
                                    _use = False

                    elif _glist[ntel]['pointings']['scheduler'] == 'G':
                        # remain galaxies for followup
                        # index same for different telescopes
                        # check if 2 galaxies are same

                        # remove input ra,dec which should be high priority
                        for _ra0,_dec0 in zip(_cralist,_cdeclist):
                            if _ra == _ra0 and _dec == _dec0:
                                _use = False

                        # remove already done for the same priority
                        # different telescopes
                        for _ntel in _tellist:
                            for _ra0,_dec0 in zip(_tellist[_ntel]['ra'],\
                                                  _tellist[_ntel]['dec']):
                                if _ra == _ra0 and _dec == _dec0:
                                    _use = False

                    if _use:
                        _tralist.append(_ra)
                        _tdeclist.append(_dec)

                if self.verbose:
                    print ('\t - %s fields remained for tel %s, after removing duplication'%\
                           (len(_tralist),tel['telescope']['name']))
                if len(_tralist)==0:
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # ra dec
                radecs = astropy.coordinates.SkyCoord(ra=_tralist*u.deg, \
                                                      dec=_tdeclist*u.deg)

                # read observatory
                observatory = _tellist[ntel]['obs']

                # Alt/az reference frame at observatory, now
                frame = astropy.coordinates.AltAz(obstime=_time, location=observatory)

                # Transform grid to alt/az coordinates at observatory, now
                altaz = radecs.transform_to(frame)

                # Where is the sun, now?
                sun_altaz = astropy.coordinates.get_sun(_time).transform_to(altaz)

                # sun, airmass constrains                      
                _limalt,_limsun = float(tel['observe']['limalt']), \
                                  float(tel['observe']['limsun'])               

                _cond = np.logical_and(np.logical_and(altaz.secz <= _limalt, \
                                altaz.secz >= 1),sun_altaz.alt <= _limsun*u.deg)
                gradecs = radecs[_cond]

                # go on or not
                if self.verbose:
                    print ('\t - %s fields visible for tel %s, after sun+airmass constrain'%\
                           (len(gradecs),tel['telescope']['name']))
                if len(gradecs) == 0: 
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # moon constrains
                # define seperation to the moon
                _limmoon = tel['observe']['limmoon']
                if _limmoon == 'auto':
                    if 'T' in str(_time):
                        _year, _month, _day = str(_time.value).split('T')[0].split('-')
                    else:
                        _year, _month, _day = str(_time.value).split()[0].split('-')
                    _date, _status, _light = moon_phase(int(_month),int(_day),int(_year))

                    ''' limitation for the moon / TBD '''
                    if _light>=80: _limmoon = 30
                    elif _light<80 and _light>=40: _limmoon = 20
                    elif _light<40 and _light>=10: _limmoon = 10
                    else:_limmoon = 1
                else:
                    try: _limmoon = float(_limmoon)
                    except: 
                        print ('!!! Warning: moon limitation wrong... set default')
                        _limmoon = 30

                _limmoon *= 3600 # to deg
                _loc = astropy.coordinates.get_body('moon', \
                        _time, observatory).transform_to(altaz)
                sep = gradecs.separation(_loc)              
                _length = len(gradecs[np.where(sep.arcsecond<_limmoon)])
                if _length>0:
                    if self.verbose:
                        print('\t - remove %s sources due to moon'%(_length))                   
                    gradecs = gradecs[sep.arcsecond>_limmoon]

                # go on or not
                if self.verbose:
                    print ('\t - %s fields visible for tel %s, after moon constrain'%\
                           (len(gradecs),tel['telescope']['name']))
                if len(gradecs) == 0: 
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # Where is the other sources inside solor system
                _solorobj = tel['observe']['limsobj'].split(',')  
                try:_limsolor = float(tel['observe']['limsobjr'])*3600 # to deg
                except:_limsolor=False
                if _limsolor:
                    for _source in _solorobj:
                        # astropy.coordinates.solar_system_ephemeris.bodies
                        if _source in ['mercury','venus','mars',\
                                       'jupiter','saturn','uranus','neptune']:
                            _loc = astropy.coordinates.get_body(_source, \
                                    _time, observatory).transform_to(altaz)
                            sep = gradecs.separation(_loc)
                            _length = len(gradecs[np.where(sep.arcsecond<_limsolor)])
                            if _length>0:
                                if self.verbose:
                                    print('\t - remove %s sources due to %s'%(_length,_source))
                                gradecs = gradecs[sep.arcsecond>_limsolor]
                        else:
                            if self.verbose:
                                print ('\t!!! Warning: unknown source %s'%_source)

                # go on or not
                if self.verbose:
                    print ('\t - %s fields visible for tel %s, after solor obj constrain'%\
                           (len(gradecs),tel['telescope']['name']))
                if len(gradecs) == 0: 
                    index = np.where(_tellist[ntel]['timelist']==_time)
                    _tellist[ntel]['timelist'] = \
                            np.delete(_tellist[ntel]['timelist'],\
                            index)
                    continue

                # ranking fields
                _order = int(tel['observe']['order'])

                # rank with prob
                if _order == 1: pass
                elif _order == 2: # from west to east
                    _idrank = np.argsort(gradecs.ra)
                    gradecs = gradecs[_idrank]

                elif _order in [3,4]: # consider slewing angle
                    if ralatest[ntel] is None and declatest[ntel] is None:
                        # initialize start field
                        if _order==3: # start from west
                            _ids = np.argmin(gradecs.ra)
                        elif _order==4:
                            # start from high ranking
                            # Note: the input radecs are already ranked with ranking
                            _ids = 0
                        ras,decs = gradecs.ra[_ids].deg,\
                                   gradecs.dec[_ids].deg
                    else:
                        ras,decs = ralatest[ntel], declatest[ntel]
                    index = pst.slew_angle(ras,decs,gradecs.ra.deg,\
                                           gradecs.dec.deg)
                    gradecs = gradecs[index]
               
                if self.verbose:
                    print ('\t - %s fields ranked for tel %s, with order=%i'%\
                           (len(gradecs),tel['telescope']['name'],_order))

                # append fields
                for _nn,(_ra,_dec) in enumerate(zip(gradecs.ra.deg, gradecs.dec.deg)):
                    if _nn<_numfield: # append numfield fields
                        _tellist[ntel]['ra'].append(_ra)
                        _tellist[ntel]['dec'].append(_dec)
                        ralatest[ntel], declatest[ntel] = _ra, _dec
                if self.verbose:
                    print ('\t%i fields for %s\t'%\
                           (len(_tellist[ntel]['ra']),\
                            tel['telescope']['name']))                
                if _glist[ntel]['pointings']['scheduler'] == 'T':
                    for ii in range(_numfield):
                        _tellist[ntel]['fovw'].append(float(_glist[ntel]['telescope']['fovw']))
                        _tellist[ntel]['fovh'].append(float(_glist[ntel]['telescope']['fovh']))
        if self.verbose:print ('\t### time finished')
        return _tellist

    def visualization1(self):

        ########## 1 - healpix show of trigger
        ''' fignum 1: 2d sky'''
        if 1 in self.showmap or 2 in self.showmap or \
           3 in self.showmap or 7 in self.showmap or \
           8 in self.showmap or 9 in self.showmap:

            # if show trigger, otherwise, show sky
            if 1 in self.showmap and len(self.optlist['tmp']['tmap'])>0: 
                map2show = self.optlist['tmp']['tmap']
                if self.verbose: print (' - show trigger sky')
            else: 
                map2show = np.zeros(12*self.nside**2)
                if self.verbose: print (' - show sky without trigger')

            # parameters for plotting
            pparams = {'hpmap':map2show,'theta':theta,'phi':phi,\
                'fignum':1,'ordering':ordering,'coord':coord,\
                'min':minv,'max':maxv,'title':title,\
                'theta_contour':self.theta_contour,'norm':norm,\
                'phi_contour':self.phi_contour,'timenow':'now',\
                'colors':color_contour,'distinfo':self.distdict,\
                'tellist':self.schlist,'figsize':figsize}
            # parameters that can be changed via interactive mode
            optparams = ['theta','phi','min','max','timenow',\
                         'title','coord']

            _pm = False
            if self.showmode == 5: 
                _pm, self.fig_2d = \
                    pst.interactive_show(pst.mollview,pparams,optparams)
            else:
                self.fig_2d = pst.mollview(pparams)
            if self.showmode == 4:
                input('pause to show sky map')
            if _pm: # update show options
                for _cc in optparams:self.optlist['arg']['show'][_cc]=_pm[_cc]

    def visualization2(self):

        self.ncolor = 0
        cnl = {1:'GLADE',2:'GWGC'}
        ########## 2 - show for galaxies
        if hasattr(self, 'gra'):
            ''' fignum 1: 2d sky
                fignum 2: dist
                fignum 3: lums
            ''' 
            if 2 in self.showmap: # 1 - (all) galaxy 2D distribution
                # focused, no interactive mode
                if len*self(gra) > 10000:
                    print ('!!! Warning: take care, too many galaxies')
                pparams = {'ra':self.gra,'dec':self.gdec,\
                    'theta':theta,'phi':phi,'coord':coord,\
                    'fignum':1,'color':color_field[self.ncolor],\
                    'label':'%s galaxies(%i)'%(cnl[self.catolog],len(self.gra))}
                self.fig_2d = pst.pointview(pparams)
                self.ncolor+=1
                if self.showmode in [4, 5]:
                    input('show all %i galaxies'%len*self(gra))
                

            if 3 in self.showmap: # 1 - (selected) galaxy 2D distribution
                for weight in self.tellist['G']:
                    tellist = self.tellist['G'][weight]
                    for _nt in tellist:
                        _tt = tellist[_nt]
                        pparams = {'ra':_tt['ra'],'dec':_tt['dec'],\
                            'theta':theta,'phi':phi,'coord':coord,\
                            'fignum':1,'color':color_field[self.ncolor],\
                            'label':'%s galaxies(%i)'%(_tt['name'],len(_tt['ra']))}
                        self.fig_2d = pst.pointview(pparams)
                        self.ncolor+=1
                if self.showmode in [4, 5]:
                    input('show selected galaxies in skymap')

            if 4 in self.showmap: # 2 - galaxy distance distribution
                dmin,dmax=min(self.gdist),max(self.gdist)
                pparams = {'distmin':dmin,'distmax':dmax,'dist':self.gdist,\
                        'fignum':2,'figsize':figsize,'color1':'k',\
                        'color2':'r','scale':'linear','nbin':10,\
                        'label':'%s galaxies(%i)'%(cnl[self.catolog],len(self.gra))}
                optparams = ['distmin','distmax','nbin','color1',\
                             'color2','scale']
                if self.showmode == 5:
                    _pm, self.fig_gd = pst.interactive_show(pst.distview,pparams,optparams)
                else:
                    self.fig_gd = pst.distview(pparams) 
                if self.showmode == 4:
                    input('galaxy distance distribution')

            if 5 in self.showmap: # 3 - galaxy liminosity distribution
                dmin,dmax=min(self.gdist),max(self.gdist)
                pparams = {'distmin':dmin,'distmax':dmax,'mag':self.gmag,\
                        'figsize':figsize,'dist':self.gdist,'fignum':3,\
                        'scale':'linear','color1':'r','color2':'grey','nbin':1,\
                        'label':'%s galaxies(%i)'%(cnl[self.catolog],len(self.gra))}
                optparams = ['distmin','distmax','nbin','color1',\
                             'color2','scale']
                if self.showmode == 5:
                    _pm, self.fig_gl = \
                        pst.interactive_show(pst.lumsview,pparams,optparams)
                else:
                    self.fig_gl = pst.lumsview(pparams)
                if self.showmode == 4:
                    input('galaxy luminosity distribution')

            if 6 in self.showmap and len(self.schlist['G']) > 0: 
                # 5 - galaxy cumulative score distribution
                _numgal = []
                for _nt,_tt in enumerate(self.schlist['G']):
                    _numgal.append(len(_tt['ra']))
                pparams = {'full':self.schlist['G'],'number':max(_numgal),\
                           'figsize':figsize,'showname':'False','nameloc':.2,\
                           'select':self.tellist['G'],'fignum':4}
                optparams = ['number','showname','nameloc']
                if self.showmode == 5:
                    _pm, self.fig_cum = \
                        pst.interactive_show(pst.cumshow,pparams,optparams)
                else:
                    self.fig_cum = pst.cumshow(pparams)
                if self.showmode == 4:
                    input('cumulative galaxy score distribution')

        ########## 3 - show for pointings
        if len(self.schlist['T']) > 0:
            if 7 in self.showmap: # 1 - all tilings
                for _nt,_tt in enumerate(self.schlist['T']):
                    pparams = {'ra':_tt['ra'],'dec':_tt['dec'],\
                               'fignum':1,'color':color_field[self.ncolor],\
                               'theta':theta,'phi':phi,'coord':coord,\
                               'fovw':_tt['telescope']['fovw'],\
                               'fovh':_tt['telescope']['fovh'],\
                               'label':'%s all (%i) tilings'%\
                               (_tt['telescope']['name'],len(_tt['ra']))}
                    self.fig_2d = pst.verticeview(pparams)
                    self.ncolor+=1
                if self.showmode in [4, 5]:
                    input('show all tilings in skymap')

            if 8 in self.showmap or 9 in self.showmap: # 2 - tilings/route
                for weight in self.tellist['T']:
                    tellist = self.tellist['T'][weight]
                    for _nt in tellist:
                        _tt = tellist[_nt]
                        if 8 in self.showmap: 
                            pparams = {'ra':_tt['ra'],'dec':_tt['dec'],\
                               'fignum':1,'color':color_field[self.ncolor],\
                               'theta':theta,'phi':phi,'coord':coord,\
                               'fovw':_tt['fovw'],'fovh':_tt['fovh'],\
                                'label':'%s selected (%i) tilings'%\
                                       (_tt['name'],len(_tt['ra']))}
                            self.fig_2d = pst.verticeview(pparams)
                            self.ncolor+=1
                        if 9 in self.showmap:
                            pparams = {'ra':_tt['ra'],'dec':_tt['dec'],\
                                'time':_tt['timelist'],'fignum':1,\
                                'color':color_field[self.ncolor],\
                                'theta':theta,'phi':phi,'coord':coord,\
                                'label':'%s routine'%_tt['name']}
                            self.fig_2d = pst.routeview(pparams)
                            self.ncolor+=1
                if self.showmode in [4, 5]:
                    input('show selected tilings in skymap')

            if 10 in self.showmap and len(self.schlist['T']) > 0: 
                # 3 - tiling cumulative score distribution
                _numtil = []
                for _nt,_tt in enumerate(self.schlist['T']):
                    _numtil.append(len(_tt['ra']))
                pparams = {'full':self.schlist['T'],'number':max(_numtil),\
                           'figsize':figsize,'showname':'False','nameloc':.2,\
                           'select':self.tellist['T'],'fignum':4}
                optparams = ['number','showname','nameloc']
                if self.showmode == 5:
                    _pm, self.fig_cum = \
                        pst.interactive_show(pst.cumshow,pparams,optparams)
                else:
                    self.fig_cum = pst.cumshow(pparams)
                if self.showmode == 4:
                    input('cumulative galaxy score distribution')

        if self.showmode == 3: input('show all plots')

    def send(self):

        # final report before telescopes scheduling
        self.figlist = []
        try:
            self.fig_2d.savefig(self.sname.replace('$telname$_','').\
                                replace('.txt','_1.png'))
            self.figlist.append(self.sname.replace('$telname$_','').\
                                replace('.txt','_1.png'))
        except: pass
        try:
            self.fig_gd.savefig(self.sname.replace('$telname$_','').\
                                replace('.txt','_2.png'))
            self.figlist.append(self.sname.replace('$telname$_','').\
                                replace('.txt','_2.png'))
        except: pass
        try:
            self.fig_gl.savefig(self.sname.replace('$telname$_','').\
                                replace('.txt','_3.png'))
            self.figlist.append(self.sname.replace('$telname$_','').\
                                replace('.txt','_3.png'))
        except: pass
        try:
            self.fig_cum.savefig(self.sname.replace('$telname$_','').\
                                 replace('.txt','_4.png'))
            self.figlist.append(self.sname.replace('$telname$_','').\
                                replace('.txt','_4.png'))
        except: pass

        # - email
        if eval(self.optlist['arg']['email']['activate']):
            for _toaddress in self.optlist['arg']['email']['to'].split(','):
                _sent = pst.sendemail(self.optlist['arg']['email']['from'],\
                            self.optlist['arg']['email']['pwd'],\
                            self.optlist['arg']['email']['smtp'],\
                            self.optlist['arg']['email']['sub'],\
                            self.optlist['arg']['email']['from'],\
                            _toaddress,self.str1,self.figlist,self.filelist)
                if self.verbose:
                    if _sent: print ('email sent successful to %s'%_toaddress)
                    else: print ('email sent failed to %s'%_toaddress) 

        # - slack
        if eval(self.optlist['arg']['slack']['activate']):
            for _toaddress in self.optlist['arg']['slack']['channel'].split(','):
                _sent = pst.slack(self.optlist['arg']['slack']['token'], \
                        _toaddress, self.str3, self.figlist+self.filelist,\
                        _msg=False)
                if self.verbose:
                    if _sent: print ('slack sent successful to %s'%_toaddress)
                    else: print ('slack sent failed to %s'%_toaddress) 

        # - SMS
        if eval(self.optlist['arg']['phone']['activate']):
            for _usr in self.optlist['arg']['phone']['to'].split(','):
                _sms = pst.phone(self.optlist['arg']['phone']['account'],\
                                 self.optlist['arg']['phone']['token'],\
                                 self.optlist['arg']['phone']['from'],\
                                 _usr,self.str2)
                if self.verbose:
                    if _sms: print ('SMS sent successful to %s'%_usr)
                    else: print ('SMS sent failed to %s\tno twilio'%_usr) 

        # - excute python codes
        #   for special functions, e.g. insert to DB
        for _tt in self.optlist:
            if _tt in ['arg','tmp']:continue
            if len(self.optlist[_tt]['scheduler']['py2'])>0:
                import subprocess
                for _py in self.optlist[_tt]['scheduler']['py2'].split(','):
                    _ff = self.schfile[self.optlist[_tt]['telescope']['name']]
                    if hasattr(self,'graceid'):_id=self.graceid
                    else:_id='Ttrigger'
                    pid = subprocess.Popen(['python',_py,_ff,_id,\
                            self.optlist[_tt]['observe']['limalt'],\
                            self.optlist[_tt]['scheduler']['dither'],\
                            self.optlist[_tt]['scheduler']['nexp'],\
                            self.optlist[_tt]['telescope']['exptime']],\
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    output,error = pid.communicate()
                    if self.verbose: 
                        if len(output)>0: print('### Info: %s'%output)
                        if len(error)>0: print('### Error: %s'%error)

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
            _flist.remove(ff0)
    return _flist

def decomposit(_list):

    if not type(_list) in [list, np.ndarray]:return _list
    _d = False
    while not _d:
        _lenlist = []
        for _l in _list:
            if not type(_l) in [list, np.ndarray]:
                _d = True
            else:
                _lenlist.append(len(_l))
        if len(_list) == len(_lenlist):
            _u = list(np.unique(_lenlist))
            if len(_u) == 1 and _u[0] == 1:
                _list = [ii[0] for ii in _list]
            else: _d = True
    return _list
