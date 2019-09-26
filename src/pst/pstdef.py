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

from astropy.io import fits
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
import scipy.stats
import pst

_pstpath = pst.__path__[0]

# for LVC GW only
# !!! for GRB, TBD
# from XML:
_strxml = ['GraceID', 'AlertType', 'Group', 'FAR', \
           'Terrestrial', 'HasNS', 'Instruments', 'EventPage' ]
# from fits header:
_strhdr = ['DISTMEAN', 'DISTSTD', 'DATE-OBS', \
           'MJD-OBS', 'CREATOR']
# decrease nside in order to be faster
_nside = 64 # for computer contour lines
_nside1 = 12 # for check dist range in each piece of sky

# for check dist range in each piece of sky
distmin,distmax,delta,frac = 0,10000,100,1e-1

#################################################

class main(object):
    
    def __init__(self,optlist):
    
        # Main procedure start
        start_time = time.time()
        
        # try cleaning all the plots if any
        plt.close('all')

        # read params        
        self.optlist = optlist
        self.str1 = self.optlist['arg']['email']['emailcontent']
        self.str2 = self.optlist['arg']['phone']['phonecontent']
        self.str3 = self.optlist['arg']['slack']['slackcontent']
        self.distdict, self.indict = {}, {}

        # define time
        self.def_time()

        # decide figure out type
#        self.dec_opt()

        # define nside
        self.nside = int(self.optlist['arg']['priorization']["nside"])        
        if len(self.optlist['tmp']['tmap'])>0:
            self.nside = int(np.sqrt(len(self.optlist['tmp']['tmap'])/12))        

        ''' 1 - Check triggers (hp map) '''
        self.triggers()

        # decide telescope schedule
        # if trigger available
        self.dec_scheduler()

        ''' 2 - Obtain galaxies '''
        self.galaxies()       

        ''' 3 - generate pointings 
        for tiling search '''
        self.pointings()

        ''' 4 - Priorization algorithm 
        for pointings/galaxies '''
        self.ranking()

        ''' 5 - Visualization '''
        self.visualization()

        ''' 6 - distributing to users/api/... '''
        self.send()

        ## done
        if self.optlist['arg']['show']['verbose']: 
            print("%Finished in %i secs"%\
                  (int(time.time()-start_time)))

    def def_time(self): # define time
        if self.optlist['arg']['observe']['obstime'] == 'now': 
            self.timenow = astropy.time.Time.now()            
        else:
            try: 
                self.timenow = astropy.time.Time(self.optlist['arg']['observe']['obstime'], \
                        scale='utc')
            except: 
                self.timenow = astropy.time.Time.now() + \
                        astropy.time.TimeDelta(float(optlist['arg']['observe']['obstime'])*60, \
                        format='sec')           
        self.jdnow = self.timenow.jd

    def dec_scheduler(self):
        # decide which strategy to use
        self.tellist, self.scheduler = [], []
        for ff in self.optlist:
            if ff in ['arg','tmp']:continue
            self.tellist.append(self.optlist[ff])
        area = np.array([len(self.indexlist[ii]) for ii in self.indexlist])            
        for _tel in self.tellist:               
            if eval(self.optlist['arg']['priorization']['trigger']): # if trigger                       
                if _tel['pointings']['scheduler'] == 'A':                  
                    _nn = max(area*self.areasingle)/\
                          float(_tel['telescope']['fovw'])/\
                          float(_tel['telescope']['fovh'])
                    if _nn <= 1000:self.scheduler.append('T')
                    else:self.scheduler.append('G')
            else: 
                if _tel['pointings']['scheduler'] in ['G','T']:
                    self.scheduler.append(_tel['pointings']['scheduler'])
                else: sys.exit('wrong scheduler input')        

    def dec_opt(self):
        ## decide savefig, plot, or no fig
        # pmet: 1.save figure and send via email/slack, will not pause
        #       2.pause and show plot whenever there's one available
        #       3.interactive show: pause and interact with plots
        #       4.no figure
        if self.optlist['tmp']['search'] == 'auto': 
            # for auto sarch: either 1 or 4
            if eval(self.optlist['arg']['show']["verbose"]): self.pmet=1
            else: self.pmet=4
        if self.optlist['tmp']['search'] in ['trigger', 'galaxy', 'normal']:
            # for manual search:
            if eval(self.optlist['arg']['show']["verbose"]):
                if eval(self.optlist['arg']['show']["interactive"]):self.pmet=3
                else: self.pmet=2
            else: self.pmet=4
        if self.pmet in [2,3]: pl.ion()
        if len(self.optlist['arg']['show']["showmap"])>0:
            self.showmap = [int(gg) for gg in \
                    self.optlist['arg']['show']["showmap"].split(',')]
        else: self.showmap = []

    def trigger_validation(self):
        # validate a trigger by using the header of fits
        ''' judge interests of the trigger '''

        # read more from the xml, complementry to the fits header                
        for ii in _strxml:
            if ii in self.optlist['tmp']['voevent']:
                self.str1 += '#\t%s:\t%s\n'%(ii,self.optlist['tmp']['voevent'][ii])
                self.str2 += '%s:%s, '%(ii,self.optlist['tmp']['voevent'][ii])
                self.str3 += '%s:%s, '%(ii,self.optlist['tmp']['voevent'][ii])
                self.indict[ii] = self.optlist['tmp']['voevent'][ii]
            else:
                self.str1 += '#\t%s:\tNot available\n'%ii
                self.str2 += '%s:None, '%ii
                self.str3 += '%s:None, '%ii
                self.indict[ii] = None

        # show params of the fits header
        for ii in _strhdr:
            if ii in self.header:            
                self.str1 += '#\t%s:\t%s\n'%(ii,self.header[ii])
                self.str2 += '%s:%s, '%(ii,self.header[ii])
                self.str3 += '%s:%s, '%(ii,self.header[ii])
                self.indict[ii] = self.header[ii]
            else:
                self.str1 += '#\t%s:\tNone\n'%ii
                self.str2 += '%s:None, '%ii
                self.str3 += '%s:None, '%ii
                self.indict[ii] = None

        # check the area of skymap     
        # get area in single pixel, unit in deg
        if False:self.areasingle = (hp.nside2resol(hp.get_nside(self.optlist['tmp']['tmap']), \
                                        arcmin=True)/60.)**2
        else:self.areasingle = (hp.nside2resol(_nside,arcmin=True)/60.)**2

        self.area={}
        for _cc in self.indexlist: # levels: .9,.5,....
            _num = len(self.indexlist[_cc])          
            if _num>0:
                self.str1 += ' %.2fsky: %.f sq deg\n'%(_cc,_num*self.areasingle)
                self.str2 += '%.2fsky: %.f sq deg, '%(_cc,_num*self.areasingle)
                self.str3 += '%.2fsky: %.f sq deg, '%(_cc,_num*self.areasingle)
                self.distdict[_cc] = ' %.f%% sky:%.f $deg^2$\n'%(_cc*100,_num*self.areasingle)
                self.indict['%.2fsky'%_cc] = _num*self.areasingle
                self.area[_cc] = _num*self.areasingle               
            else:
                self.str1 += ' %.2f%sky:None\n'%(_cc)
                self.str2 += '%.2fsky:None, '%(_cc)
                self.str3 += '%.2fsky:None, '%(_cc)
                self.distdict[_cc] = ' %.f%% sky:None\n'%(_cc*100)
                self.indict['%.2fsky'%_cc] = None
                self.area[_cc] = None
        self.str1 += '\n'
        self.str2 = self.str2[:-2]
        self.str3 = self.str3[:-2]

    def triggers(self):
        # check if need trigger

        if eval(self.optlist['arg']['priorization']['trigger']): # read if yes

            # read and print some values from the FITS header.
            self.header = dict(self.optlist['tmp']['header'])

            # obtain containment contour of GW PDF
            probs = hp.pixelfunc.ud_grade(self.optlist['tmp']['tmap'], 
                    _nside) #reduce nside to make it faster
            probs = probs/np.sum(probs)

            levels = [float(kk) for kk in \
                self.optlist['arg']['show']["contours"].split(',')]
            self.theta_contour, self.phi_contour, self.indexlist = \
                                pst.compute_contours(levels,probs)           

            # validate skymap
            self.trigger_validation()
           
            # immediately report before telescopes scheduling
            # - email
            if eval(self.optlist['arg']['email']['sendemail']):
                for _toaddress in self.optlist['arg']['email']['emailto'].split(','):
                    _sent = pst.sendemail_1(self.optlist['arg']['email']['email'],\
                                self.optlist['arg']['email']['emailpass'],\
                                self.optlist['arg']['email']['emailsmtp'],\
                                self.optlist['arg']['email']['emailsub'],\
                                self.optlist['arg']['email']['email'],\
                                _toaddress,self.str1+\
                                '\n\nimmediate report, scheduler still ongoing...')
            # - slack
            if eval(self.optlist['arg']['slack']['activate']):
                for _usr in self.optlist['arg']['slack']['channel'].split(','):
                    _slack = pst.slack(self.optlist['arg']['slack']['slack_bot_token'], \
                                       _usr, self.str3)
                    if _slack: print ('slack sent successful to %s'%_usr)
                    else: print ('slack sent failed to %s\tno slackclient'%_usr) 

            # - SMS
            if eval(self.optlist['arg']['phone']['activate']):
                for _usr in self.optlist['arg']['phone']['to'].split(','):
                    _sms = pst.phone(self.optlist['arg']['phone']['account'],\
                                     self.optlist['arg']['phone']['token'],\
                                     self.optlist['arg']['phone']['from'],\
                                     _usr,self.str2)
                    if _sms: print ('SMS sent successful to %s'%_usr)
                    else: print ('SMS sent failed to %s\tno twilio'%_usr) 

            # - insert to DB               
            if False: #eval(optlist['arg']['database']['activate']):
                import sqlconn
                sqlconn.insert_values('ligoevents',indict)
        else: self.theta_contour, self.phi_contour = None, None

        """
        if 1 in self.showmap:
            map2show = self.optlist['tmp']['tmap']
        else:
            map2show = np.zeros(12*self.nside**2)

        # healpix show of trigger
        ''' fignum 1: 2d sky'''
        # parameters for plotting
        pparams = {'hpmap':map2show,\
                   'title':self.optlist['arg']['show']['title'],\
                   'theta':self.optlist['arg']['show']["theta"],\
                   'phi':self.optlist['arg']['show']["phi"],\
                   'fignum':1,'ordering':self.optlist['arg']['show']["ordering"],\
                   'coord':self.optlist['arg']['show']["coord"],\
                   'norm':self.optlist['arg']['show']["norm"],\
                   'figsize':self.optlist['arg']['show']["figsize"],\
                   'min':self.optlist['arg']['show']["min"],\
                   'max':self.optlist['arg']['show']["max"],\
                   'theta_contour':self.theta_contour,\
                   'phi_contour':self.phi_contour,\
                   'label':self.optlist['arg']['show']["label"],\
                   'colors':self.optlist['arg']['show']["colors"],\
                   'distinfo':self.distdict,'tellist':self.tellist,\
                   'timenow':self.timenow}
        # parameters that can be changed via interactive mode
        optparams = ['theta','phi','min','max','norm',\
                     'title','figsize','ordering','coord']

        _pm = False
        if self.pmet==3: 
            _pm, fig_2d = pst.interactive_show(pst.mollview,pparams,optparams)            
        elif self.pmet in [1,2]:            
            fig_2d = pst.mollview(pparams)
            if self.pmet == 2:input('show trigger in skymap')
        if _pm: # update show options
            for _cc in optparams:self.optlist['arg']['show'][_cc]=_pm[_cc]
        """

    def galaxies(self):
        # check if need galaxies

        if eval(self.optlist['arg']["priorization"]["mass"]) or \
           eval(self.optlist['arg']["priorization"]["dist"]) or \
           eval(self.optlist['arg']["priorization"]["number"]) or \
           'G' in self.scheduler:

            # query galaxies via vizier
            if int(self.optlist['arg']['galaxies']["catalog"])==1:self.catname='GLADE'
            elif int(self.optlist['arg']['galaxies']["catalog"])==2:self.catname='GWGC'
            else:self.catname='???'

            limdist,limrag,limdecg,limmag = \
            [float(zz) for zz in self.optlist['arg']['galaxies']["limdist"].split(',')],\
            [float(zz) for zz in self.optlist['arg']['galaxies']["limra"].split(',')],\
            [float(zz) for zz in self.optlist['arg']['galaxies']["limdec"].split(',')],\
            [float(zz) for zz in self.optlist['arg']['galaxies']["limmag"].split(',')]

            if self.optlist['tmp']['distmu'] is not None and \
               self.optlist['tmp']['distsigma'] is not None and \
               self.optlist['tmp']['distnorm'] is not None:
                # if 3D GW loc available
                # query the min to max distance from catalog
                # meanwhile, get dist range for every piece of sky 
                # (use small nside, since it would not change too much)
                # the dist range of each piece will be then included,
                # together with selected galaxies (gaussian), calculated
                # for a score, contributing to the final list                
                self.tpiece,self.ppiece = hp.pix2ang(_nside1, np.arange(12*_nside1**2))                
                self.dminlist, self.dmaxlist = pst.gwdist(self.nside,\
                    self.optlist['tmp']['distmu'], self.optlist['tmp']['distsigma'], \
                    self.optlist['tmp']['distnorm'], self.tpiece,self.ppiece)
            else: 
                self.dminlist, self.dmaxlist = np.array([min(limdist)]), \
                                               np.array([max(limdist)])
            dmin = min(self.dminlist[np.where(self.dminlist!=-99)])
            dmax = max(self.dmaxlist[np.where(self.dmaxlist!=-99)])

            self.gid,self.gname,self.gra,self.gdec,self.gmag,self.gdist = \
                    pst.galaxies(catalog=int(self.optlist['arg']['galaxies']["catalog"]),\
                    limra=limrag, limdec=limdecg, limdist=[dmin,dmax],\
                    size=int(self.optlist['arg']['galaxies']["size"]),\
                    limmag=limmag,filtro=self.optlist['arg']['galaxies']["filter"],\
                    verbose = eval(self.optlist['arg']['show']["verbose"]),\
                    cachemode=int(self.optlist['arg']['galaxies']['cachemode']),\
                    cachefile=self.optlist['arg']['galaxies']['cachefile'])
         
            """
            # healpix show of trigger
            ''' fignum 1: 2d sky
                fignum 2: dist
                fignum 3: lums
            ''' 
            # 1 - 2D distribution
            # focused, no interactive mode
            pparams = {'ra':self.gra,'dec':self.gdec,\
                       'theta':self.optlist['arg']['show']["theta"],\
                       'phi':self.optlist['arg']['show']["phi"],\
                       'fignum':1,'color':'r',\
                       'coord':self.optlist['arg']['show']["coord"],\
                       'label':'%s galaxies(%i)'%(self.catname,len(self.gra))}
            if 2 in self.showmap:
                fig_2d = pst.pointview(pparams)
                if self.pmet in [2,3]:
                    input('show galaxies in skymap')

            # 2 - distance distribution
            pparams = {'distmin':dmin,'distmax':dmax,'dist':self.gdist,'fignum':2,\
                       'color1':'k','color2':'r','scale':'linear','nbin':10,\
                       'label':'%s galaxies(%i)'%(self.catname,len(self.gra))}
            optparams = ['distmin','distmax','nbin','color1','color2','scale']
            if self.pmet==3 and 3 in self.showmap:
                _pm, fig_gd = pst.interactive_show(pst.distview,pparams,optparams)
            elif self.pmet in [1,2] and 2 in self.showmap:
                fig_gd = pst.distview(pparams)    
            if self.pmet==2 and 3 in self.showmap:
                input('galaxy distance distribution')

            # 3 - liminosity distribution
            pparams = {'distmin':dmin,'distmax':dmax,'mag':self.gmag,\
                       'dist':self.gdist,'fignum':3,'scale':'linear',\
                       'color1':'r','color2':'grey','nbin':1,\
                       'label':'%s galaxies(%i)'%(self.catname,len(self.gra))}
            optparams = ['distmin','distmax','nbin','color1','color2','scale']
            if self.pmet==3 and 3 in self.showmap:
                _pm, fig_gl = pst.interactive_show(pst.lumsview,pparams,optparams)
            elif self.pmet in [1,2] and 3 in self.showmap:
                fig_gl = pst.lumsview(pparams)
            if self.pmet==2 and 3 in self.showmap:
                input('galaxy luminosity distribution')

            '''
            # convolve galaxy info into score
            # from arcsec to radians
            try:self.gradius = float(optlist['arg']['galaxies']["radius"])*2*mt.pi/60/360
            except:self.gradius=False

            if eval(self.optlist['arg']["priorization"]["mass"]): # lums
                self.lums = 10**((-1)*(self.gmag/2.5))
                self.glumsmap = pst.make_hpfitsmap(self.gra,self.gdec,\
                                self.lums,self.nside,self.gradius)               
                # plot
                pparams = {'hpmap':self.glumsmap,'title':'galaxy dist',\
                   'theta':self.optlist['arg']['show']["theta"],\
                   'phi':self.optlist['arg']['show']["phi"],\
                   'fignum':1,'ordering':self.optlist['arg']['show']["ordering"],\
                   'coord':self.optlist['arg']['show']["coord"],\
                   'norm':self.optlist['arg']['show']["norm"],\
                   'figsize':self.optlist['arg']['show']["figsize"],\
                   'min':self.optlist['arg']['show']["min"],\
                   'max':self.optlist['arg']['show']["max"],\
                   'theta_contour':self.theta_contour,\
                   'phi_contour':self.phi_contour,\
                   'label':self.optlist['arg']['show']["label"],\
                   'colors':self.optlist['arg']['show']["colors"],\
                   'distinfo':self.distdict,'tellist':self.tellist,\
                   'timenow':self.timenow}
                optparams = ['theta','phi','min','max','norm',\
                             'title','figsize','ordering','coord']
                _pm = False
                if self.pmet==3: 
                    _pm, fig_2d = pst.interactive_show(pst.mollview,pparams,optparams)            
                elif self.pmet in [1,2]:            
                    fig_2d = pst.mollview(pparams)
                    if self.pmet == 2:input('show galaxies map')
                if _pm: # update show options
                    for _cc in optparams:self.optlist['arg']['show'][_cc]=_pm[_cc]
            '''
            """

    def pointings(self):
        # generate pointings
        _idc,_ral,_decl,fig_np = pst.pointings(ralist=_ralist,declist=_declist,\
                        limdec=limdecp,limra=limrap,fovh=float(optlist[_tt]['telescope']['fovh']),\
                        fovw=float(optlist[_tt]['telescope']['fovw']),\
                        obx=int(optlist[_tt]['pointings']['ob'].split('*')[0]),\
                        oby=int(optlist[_tt]['pointings']['ob'].split('*')[1]),\
                        shifth=0.,shiftw=0.)


#    def ranking(self):
#            if eval(self.optlist['arg']["priorization"]["dist"]): # dist
#                ss
#            if eval(self.optlist['arg']["priorization"]["number"]): # counts
#                ss
#        self.pmap = np.zeros(hp.nside2npix(self.nside))+1
        # cross map with trigger map
#        self.pmap *= self.optlist['tmp']['tmap']


def choose(_dicti):
    done, _dict = False, _dicti   
    while not done:        
        answ = input('what to show: %s or [a]ll or [o]uter layer or [q]uit?'%_dict.keys())
        if answ == 'q': done=True
        elif answ == 'o':_dict = _dicti
        elif answ in _dict.keys():
            try:os.system('clear')
            except:pass
            print ('>'*5,'go to %s'%answ)            
            _dict = _dict[answ]                 
        elif answ == 'a':
            for _key in _dict:              
                print ('\t--> %s : %s \n'%(_key,_dict[_key]))
        else: 
            try:os.system('clear')
            except:pass
            print ('Error: wrong input...')
        try: _dict.keys()
        except: 
            print ('\t %s'%_dict)
            done = True

def build_hp_map(v,mapname,nside,_coord='C'):
    # generate healpix fits map
    nside = int(nside)

    # read infos
    try:
        _nra,_ndec,_loc = float(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Value2.C1),\
                          float(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Value2.C2),\
                          float(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Error2Radius)
    except:
        print ('### Error: no ra,dec,radius found in voevent!!!')
        return [],None

    if _nra == 0 and _ndec == 0:
        print ('### Error: no ra,dec reported in voevent!!!')
        return [],None

    if _loc == 0:
        print ('### Error: no radius reported in voevent!!!')
        return [],None
                
    try:
        _timeobs = str(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Time.TimeInstant.ISOTime)        
    except:
        print ('### Error: no timeobs found in voevent!!!')
        return [],None

    _timeobs1 = astropy.time.Time(_timeobs, format='isot', scale='utc')
    _mjdobs = _timeobs1.mjd

    try:_object = v.Why.Inference.Name
    except:_object = 'unKnown'           

    gradius = _loc*2*mt.pi/360 # from deg to radians
    _pmap = np.zeros(hp.nside2npix(nside))
    _index = pst.DeclRaToIndex(_ndec,_nra,nside)
    _pmap[_index]+=1
    _pmap=hp.sphtfunc.smoothing(_pmap,fwhm=gradius)
        
    _pmap = _pmap/sum(_pmap)       
    hlist = [('CREATOR','PSTOOLS'),
             ('OBJECT',_object),
             ('NSIDE',nside),
             ('MJD-OBS',_mjdobs),
             ('DATE-OBS',_timeobs)]
    #
    hp.write_map(mapname,_pmap, coord=_coord, extra_header=hlist, overwrite=True)
    return _pmap, hlist

def get_hp_map(_fits,verbose=False,nside=1024):
    # get healpix fits

    params = []
    try:   # if it's healpix format
        try:  # read 3D trigger healpix map (for GW)            
            (tmap, distmu, distsigma, distnorm), header = hp.read_map(_fits, \
                field=[0, 1, 2, 3],h=True, verbose=verbose)
        except: # !!!!  need to be tested for the other trigger types: GRB/AMON/ect
            tmap, header = hp.read_map(_fits, h=True, verbose=verbose)
            distmu, distsigma, distnorm = None, None, None

    except:  # if it's xml format, which can be used to extract the fits url
        try: import voeventparse
        except: sys.exit('### Error: pls installed voenet-parse first...')
        try: 
            with open(_fits, 'rb') as f:root = voeventparse.load(f)
        except: 
            sys.exit('### Error: %s not readable. Why:\n'%_fits+\
                     '\t1. healpix header missing END card;\n'+\
                     '\t2. not healpix fits;\n'+\
                     '\t3. not xml including fits url')
        params = {elem.attrib['name']:
                  elem.attrib['value']
                  for elem in root.iterfind('.//Param')}
        try: mapurl = params['skymap_fits']
        except: mapurl = False
        if mapurl:   # download map
            try: _fits = pst.get_skymap(mapurl,os.path.basename(root.attrib['ivorn']))
            except: sys.exit('### Error: no ivorn in %s, TB checked'%_fits)
            try:
                (tmap, distmu, distsigma, distnorm), header = hp.read_map(_fits, \
                    field=[0, 1, 2, 3],h=True, verbose=verbose)
            except:
                tmap, header = hp.read_map(_fits, h=True, verbose=verbose)
                distmu, distsigma, distnorm = None, None, None
        else:    # generate map
            print ('### Warning: no skymap_fits in %s, try bulding...'%_fits)
            tmap, header = pst.build_hp_map(root,os.path.basename(_fits).replace('.xml','.fits'),\
                    nside,_coord='C')
            
            if len(tmap) > 0:
                print ('### Warning: genearted a fits,'+\
                       ' %s'%os.path.basename(_fits).replace('.xml','.fits'))
                distmu, distsigma, distnorm = None, None, None
            else:sys.exit('### Error: failed to build fits')
    return tmap, distmu, distsigma, distnorm, header, params

def gwdist(nside, distmu, distsigma, distnorm, thetal, phil):

    rminlist,rmaxlist = [],[]
    for theta,phi in zip(thetal,phil):       
        ipix = hp.ang2pix(nside, theta, phi)
        r = np.linspace(distmin,distmax,num=delta)   
        dp_dr = r**2 * distnorm[ipix] * scipy.stats.norm(\
                    distmu[ipix], distsigma[ipix]).pdf(r)
        rl = r[np.where(dp_dr>frac*max(dp_dr))]        
        if len(rl)>0: 
            rminlist.append(min(rl))
            rmaxlist.append(max(rl))
        else: 
            rminlist.append(-99)
            rmaxlist.append(-99)
    return np.array(rminlist), np.array(rmaxlist)

def get_skymap(skymap_url,graceid='Stest',_dir='/tmp/'):
    """
    Look up URL of sky map in VOEvent XML document,
    download sky map, and parse FITS file.
    """  

    filename = _dir + graceid + '_' + os.path.basename(skymap_url)
    try: 
        import wget
        opt = 1
        print ('using wget for downloading')
    except:
        import requests,tempfile,shutil
        opt = 2
        print ('no wget, using requests for downloading')

    if opt == 1:
        # wget
        import wget
        if os.path.exists(filename):os.remove(filename)
        wget.download(skymap_url, out=filename)

    elif opt == 2:
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

    # Done!
    return filename

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
