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
#from wxpy import Bot
#import sqlconn

from astropy.io import fits
from astropy.table import Table
import astropy.coordinates
import astropy.time
import astropy.units as u
from astroquery.vizier import Vizier
import scipy.stats
#from slackclient import SlackClient

#from scp import SCPClient
#from pst import scheme,priorization,pstplot,scheduler,configure,link
import pst
_pstpath = pst.__path__[0]

#################################################

class main(object):
    
    def __init__(self,mapname,optlist):
    
        # Main procedure start
        start_time = time.time()
        
        # try cleaning all the plots if any
        plt.close('all')

        # read params
        self.mapname = mapname
        self.optlist = optlist

        # define plots setting
        self.colorlist, self.ncolor = ['b','g','k','y','c','m'], 0
        
        # define time
        self.def_time()

        # decide figure out type
        self.dec_opt()

        ''' 1- Priorization generation algorithm'''
        # define a blanket Priorization map
        self.nside = int(self.optlist['arg']['priorization']["nside"])
        self.pmap = np.zeros(hp.nside2npix(self.nside))+1

        ''' 1.1 - trigger Priorization
                   include trigger prob if any '''
        self.tprocess()

        ''' 1.2 - mass Priorization 
                   include galaxies' effect if any '''

        ''' 2 - generate pointings       
        For different strategies: one for tiling (if any) and one for galaxy (if any)
        !!!! important:
        for tiling, valid only for multiple telescopes with same or similar FoV
        if not, try to use OB mode to make them have the same FoV
        Otherwise, I stopped'''

        ''' 3 - Ranking pointings '''
        
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

    def dec_opt(self):
        ## decide savefig, plot, or no fig
        # pmet: 1.if send email: savefig for specific plots
        #       2.normal verbose: input for all plots
        #       3.inter verbose: inter_plot for all plots
        #       4.no fig
        if self.optlist['arg']['search'] == 'auto':
            if eval(self.optlist['arg']['plot']["verbose"]): self.pmet=1
            else: self.pmet=4
        if 'man' in self.optlist['arg']['search']:
            if eval(self.optlist['arg']['plot']["verbose"]):
                if eval(self.optlist['arg']['plot']["interactive"]):self.pmet=3
                else: self.pmet=2
            else: self.pmet=4
        if eval(self.optlist['arg']['email']["sendemail"]) or \
           eval(self.optlist['arg']['wechat']["activate"]) or \
           eval(self.optlist['arg']['phone']['activate']):  self.pmet=1

        if self.pmet in [2,3]: pl.ion()
        if len(self.optlist['arg']['plot']["showmap"])>0:
            self.showmap = self.optlist['arg']['plot']["showmap"].split(',')
        else: self.showmap = []

    def tprocess(self):
        if eval(self.optlist['arg']['priorization']['trigger']):

            # in case there're trigger healpix map input
            if not self.mapname is None:
                            
                try:  
                    # read 3D trigger healpix map (for GW)
                    # need to be tested for the other trigger types: GRB/AMON/ect
                    (self.maptrigger, self.distmu, self.distsigma, \
                     self.distnorm), self.header = \
                        hp.read_map(self.mapname, field=[0, 1, 2, 3], h=True, \
                        verbose=eval(self.optlist['arg']['plot']["verbose"]))
                except:
                    # report there's an error in the process of reading healpix map
                    if eval(optlist['arg']['slack']['activate']):
                        # send message via phone immediately
                        slack_client = SlackClient(optlist['arg']['slack']['slack_bot_token'])
                        for _channel in optlist['arg']['slack']['channel'].split(','):
                            slack_client.api_call("chat.postMessage", channel=_channel,
                                text='new alert, however, wrong healpix fits format!!', \
                                as_user=True)  
                    logging.info('wrong fits format for healpy!!') 
                    return

            # read and print some values from the FITS header.
            header = dict(header)

            # validate skymap
            _ilist,_alist,_ld,_str,_date_obs,_mjd_obs = trigger_validation(maptrigger, header)
            optlist['arg']['email']['emailcontent'] += _str

            # read areas: 50 68 90 99
            index1,index2,index3,index4=_ilist.values()
            _a50,_a68,_a90,_a99=_alist.values()                      

            # voevent params
            if 'voevent' in optlist['arg'] and header['CREATOR']!='PSTOOLS':
                # for LVC
                params = optlist['arg']['voevent']

                if True:
                    if optlist['arg']['email']['role'] == 'test':_trigger = False
                    else:_trigger = True
                else:_trigger = True
                    
                # for LVC: params of VOEvent to show in email
                _strremain = ['Terrestrial', 'Group', 'EventPage', 'GraceID', \
                              'HasRemnant', 'skymap_png', 'skymap_fits', 'AlertType', \
                              'HasNS', 'BBH', 'BNS', 'Instruments', 'NSBH']
                for ii in params:
                    if ii in _strremain:optlist['arg']['email']['emailcontent']+='#\t%s:\t%s\n'%(ii,params[ii])

                # FAR
                FAR = params['FAR']
                FAR = float(FAR)*360*24*3600 # Hz to per year
                optlist['arg']['email']['emailcontent']+='#\tFAR (per year):\t%.4e\n'%FAR

                # params for alert validation
                grace_id = params['GraceID']
                mapurl = params['skymap_fits']
                BBH = float(params['BBH'])
                BNS = float(params['BNS'])
                HasNS = float(params['HasNS'])
                HasRemnant = float(params['HasRemnant'])
                NSBH = float(params['NSBH'])
                Terrestrial = float(params['Terrestrial'])
                Group = params['Group']
                EventPage = params['EventPage']
                Instruments = params['Instruments']

                #The first selection criteria are about the likelyhood that the event is real and astrophysical. So we will take into consideration only events with a False Alarm Rate (FAR) < 1 event per year. Then we will exclude all the events for which the probability of being a terrestrial event is equal or higher than 90%. 
                
                if FAR < 1 and Terrestrial<0.9:_trigger *= True
                else:_trigger *= False
                
                # !!! WHEN TO TRIGGER ???
                ''' For VST:
                PROB_NS > 0.1
                90% skymap area < 100 deg2
                or
                - PROB_NS < 0.1
                - 90% skymap area < 50 deg2
                - Distance < 100 Mpc
                '''
                if float(HasNS)>.1 and float(HasRemnant)>.1:
                    # interesting for all BNS NS-BH
                    _trigger*=True     
                elif float(HasNS)<.1 and _a90<200:
                    if not _ld is None:
                        Dmean = float(_ld.split('+/-')[0])
                        if Dmean<500:_trigger*=True                
                        else:_trigger*=False
                    else:_trigger*=True    
                else:_trigger*=False

                # mysql insert
                if eval(optlist['arg']['database']['activate']):
                    if _ld is None:Dmean,Dvar='None','None'
                    else:Dmean,Dvar=_ld.split('+/-')[0],_ld.split('+/-')[1]
                    if 'Preliminary' in optlist['arg']['email']['emailcontent']:_stage='Preliminary'
                    elif 'Initial' in optlist['arg']['email']['emailcontent']:_stage='Initial'
                    elif 'Update' in optlist['arg']['email']['emailcontent']:_stage='Update'
                    else:_stage='Null'

                    sqlconn.insert_values('ligoevents',\
                                          {'GraceID':grace_id,\
                                           'JD':_mjd_obs,\
                                           'stage':_stage,\
                                           'type':Group,\
                                           'Dmean':Dmean,\
                                           'Dvar':Dvar,\
                                           'FAR':FAR,\
                                           'BNS':BNS,\
                                           'BBH':BBH,\
                                           'NSBH':NSBH,\
                                           'HasNS':HasNS,\
                                           'HasRemnant':HasRemnant,\
                                           'Terrestrial':Terrestrial,\
                                           'time':str(_date_obs),\
                                           'loc50':str(_a50),\
                                           'loc68':str(_a68),\
                                           'loc90':str(_a90),\
                                           'loc99':str(_a99),\
                                           'detector':Instruments,\
                                           'url':mapurl
                                       })
                  
                # email subject
                if _trigger:  _comments = '[interesting]'
                else:  _comments = '[not interesting]'
                optlist['arg']['email']['emailsub'] += _comments

                # slack
                if eval(optlist['arg']['slack']['activate']):
                    # send message via phone immediately   
                    optlist['arg']['phone']['phonecontent'] += '`%s`\n a %s event detected at %sUT\n *FAR=%.2e per year, HasNS=%s, HasRemnant=%s, localization=%.2f-%.2f-%.2f-%.2f (50-68-90-99) sq deg, Distance=%s Mpc, Terrestrial=%.5e*\nGraceDB: %s\ntelescope scheduler working...'%(grace_id,Group,_date_obs,FAR,HasNS,HasRemnant,_a50,_a68,_a90,_a99,_ld,Terrestrial,EventPage)
                    slack_client = SlackClient(optlist['arg']['slack']['slack_bot_token'])
                    for _channel in optlist['arg']['slack']['channel'].split(','):
                        slack_client.api_call("chat.postMessage", channel=_channel,
                                        text=optlist['arg']['phone']['phonecontent'], as_user=True)

                # send phone message
                if eval(optlist['arg']['phone']['activate']) and \
                   optlist['arg']['email']['role'] == 'observation':                       
                    optlist['arg']['phone']['phonecontent'] += '%s, a %s event detected at %sUT %s: FAR=%.2e per year. HasNS=%s. HasRemnant=%s. localization (90 per)=%.2f sq degree. Distance=%s Mpc. Terrestrial=%.5e. %s'%(grace_id,Group,_date_obs,_comments,FAR,HasNS,HasRemnant,_a90,_ld,Terrestrial,EventPage)
                    for _to in optlist['arg']['phone']['to'].split(','):
                        if eval(optlist['arg']['phone']['activate']):pst.link.phone(optlist['arg']['phone']['account'],optlist['arg']['phone']['token'],optlist['arg']['phone']['from'],_to,optlist['arg']['phone']['phonecontent'])

            elif 'voevent' in optlist['arg'] and header['CREATOR']=='PSTOOLS':

                # for other alerts
                params = optlist['arg']['voevent']
                for ii in params:optlist['arg']['email']['emailcontent']+='#\t%s:\t%s\n'%(ii,params[ii])

            else:
                # update graceid via header: MS190402o -> M328682
                # choose whatever you want
                try:grace_id = header['OBJECT']
                except:grace_id='GWxyz'

            # update nside
            if nside != int(header['NSIDE']):
                nside = int(header['NSIDE'])
                _info = 'Warning: for trigger search, change nside to %i, according to trigger healpix map!'%nside
                if eval(optlist['arg']['plot']["verbose"]):print(_info)
                logging.info(_info)

            # read time of GW
            _jdgw = astropy.time.Time(header['MJD-OBS'], format='mjd') + 2400000.5
            timegw = _jdgw.utc.datetime

            # read distance
            # for offline search, there's no ld available
            try:
                Dmean = header['DISTMEAN']
                Dvar = header['DISTSTD']
                _ld = str(Dmean) + ',' + str(Dvar)    
                _info = 'Distance = %s+/-%s'%(Dmean,Dvar)
            except:
                _ld = None
                _info = 'unmodelled trigger without distance'

            if eval(optlist['arg']['plot']["verbose"]):print(_info)           
            logging.info(_info)

            # for Plotting
            # healpix show
            ''' fignum 1'''
            pparams = {'hpmap':maptrigger,'title':'trigger sky map',\
                       'rot_phi':float(optlist['arg']['plot']["rot_phi"]),\
                       'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'fignum':0,\
                       'ordering':eval(optlist['arg']['plot']["ordering"]),\
                       'coord':optlist['arg']['plot']["coord"],\
                       'norm':str(optlist['arg']['plot']["norm"])}
            optparams = ['rot_theta','rot_phi']

            if pmet==3 and 'trigger' in _showmap:pstplot.interactive_show(pstplot.mollview,pparams,optparams)
            if pmet in [1,2] and 'trigger' in _showmap:
                if all(i for i in _alist.values()):fig_trigger = pstplot.contourview(pparams)
                else:fig_trigger = pstplot.mollview(pparams)
            '''
            if pmet==2 and 'trigger' in _showmap:                
                pparams = {'ra':[255.58],'dec':[-12.4856],'fignum':0,'color':'k','rot_phi':float(optlist['arg']['plot']["rot_phi"]),'rot_theta':float(optlist['arg']['plot']["rot_theta"]),'coord':str(optlist['arg']['plot']["coord"]),'ms':2,'label':'UVOT candidate'}
                pstplot.pointview(pparams)       
                input('trigger map')
            '''

            ## cross map with trigger
            # length = 12*nside**2            
            _pmap = maptrigger
            _info = "%i sec to finish trigger priorization"%int(time.time()-start_time)
            if eval(optlist['arg']['plot']["verbose"]):print(_info)
            logging.info(_info)
    else:grace_id = 'no_astro'    

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

def gwdist(prob, distmu, distsigma, distnorm,ra,dec):
     
    distmin,distmax,delta,frac = 0,10000,1000,.1
    npix = len(prob)                  
    nside = hp.npix2nside(npix)

    rminlist,rmaxlist = [],[]
    for _ra,_dec in zip(ra,dec):
        theta = 0.5 * np.pi - np.deg2rad(_dec)
        phi = np.deg2rad(_ra)
        ipix = hp.ang2pix(nside, theta, phi)
        r = np.linspace(distmin,distmax,num=delta)   
        dp_dr = r**2 * distnorm[ipix] * scipy.stats.norm(\
                    distmu[ipix], distsigma[ipix]).pdf(r)
        rl = r[np.where(dp_dr>frac*max(dp_dr))]        
        if len(rl)>0: 
            rminlist.append(min(rl))
            rmaxlist.append(max(rl))
        else: 
            rminlist.append(None)
            rmaxlist.append(None)
    return rminlist, rmaxlist

def get_skymap(skymap_url,graceid,_dir):
    """
    Look up URL of sky map in VOEvent XML document,
    download sky map, and parse FITS file.
    """  
    import requests,tempfile,shutil

    if False:
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

    else:
        # wget
        os.system(' '.join(['wget',skymap_url,'-O',_dir + graceid+'_'+os.path.basename(skymap_url)]))

    # Done!
    return _dir + graceid +'_'+os.path.basename(skymap_url)

def trigger_validation(skymap, header):

    # validate a trigger by using the header of fits
    ''' judge interests of the trigger '''

    # read and print some values from the FITS header.
    ''' need to be noticed: sometimes failed'''
    header = dict(header)

    # read distance
    try:
        Dmean = header['DISTMEAN']
        Dvar = header['DISTSTD']
        _ld = '%.2f+/-%.2f'%(Dmean,Dvar)
        _str = '#\tDist:\t%s Mpc\n'%_ld
    except:
        _ld = None
        _str = '#\tDist:\tNot available\n'

    # show params of the fits header
    _strremain = ['DATE-OBS', 'CREATOR', 'MJD-OBS']
    for ii in header:
        if ii in _strremain: _str+='#\t%s:\t%s\n'%(ii,header[ii])

    # check the area of skymap
    ilist,hpx = contour(skymap)

    # get area in single pixel
    _areasingle = (hp.nside2resol(hp.get_nside(hpx), arcmin=True)/60.)**2

    # for each contour
    _str += '#\tSky localization:\t'
    _area={}
    for _cc in ilist:
        index1=ilist[_cc]
        if len(index1)>0:
            _str += ' %.2f sq. deg [%i%%]'%(len(index1)*_areasingle,_cc*100)
            _area[_cc]=len(index1)*_areasingle
        else:
            _str += ' NULL [%i%%]'%(_cc*100)
            _area[_cc]=None
    _str+='\n'
    return ilist,_area,_ld,_str,header['DATE-OBS'],header['MJD-OBS']

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
