import os
import time
import re
from slackclient import SlackClient
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Circle, FancyArrowPatch
import matplotlib.patheffects as PathEffects
import pst
from pst import sqlconn,p2api
from astroquery.eso import Eso
from astroquery.vizier import Vizier
from astropy.io import ascii
from astroquery.simbad import Simbad
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, get_sun, get_moon
from astropy.visualization import PercentileInterval, AsinhStretch
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.time import Time,TimeDelta
import datetime
import pandas as pd
import numpy as np
import healpy as hp
import requests
import warnings
import wget
from PIL import Image
from io import BytesIO
import urllib
import pprint
from scipy.stats import norm
import sys
import subprocess
import psutil
import shlex
import math as mt
import sewpy
from scipy.spatial import distance
import logging
logging.basicConfig(filename='slack.log', \
                    level=logging.INFO,\
                    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
                    datefmt='%Y-%m-%d %H:%M:%S')
for key in logging.Logger.manager.loggerDict:logging.getLogger(key).setLevel(logging.CRITICAL)
_pstpath = pst.__path__[0]
print 'pstools version: ',_pstpath

#  define environ
botid = "BOT_ID_grawita"
token = 'SLACK_BOT_TOKEN_grawita'

if False:
    #check bot id
    BOT_NAME = 'gwhelp'
    slack_client = SlackClient(os.environ.get(token))
    api_call = slack_client.api_call("users.list")
    if api_call.get('ok'):
        # retrieve all users so we can find our bot
        users = api_call.get('members')
        for user in users:
            if 'name' in user and user.get('name') == BOT_NAME:
                print("Bot ID for '" + user['name'] + "' is " + user.get('id'))
    else:print("could not find bot user with the name " + BOT_NAME)
    sys.exit()

BOT_ID = os.environ.get(botid)
AT_BOT = "<@" + BOT_ID + ">"
slack_client = SlackClient(os.environ.get(token))

if False:
    # chennel
    channel_list = {'CC9QKJ3PX':'general(OAPD)',\
                    'DJBRPHUD8':'GWalert(OAPD)',\
                    'CJTH4UBKJ':'general(hmt)',\
                    'DJTHHB4PJ':'gwbot(hmt)',\
                    'DJWKK2E4U':'gwhelp(GRAWITA)',\
                    'CJZB7F7DL':'alerts(GRAWITA)'}
    msg1 = "Since running some codes is time consuming, in order to be mroe efficient,"+\
           " I will response only to some specific channel, see below: \n"
    for _schannel in channel_list:msg1 += '\t\t - %s \n'%channel_list[_schannel]

if False:
    # userlist
    _usrl = {}
    _userall = slack_client.api_call("users.list")

    for _usr in _userall['members']:_usrl[_usr['id']] = _usr['name']
    if False: print(_usrl)
    _susers = ['enrico.cappellaro',\
               'sheng.yang']
    msg2 = "Since running some codes is time consuming, in order to be mroe efficient,"+\
           " I will response only to the person in charge of this week, see below: \n"
    for _userweek in _susers:msg2 += '\t\t - %s \n'%_userweek

class handle_command(object):
    """
        Receives commands directed at the bot and determines if they
        are valid commands. If so, then acts on the commands. If not,
        returns back what it needs for clarification.
    """
    def __init__(self,command, channel):
        self.cmdlist = {'general':['help'],\
                        'other tools':['visibility','extinction','gwdist','findingchart'],\
                        'database':['url','check','sql','tmap'],\
                        'query catlog':['GLADE','asteroid','simbad','TNS',\
                                        'ESO','NED','VST'],\
                        'query image':['DSS','Panstarrs','DES'],\
                        'survey tool':['galaxies','pointings','vismap'],\
                        'GW monitor':['status','On','Off']}
        self.cmdlist1 =  [x for xs in self.cmdlist.values() for x in xs]        
        self.command = command
        self.channel = channel
        qresult = sqlconn.query(['select * from ligoevents'])        
        self.tlist=[]
        for item in qresult:self.tlist.append(item['GraceID'])
        self.response,self.attachments,self.fname = "Not sure what you mean. "+\
                            "Use the *help* command to see the defined commands and usage.\n"+\
                            "and then use:`cmd:arg1=val1;arg2=val2...`",False,False

        if self.command.strip() == 'help': 
            self.helper_general()    
        elif self.command.strip() == 'status': 
            self.status()
        elif self.command.strip() == 'On': 
            self.start()
        elif self.command.strip() == 'Off': 
            self.terminate()
        else:
            if len(self.command.split(':')) > 1:                              
                self.cmd = False
                for _ccmd in self.cmdlist1:
                    if self.command.startswith('%s:'%_ccmd): self.cmd = _ccmd               
                if self.cmd:                   
                    self.argv = self.command.split('%s:'%self.cmd)[1]

                    # for cmd with no =
                    if self.cmd == 'help': self.helper()
                    elif self.cmd == 'check': self.check()
                    elif self.cmd == 'url': self.url()
                    elif self.cmd == 'sql': self.sql()
                    elif self.cmd == 'tmap': self.tmap()
                    elif self.cmd == 'VST': self.VST()

                    # for cmd with key=val
                    elif self.cmd in self.cmdlist1:                       
                        if True:
#                        try:
                            # read args
                            self.args = dict(e.strip().split('=') for e in self.argv.split(';'))                           
                            self.read_args()                           

                            # tasks                           
                            if self.cmd == 'galaxies': self.galaxies()
                            if self.cmd == 'vismap': self.vismap()
                            if self.cmd == 'visibility': self.visibility()
                            if self.cmd == 'gwdist': self.gwdist()                            
                            if self.cmd == 'GLADE':self.GLADE()
                            if self.cmd == 'extinction':self.extinction()
                            if self.cmd == 'asteroid':self.asteroid()
                            if self.cmd == 'simbad':self.simbad()
                            if self.cmd == 'TNS':self.TNS()
                            if self.cmd == 'DSS':self.DSS()
                            if self.cmd == 'Panstarrs':self.Panstarrs()
                            if self.cmd == 'DES':self.DES()
                            if self.cmd == 'findingchart':self.findingchart()
                            if self.cmd == 'ESO':self.ESO()
                            if self.cmd == 'NED':self.NED()
                            if self.cmd == 'Skymappercat':self.Skymappercat()
                            if self.cmd == 'Skymapperimage':self.Skymapperimage()
#                        except: self.response = 'check parameters for `%s`'%self.cmd                    
                    else: self.response = 'command %s not defined!'%self.command.split(':')[0]
                else: self.response = 'command %s not defined!'%self.command.split(':')[0]

        # output
        if self.attachments:
            slack_client.api_call("chat.postMessage", channel=channel, \
                                  text=self.response, attachments=self.attachments)
        elif self.fname:
            with open(self.fname, 'rb') as file_content:
                slack_client.api_call(
                    "files.upload",
                    channels=channel,
                    file=file_content,
                    title=self.response,
                    username='mybot')
        else:
            slack_client.api_call("chat.postMessage", channel=channel,
                                  text=self.response, as_user=True)        

    #########################
    def read_args(self):

        try:self.ra = str(self.args['ra'])                         
        except:self.ra = False            
       
        try:self.dec = str(self.args['dec'])
        except:self.dec = False            

        try:self.ral = str(self.args['ral'])                         
        except:self.ral = False            
       
        try:self.decl = str(self.args['decl'])
        except:self.decl = False   

        try:self.magl = str(self.args['magl'])                         
        except:self.magl = False            
       
        try:self.distl = str(self.args['distl'])
        except:self.distl = False 

        try:self.name = str(self.args['name'])
        except:self.name = False           

        try:self.lon = str(self.args['lon'])                         
        except:self.lon = False           

        try:self.lat = str(self.args['lat'])
        except:self.lat = False            

        try:self.alt = str(self.args['alt'])
        except:self.alt = False           

        try:self.url = str(self.args['url'])
        except:self.url = False           

        try:self.tdiff = float(self.args['timediff'])
        except:self.tdiff = False           

        try:self.time = str(self.args['time'])
        except:self.time = False            

        try:self.radius = float(self.args['radius'])
        except:self.radius = False        

        try:self.mag = float(self.args['mag'])
        except:self.mag = False            

        try:self.instr = str(self.args['instrument'])
        except:self.instr = False           

        try:self.number = float(self.args['number'])
        except:self.number = False

        try:self.area = float(self.args['area'])
        except:self.area = False

        try:self.img = self.args['image']
        except:self.img = None

        if self.ra and self.dec:
            self.radeg, self.decdeg, self.rahms, self.decdms = readradec(self.ra,self.dec)

    def helper_general(self):
        self.response = '\n*command list*:\n'                        
        for _ccmdk in self.cmdlist:
            self.response += ' - *%s*: \t'%_ccmdk
            for _ccmd in self.cmdlist[_ccmdk]:
                self.response += '`%s`\t'%_ccmd
            self.response += '\n\n'
        self.response += "\n\n"
        self.response += "\nuse `help:command` to see more detailed tutorials for specific command\n"
        self.response += "\nuse `help:all` to see all\n"
        
    def verification(self):
        self.pid = False
        for pid in psutil.pids():
            try:p = psutil.Process(pid)
            except:continue
            if p.name() == "python" and len(p.cmdline()) > 1 and "pstools.py" in p.cmdline()[1]:                                
                self.pid = pid

    def status(self):
        self.verification()
        if self.pid:self.response = 'GW monitor status: *On*'          
        else: self.response = 'GW monitor status: *Off*'        

    def start(self):
        self.verification()        
        if self.pid:
            self.response = 'Skip\nGW monitor already *On*'
        else:
            if True:
                proc = subprocess.Popen(['python', '%s/../../bin/pstools.py'%_pstpath, '-g'])
            else:
                proc = subprocess.Popen(['python', '%s/../../bin/pstools.py'%_pstpath, '-g','--server','local'])
            self.response = 'Start\nGW monitor *On*'

    def terminate(self):
        self.verification()        
        if self.pid:
            process = psutil.Process(self.pid)
            try:
                process.terminate()  #or p.kill()
                self.response = 'Closed\nGW monitor status *Off*'    
            except:
                self.response = 'try again'    
        else:self.response = 'GW monitor status already *Off*'

    def helper(self):       
        radecformat = "\t`ra/dec/radius format: deg or hms`\n"
        radecformat += "\t\t`cmd: ra=00h42.5m;dec=+41d12m`\n"
        radecformat += "\t\t`cmd: ra=00:42:30;dec=41:12:00`\n"
        radecformat += "\t\t`cmd: ra=00 42 30;dec=41 12 00`\n"
        radecformat += "\t\t`cmd: ra=42.00;dec=41.00`\n"

        self.response = ''        
        if self.argv in self.cmdlist1 or self.argv=='all':
            if self.argv == 'all':self.argv = self.cmdlist1
            if 'help' in self.argv:
                self.response += "*help: <command>*: \tshow helper of input command\n"
            if 'status' in self.argv:
                self.response += "*status*: \tcheck status of GW listener\n"
            if 'On' in self.argv:
                self.response += "*On*: \tTurn on GW listener\n"
            if 'Off' in self.argv:
                self.response += "*Off*: \tTurn off GW lintener\n"
            if 'findingchart' in self.argv:
                self.response += "*findingchart: ra=<ra> ; dec=<dec> ; image=<fits, default=DSS>*:"
                self.response += "\tcreate finding chart for specific positions; "
                self.response += "if image not given, use DSS \n\te.g. `findingchart: ra=10;dec=-10`\n"
            if 'VST' in self.argv:
                self.response += "*VST: <GraceID>*:"
                self.response += "\tcheck VST accomplished OBs from P2 page, search "
                self.response += "also if they're observed or not \n\te.g. `VST: S190510g`\n"
            if 'check' in self.argv:
                self.response += "*check: <GraceID>*:"
                self.response += "\tcheck trigger from our database \n\te.g. `check: S190510g`\n"
            if 'sql' in self.argv:
                self.response += "*sql:<mysql commadn>*: \texcute sql in our GW database\n"
                self.response += "\te.g. `sql: select * from ligoevents limit 10`\n"
                self.response += "\t\t`sql: select * from ligoevents where Dmean<100`\n"
            if 'url' in self.argv:
                self.response += "*url:<GraceID>*: \turl of our GW page on the trigger\n"
                self.response += "\te.g. `url: S190510g`\n"
            if 'tmap' in self.argv:
                self.response += "*tmap:<GraceID>*: \t2D localization map of the trigger\n"
                self.response += "\te.g. `tmap: S190510g`\n"
            if 'galaxies' in self.argv:
                self.response += "*galaxies:url=<fits url>;number=<num>;number=<area, choices=.5/.68/.9/.99>*: \tselect a number of galaxies\n"
                self.response += "\te.g. `galaxies: number=50;area=.68;url=https://gracedb.ligo.org/"
                self.response += "api/superevents/S190408an/files/bayestar.fits.gz`\n"
            if 'gwdist' in self.argv:
                self.response += "*gwdist:ra=<ra> ; dec=<dec> ; distl=<dist range: [mpc,mpc]> ; url=<fits url>*: "
                self.response += "\tGW distance estimation for specific direction\n"
                if self.argv=='gwdist':self.response += radecformat
                self.response += "\te.g. `gwdist:ra=00:42:30;dec=41:12:00;distl=[2000,4000];"
                self.response += "url=https://gracedb.ligo.org/api/superevents/S190408an/files/bayestar.fits.gz`\n"
            if 'vismap' in self.argv:
                self.response +=  "*vismap:lon=<lon, deg> ; lat=<lat, deg> ; alt=<alt, m> ; "
                self.response += "timediff=<tdiff: time difference from now, "
                self.response += "\te.g. +10, 10 hours later> ; url=<fits url>*: \t"
                self.response += "read fits from url link, show GW 50,68,90,99% contour, "
                self.response += "and the 2D visible region for specific observatory, "
                self.response += "considering constrains of sun, moon, airmass, etc\n"
                self.response += "\te.g. `vismap: lat=-24.625;lon=-70.4033;alt=2635;"
                self.response += "time=+10;url=https://gracedb.ligo.org/api/"
                self.response += "superevents/S190408an/files/bayestar.fits.gz`\n"
            if 'GLADE' in self.argv:
                self.response += "*GLADE: ral= <ra range: [deg,deg]>;"
                self.response += "decl=<dec range: [deg,deg]> ; "
                self.response += "magl=<abs mag range: [mag,mag]>; "
                self.response += "distl=<dist range: [mpc,mpc]>*: \tcheck GLADE galaxies\n"
                self.response += "\te.g. `GLADE: ral=[12,13];decl=[20,29];distl=[0,60];magl=[-18,-22]`\n"
            if 'visibility' in self.argv:
                self.response += "*visibility: ra=<ra> ; dec=<dec> ; name=<name, default=GWC1>*:"
                self.response += " \tcheck visibility\n"
                if self.argv=='visibility':self.response += radecformat
            if 'extinction' in self.argv:
                self.response += "*extinction: ra=<ra> ; dec=<dec>*: \t check extinction\n"           
                if self.argv=='extinction':self.response += radecformat
            if 'asteroid' in self.argv:
                self.response += "*asteroid: ra=<ra> ; dec=<dec> ; time=<time, 'yyyy-mm-dd.dd'> ; "
                self.response += "radius=<rad, arcsec> ; mag=<mag cut, V>*: \tcheck minor plantes\n"
                if self.argv=='asteroid':self.response += radecformat
                self.response += "\te.g. `asteroid:ra=17 58 52.23;dec=+57 55 12.4;time=2019-05-05.41;radius=10;mag=0`\n"
            if 'simbad' in self.argv:
                self.response += "*simbad: ra=<ra> ; dec=<dec> ; radius=<rad, arcsec>*: "
                self.response += "\tcheck Simbad catalog\n"
                if self.argv=='simbad':self.response += radecformat
                self.response += "\te.g. `simbad:ra=00h42.5m;dec=+41d12m;radius=50`\n"
            if 'TNS' in self.argv:
                self.response += "*TNS: ra=<ra> ; dec=<dec> ; radius=<rad, arcsec> ; "
                self.response += "name=<name, default=GWT1>*: \tcheck TNS catalog\n"
                if self.argv=='TNS':self.response += radecformat
                self.response += "\te.g. `TNS:ra=07:55:00.889;dec=-76:24:43.11;radius=10`\n"
            if 'DSS' in self.argv:
                self.response += "*DSS: ra=<ra> ; dec=<dec> ; name=<name, default=GWT1>*: \t"
                self.response += "show DSS cutout image \n"
                if self.argv=='DSS':self.response += radecformat
                self.response += "\te.g. `DSS:ra=07:55:00.889;dec=-76:24:43.11;name=AT19xyz`\n"
            if 'Panstarrs' in self.argv:
                self.response += "*Panstarrs: ra=<ra> ; dec=<dec> ; name=<name, default=GWT1>*: \t"
                self.response += "show Panstarrs cutout image \n"
                if self.argv=='Panstarrs':self.response += radecformat
                self.response += "\te.g. `Panstarrs:ra=00h42.5m;dec=+41d12m;name=AT19xyz`\n"
            if 'DES' in self.argv:
                self.response += "*DES: ra=<ra> ; dec=<dec>*: \t"
                self.response += "check if available DES reference images\n"
                if self.argv=='DES':self.response += radecformat
            if 'ESO' in self.argv:
                self.response += "*ESO: ra=<ra> ; dec=<dec> ; radius=<rad, arcsec> ; "
                self.response += "instrument=<instr>*: \tcheck ESO archive \n"
                if self.argv=='ESO':self.response += radecformat
            if 'NED' in self.argv:
                self.response += "*NED: ra=<ra> ; dec=<dec> ; radius=<rad, arcsec>*: \t"
                self.response += "check NED galaxies \n"
                if self.argv=='NED':self.response += radecformat
            if 'Skymappercat' in self.argv:
                self.response += "*Skymappercat: ra=<ra> ; dec=<dec> ; radius=<rad, arcsec>*: \t"
                self.response += "check Skymapper catalog \n"
                if self.argv=='Skymappercat':self.response += radecformat
            if 'Skymapperimage' in self.argv:
                self.response += "*Skymapperimage: ra=<ra> ; dec=<dec> ; radius=<rad, arcsec>*: \t"
                self.response += "check Skymapper images \n"
                if self.argv=='Skymapperimage':self.response += radecformat
            if 'contact' in self.argv:
                self.response += "*contact: <message>*: \tWrite to me if you found sth wrong"
            if type(self.argv) is list:
                self.response += 'Ra Dec format:\n\n'
                self.response += radecformat                
        else:
            self.response = '%s not found in helper'%self.argv

    def check(self):     
        self.response = 'trigger %s not found in our DB'%self.argv
        for trigger in self.tlist:           
            ck = re.search(trigger,self.command,re.IGNORECASE)
            if ck:
                self.response ='Results: \n'               
                _command = 'select * from ligoevents'+\
                           ' where GraceID="%s"'%ck.group(0)                                
                _result = sqlconn.query([_command])                                                
                _nn = 0
                for _rr in _result:
                    _nn+=1
                    self.response+=str(_nn)
                    for _key in _rr:
                        self.response += ' %s:%s \t'%(_key,_rr[_key])
                    self.response+='\n\n'
    def galaxies(self):
        self.response=''
        outf = '/tmp/tmphp.fits'
        if self.url:
            if not self.number:
                self.number = 50               
                self.response+='input num wrong, use 50\n' 
            if not self.area:
                self.area = .68
                self.response+='input area wrong, use .68\n' 

            try:
                wget.download(self.url,out=outf)
                if not os.path.exists(outf):
                    self.response += 'wrong with url\n'
                else:
                    try:
                        (maptrigger, distmu, distsigma, distnorm), header = \
                                        hp.read_map(outf, field=[0, 1, 2, 3], h=True)
                        header = dict(header)
                        # read distance
                        try:
                            Dmean = header['DISTMEAN']
                            dsigma = header['DISTSTD']
                            self.response += 'dist in fits header: %.2f \n'%Dmean
                            npix = len(maptrigger)                  
                            nside = hp.npix2nside(npix)
                            _ilist,_alist,_ld,_str,_date_obs,_mjd_obs = pst.trigger_validation(maptrigger, header)
                        
                            _id0,_gname0,_ra0,_dec0,_mag0,_dist0,fig_g1,fig_g3,fig_g4 = \
                            pst.scheme.galaxies(catname="GLADE",\
                            limra=[0,360],limdec=[-90,90],interactive=False,\
                            distmin=0,distmax=2000,absmag=-18,outfits=False,\
                            size=-1,mcolor='k',mcolor2='r',mcolor3='grey',\
                            gcolor='grey',gcolor2='g-',gcolor3='r-',_dir='./',\
                            rot_theta=0,rot_phi=0,nside=nside,ordering=False,\
                            coord='C',_showmap='trigger,cum',norm = 'hist',\
                            verbose = True,pmet=1,cachefile='tmp_glade.npz',\
                            ypos1=2,ypos2=-2,nbindist=2000,nbinlums=1,figset=[2,3,4,5])

                            # ra,dec cut
                            _indexg = pst.DeclRaToIndex(_dec0,_ra0,nside)
                            _indexgin = np.in1d(_indexg, _ilist[self.area])
                            _id0,_gname0,_ra0,_dec0,_mag0,_dist0 = _id0[_indexgin],_gname0[_indexgin],\
                                                                   _ra0[_indexgin],_dec0[_indexgin],\
                                                                   _mag0[_indexgin],_dist0[_indexgin]
                            # dist cut
                            rminlist,rmaxlist = pst.gwdist(maptrigger, distmu, distsigma, distnorm,_ra0,_dec0)
                            _iddc,_distscore = [],[]
                            for _ii in range(len(_id0)):
                                _dmin,_dmax,_dgal = rminlist[_ii],rmaxlist[_ii],_dist0[_ii]                                  
                                if _dmin is None or _dmax is None:continue
                                if np.logical_and(_dgal>_dmin,_dgal<_dmax):
                                    dmean = (_dmin+_dmax)/2.
                                    _iddc.append(_ii)
                                    _distscore.append(np.e**(-(_dgal-dmean)**2/2./dsigma**2))
                            _iddc = np.array(_iddc)
                            
                            if not _iddc is None:
                                _id0,_gname0,_ra0,_dec0,_mag0,_dist0 = _id0[_iddc],_gname0[_iddc],_ra0[_iddc],_dec0[_iddc],_mag0[_iddc],_dist0[_iddc]
                                self.response += '%i galaxies found depends on settings\n'%len(_id0)
                                self.response += 'Index Ra Dec Bmag Dist Score \n'

                                # 2d priorization
                                gradius = 5*2*mt.pi/60/360
                                _idc,_rac,_decc,_score = pst.priorization.calprob_gal(maptrigger,_ra0,_dec0,_id0,radius=gradius)    

                                # 3d + mass                                
                                _score = _score*_distscore*10**((-1)*(_mag0/2.5))

                                # score normalization
                                _score = _score/sum(_score)    

                                # sort       
                                idx=np.argsort(np.asarray(_score))[::-1]
                                _rac,_decc,_idc,_score=_rac[idx],_decc[idx],_idc[idx],_score[idx]                

                                # output
                                _oo = 0
                                for _ra00,_dec00,ii,_score00 in zip(_rac,_decc,_idc,_score):
                                    if _oo>=self.number:continue
                                    _oo+=1                                    
                                    jj = np.where(_id0 == ii)
                                    _id00,_gname00,_mag00,_dist00 = _id0[jj],_gname0[jj],_mag0[jj],_dist0[jj]
                                    self.response += '%i \t %s \t %.2f \t %.2f \t %.2f \t %.2f \t %.2e \n'%(_oo,_gname00[0],_ra00[0],_dec00[0],_mag00[0],_dist00[0],_score00)
                            else:
                                self.response += 'No galaxies found after dist,ra,dec,mag cut'
                        except Exception as e:self.response += str(e)                      
                    except: 
                        self.response+='wrong fits format for healpy!!'
            except Exception as e:self.response += str(e)
        else:self.response += 'url is needed!'
 

    def VST(self):
        self.response = 'trigger %s not found in our DB'%self.argv
        # acontainer for ToO, b for follow-up
        enviroment, username, password, acontainerid, bcontainderid, instr = \
                        'production','saberyoung','Ys19900615',2229852,2229505, 'omegacam'
        _size = '01:00:00' #VST size
        _stime,_etime = '2000-01-01','2020-01-01'       
        for trigger in self.tlist:
            ck = re.search(trigger,self.command,re.IGNORECASE)
            if ck:
                self.response ='Results: \n'
                trigger = ck.group(0)

                # connect to astropy.eso
                eso = Eso()

                # need password at the first time
                eso.login(username, store_password=True)               

                # connect to p2 page
                p = pprint.PrettyPrinter(indent=4)

                # login
                api = p2api.ApiConnection(enviroment,username,password)                

                # check items via p2
                _coblist = {}
                for nn,_cid in zip(['ToO','follow'],[acontainerid, bcontainderid]):                        
                    items, itemsVersion = api.getItems(_cid)
                    for _item in items:
                        # First level: Folder
                        if _item['itemType']=='Folder':
                            _trigger = _item['name']
                            if not trigger in _trigger:continue
                            _coblist[_trigger] = {}
                            subitems, subitemsVersion = api.getItems(_item['containerId'])
                            # second level: Concatenation
                            for _subitem in subitems:               
                                if _subitem['itemType']=='Concatenation':                    
                                    # third level: OB
                                    obs, obsVersion = api.getItems(_subitem['containerId'])
                                    for _ob in obs:
                                        # get OB
                                        ob, obVersion = api.getOB(_ob['obId'])                       
                                        if ob['obStatus'] == 'C':
                                            # if OB done                                
                                            # check eso archive
                                            table = eso.query_instrument(instr, \
                                                    column_filters={'ob_id':_ob['obId']}, \
                                                    columns=['prog_type'])                               
                                            _tlist,flist,slist=[],[],[]
                                            for _table in table:
                                                info,filtro,fwhm = _table['DP.ID'].split('.')[1],\
                                                                   _table['INS.FILT1.NAME'],\
                                                                   _table['DIMM Seeing-avg'].split()[0]
                                                _tlist.append(Time(info))
                                                flist.append(filtro)
                                                slist.append(float(fwhm))

                                            t = (_tlist[0]-_tlist[1])/2+_tlist[0]
                                            t.format = 'datetime'
                                            _mjd = t.mjd
                                            _y,_m,_d,_h = t.value.year,t.value.month,t.value.day,t.value.hour
                                            if _h>12:_date = '%s-%.2i-%.2i'%(_y,_m,_d)
                                            else:_date = '%s-%.2i-%.2i'%(_y,_m,_d-1)

                                            _filtro = np.unique(flist)[0]
                                            _seeing = np.mean(slist)
                             
                                            _coblist[_trigger][ob['name']] = \
                                                    [ob['target']['ra'], ob['target']['dec'], \
                                                    _mjd, _date, _filtro, _seeing]

                    for epoch in _coblist:
                        if epoch == trigger:self.response+='epoch 1\n'
                        else:self.response+=('epoch %s\n'%re.search('_\d{1,2}',epoch).group(0).split('_')[1])
                        for _name in _coblist[epoch]:
                            _ra,_dec, _mjd, _date, _filtro, _seeing = _coblist[epoch][_name]
                            pointing = re.search('p\d{1,3}',_name).group(0).split('p')[1]
                            _p = SkyCoord('%s %s'%(_ra,_dec),unit=(u.hourangle, u.deg))                                 
                            _rad,_decd = _p.ra.deg,_p.dec.deg
                            pra,pdec = str('%.1f'%_rad).replace('.','').replace('-',''),\
                                       str('%.1f'%_decd).replace('.','').replace('-','')
                            self.response += '%s %s %s(%s,%s) %s %s\n'%(trigger,_filtro,pointing,pra,pdec,_date,_mjd)

    def sql(self):              
        _result = sqlconn.query([self.argv])                   
        if _result:
            self.response ='Results: \n'  
            _nn = 0
            for _rr in _result:
                _nn+=1
                self.response+=str(_nn)
                for _key in _rr:
                    self.response += ' %s:%s \t'%(_key,_rr[_key])
                self.response+='\n\n'
        else:self.response ='No records found'

    def url(self):
        self.response = '%s not trigger found in our DB'%self.argv
        for trigger in self.tlist:
            ck = re.search(trigger,self.command,re.IGNORECASE)
            if ck:self.response = 'http://sngroup.oapd.inaf.it/gw/pic/%s/'%ck.group(0)   

    def findingchart(self):
        self.response=''
        self.fov = 5.
        self.fov_inset = 1.
        self.rt1,self.rt21,self.rt22,self.nlim = 2,2,200,4

        # sex params
        self.threshold = 3
        self.seeing = 5
        self.pixel_scale = 1
        self.fwhm = float(self.seeing)/self.pixel_scale
        self.phot_aperture = '10,15,20'
        self.analysis_thresh = 3
        self.detect_minarea = 2
        self.mag_zeropoint = 30
        self.satur_level = 30000
        self.back_size = 64
        if self.ra and self.dec:
            radeg, decdeg, rahms, decdms = readradec(self.ra,self.dec)
            radeg, decdeg = float(radeg), float(decdeg)
            # size
            fov = self.fov/60 # convert to degrees
            fovinset = self.fov_inset/60
            # image
            if self.img is None:
                self.response += 'use DSS image as the finding chart\n'
                img = download_dss(rahms,decdms,fov*60,fov*60,True,True)
            if img is None: 
                self.response+='wrong image'           
            else:
                # zscale           
                _z1,_z2 = zscale(fits.getdata(img))
                _z1inset = _z1
                _z2inset = _z2

                # wcs
                X, hdr = fits.getdata(img, header=True)
                height, width = X.shape
                if not self.name:
                    if 'object' in hdr: self.name = hdr['object']
                    else:self.name = 'GWT1'
                left = col2RA(0.5,hdr)
                right = col2RA(width+0.5,hdr)
                bottom = row2dec(0.5,hdr)
                top = row2dec(height+0.5,hdr)

                # out region?
                if min(left,right)<=radeg<=max(left,right) \
                   and min(bottom,top)<=decdeg<=max(bottom,top):                    

                    # select the object
                    sextab = self.runsex(img)
                    _coord=[]
                    for _ra,_dec in zip(sextab['X_WORLD'],sextab['Y_WORLD']):
                        _coord.append([_ra,_dec])

                    # distance                    
                    _dist = distance.cdist([[radeg,decdeg]], _coord, 'euclidean')[0]*60*60
                    _index1 = np.where(_dist<=self.rt1)    
                    if len(_index1)>0:
                        self.response += 'target found'
                        self.rt21 = _dist[np.argsort(_dist)][0]        
                    else:self.response += 'no target'
                    _index2 = np.logical_and(_dist>self.rt21,_dist<self.rt22)
                    if len(_index2)>0:
                        _rar,_decr,_distr = sextab['X_WORLD'][_index2],\
                                            sextab['Y_WORLD'][_index2],\
                                            _dist[_index2]
                        _rar,_decr,_distr = _rar[np.argsort(_distr)][:self.nlim],\
                                            _decr[np.argsort(_distr)][:self.nlim],\
                                            _distr[np.argsort(_distr)][:self.nlim]       
                    else:self.response += 'no rference star'

                    # plot
                    # target
                    self.fname = '/tmp/tmpFC.png'
                    plt.rcParams['figure.figsize'] =12, 7
                    plt.rcParams['figure.facecolor']='white'
                    ax = plt.axes([.09,.103,.5,.85])
                    ax.imshow(X, cmap='gray_r', aspect='equal', interpolation='nearest', vmin=_z1, vmax=_z2, origin='lower', extent=(left, right, bottom, top))
                    plt.xlabel('R.A. (degrees)')
                    plt.ylabel('Dec. (degrees)')
                    ax.get_xaxis().get_major_formatter().set_useOffset(False)
                    ax.get_yaxis().get_major_formatter().set_useOffset(False)

                    plt.plot(radeg, decdeg, marker='o', mfc='none', ms=25, mew=2, mec='r', label='SN '+rahms+' '+decdms, linestyle='none')
                    plt.xlim(radeg + fov/2, radeg - fov/2)
                    plt.ylim(decdeg - fov/2, decdeg + fov/2)

                    # draw the compass
                    xcomp = radeg - 0.9*fov/2
                    ycomp = decdeg - 0.9*fov/2
                    plt.arrow(xcomp, ycomp, 0, 1/60., fc='k', ec='none', length_includes_head=True)
                    plt.arrow(xcomp, ycomp, 1/60., 0, fc='k', ec='none', length_includes_head=True)
                    plt.text(xcomp, ycomp+1/60., 'N', ha='center', va='bottom', size='large', weight='bold')
                    plt.text(xcomp+1/60., ycomp, 'E', ha='right', va='center', size='large', weight='bold')
                    plt.text(xcomp+1/120., ycomp+0.001, "1'", ha='center', va='bottom', size='large', weight='bold')

                    # find and draw the reference stars  
                    simbolo=['d','*','p','s']
                    colore=['b','g','m','c','w']

                    for nn,(_ra,_dec) in enumerate(zip(_rar,_decr)):                       
                        _radeg, _decdeg, _rahms, _decdms = readradec(str(_ra),str(_dec))               
                        plt.plot(float(_radeg), float(_decdeg), marker=simbolo[nn], mfc='none', ms=25, mew=2, mec=colore[nn], label='#'+str(nn+1)+' '+_rahms+' '+_decdms, linestyle='none')
        
                        # offset
                        mm1 = SkyCoord(_rahms, _decdms,unit=(u.hourangle,u.deg))
                        mm2 = SkyCoord(rahms, decdms,unit=(u.hourangle,u.deg))        
                        ra_offset = (mm2.ra - mm1.ra) * np.cos(mm1.dec.to('radian'))
                        dec_offset = (mm2.dec - mm1.dec)
                    
                        AB = ra_offset.to('arcsec').value
                        CD = dec_offset.to('arcsec').value                    
                        if float(AB)<0:  AB1='W'
                        else:            AB1='E'
                        if float(CD)<0:  CD1='S'
                        else:            CD1='N'
                        plt.figtext(0.7, 0.7-nn*.06, 'From #'+str(nn+1)+' to SN', color=colore[nn], backgroundcolor='white', size='large')
                        plt.figtext(0.7, 0.67-nn*.06, '{:4.2f}"{} {:4.2f}"{}'.format(abs(AB), str(AB1), abs(CD), str(CD1)), color=colore[nn], backgroundcolor='white', size='large')
                        plt.legend(numpoints=1, markerscale=.6, loc=(1.1,.8), fancybox=True)

                    # draw the inset
                    ax = plt.axes([.62,.1,.35,.35])
                    if _z1inset is None: _z1inset = _z1
                    if _z2inset is None: _z2inset = _z2
                    ax.imshow(X, cmap='gray_r', aspect='equal', interpolation='nearest', vmin=_z1inset, vmax=_z2inset, origin='lower', extent=(left, right, bottom, top))
                    plt.xlabel('R.A. (degrees)')
                    plt.xlabel('Dec. (degrees)')
                    plt.title(self.name)
                    ax.get_xaxis().get_major_formatter().set_useOffset(False)
                    ax.get_yaxis().get_major_formatter().set_useOffset(False)
                    ax.plot(radeg, decdeg, marker='o', mfc='none', ms=25, mew=2, mec='r')
                    plt.xlim(radeg + fovinset/2, radeg - fovinset/2)
                    plt.ylim(decdeg - fovinset/2, decdeg + fovinset/2)
                    ax.get_xaxis().set_major_locator(ticker.MultipleLocator(round(fovinset/3,3)))
                    plt.savefig(self.fname,dpi=120,format='PNG')
                    plt.close('all')
                else:self.response += 'ra dec out of region'

    def gwdist(self):

        self.response=''
        outf = '/tmp/tmphp.fits'
        if os.path.exists(outf):os.remove(outf)
        if self.ra and self.dec and self.url:
            if not self.distl:
                self.distl = [0,1000]
                self.response+='no distl fround ,use %s\n'%str(self.distl)
            try:distmin,distmax=min(eval(self.distl)),max(eval(self.distl))               
            except:
                distmin,distmax=0,1000
                self.response+='input dist format wrong, use 0,1000\n' 
            try:
                wget.download(self.url,out=outf)
                if not os.path.exists(outf):
                    self.response += 'wrong with url'
                else:
                    self.fname = '/tmp/tmpvismap.png'
                    self.response += 'GW distance distribution along the sight\n'
                    if os.path.exists(self.fname):os.remove(self.fname)

                    prob, distmu, distsigma, distnorm = \
                                hp.read_map(outf, field=[0, 1, 2, 3])
                    npix = len(prob)                  
                    nside = hp.npix2nside(npix)
                    theta = 0.5 * np.pi - np.deg2rad(float(self.decdeg))
                    phi = np.deg2rad(float(self.radeg))
                    ipix = hp.ang2pix(nside, theta, phi)

                    r = np.linspace(distmin,distmax)                              
                    dp_dr = r**2 * distnorm[ipix] * norm(\
                            distmu[ipix], distsigma[ipix]).pdf(r)

                    # make plot
                    plt.plot(r, dp_dr)
                    plt.xlabel('distance (Mpc)')
                    plt.ylabel('prob Mpc$^{-1}$')
                    plt.title('GW distance distribution along the sight')                                                 
                    plt.savefig(self.fname,dpi=120,format='PNG')
                    plt.close('all')
            except Exception as e:self.response += str(e)
        else:self.response += 'ra and dec and url are needed!'

    def vismap(self):

        self.response=''
        if self.tdiff:
            timenow = Time.now() + TimeDelta(self.tdiff*3600, format='sec') 
        else:
            timenow = Time.now()
            self.response += 'No timediff given, show current time\n'
        _sun,_airmass = -18, 2.0
        outf = '/tmp/tmphp.fits'
        if os.path.exists(outf):os.remove(outf)
        if self.lon and self.lat and self.alt and self.url:
            
            # GW healpix plot
            try:
                wget.download(self.url,out=outf)
                if not os.path.exists(outf):
                    self.response += 'wrong with url'
                else:                    
                    self.fname = '/tmp/tmpvismap.png'
                    self.response += 'visibility map'
                    if os.path.exists(self.fname):os.remove(self.fname)

                    ######## GW map                   
                    fig = plt.figure(figsize=(15, 10))
                    hp.graticule()                   
                    _ilist,hpx = pst.contour(outf)  
                    nside = hp.get_nside(hpx)
                    # get area in single pixel
                    _areasingle = (hp.nside2resol(nside, arcmin=True)/60.)**2
                    label1,label2,label3,label4=_ilist.keys()
                    index1,index2,index3,index4=_ilist.values()                      
                    color1='k'
                    color2='y'
                    color3='r'
                    color4='g'
                    _label='GW loc'
                    locl={}
                    for _index,_color,_cont in zip([index4,index3,index2,index1],\
                                                   [color1,color2,color3,color4],\
                                                   [label4,label3,label2,label1]):
                        if len(_index)==0:continue
                        locl[_cont] = len(_index)*_areasingle                                                         
                        theta,phi = hp.pix2ang(nside,_index)                       
                        try:                           
                            hp.projplot(theta[0],phi[0],_color,coord='C',\
                                        label='%s %s'%(_label,_cont))
                            hp.projplot(theta,phi,_color,coord='C')
                        except:                            
                            hp.projplot(theta,phi,_color,coord='C',\
                                        label='%s %s'%(_label,_cont))
                  
                    ########### plot coord
                    _ralist,_declist = [0,45,90,135,180,225,270,315],\
                                       [0,30,60,-30,-60]
                    for _ra in _ralist:
                        for _dec in _declist:
                            # select some special points
                            theta1,phi1 = pst.RadecToThetaphi(_ra,_dec) 
  
                            # visualization
                            if _ra == 0:
                                hp.projtext(theta1,phi1, str(_dec), coord='C')
                            elif _dec == 0:
                                hp.projtext(theta1,phi1, str(_ra), coord='C')
                            else:
                                hp.projtext(theta1,phi1, str(_ra)+','+str(_dec), coord='C')
                   
                    ###########plot the sky
                    # define telescope                   
                    observatory = EarthLocation(lat=float(self.lat)*u.deg, \
                                                lon=float(self.lon)*u.deg, \
                                                height=float(self.alt)*u.m)

                    # plot the horizon           
                    _smlabel=True
                    for _timeplus in np.arange(0,360,1):
                        newAltAzcoordiantes = SkyCoord(alt = 0*u.deg, \
                                        az = 0*u.deg + _timeplus*u.deg, \
                                        obstime = timenow, \
                                        frame = 'altaz', \
                                        location = observatory)
                       
                        # transform to theta phi
                        _htheta,_hphi = pst.RadecToThetaphi(newAltAzcoordiantes.icrs.ra.deg, \
                                                            newAltAzcoordiantes.icrs.dec.deg)

                        if _smlabel:
                            hp.projplot(_htheta,_hphi,'.', color = 'b', \
                                        coord='C', ms = 2, label='horizon')
                            _smlabel=False
                        else:
                            hp.projplot(_htheta,_hphi,'.', color = 'b', \
                                        coord='C', ms = 2)                   
                  
                    # plot the galactic plane        
                    _hral = np.arange(0,360,10)
                    _hdecl = np.zeros(len(_hral))
                    _smlabel=True
                    for _hra,_hdec in zip(_hral,_hdecl):

                        # from galactic coordinates to equatorial
                        _hradecs = SkyCoord(l=_hra*u.deg, b=_hdec*u.deg, frame='galactic')                

                        # transform to theta phi
                        _htheta,_hphi = pst.RadecToThetaphi(_hradecs.icrs.ra.deg, _hradecs.icrs.dec.deg)
            
                        if _smlabel:
                            hp.projplot(_htheta,_hphi,'x', color = 'k', coord='C', ms = 10, label='galactic plane')
                            _smlabel=False
                        else:hp.projplot(_htheta,_hphi,'x', color = 'k', coord='C', ms = 10)

                    # plot the sun, moon
                    _smlabel=True                    
                    _sra,_sdec = get_sun(timenow).ra.deg,get_sun(timenow).dec.deg
                    _stheta,_sphi = pst.RadecToThetaphi(_sra,_sdec)              
                    _mra,_mdec = get_moon(timenow).ra.deg,get_moon(timenow).dec.deg
                    _mtheta,_mphi = pst.RadecToThetaphi(_mra,_mdec)                                
                    hp.projplot(_stheta,_sphi,'o', fillstyle='left', color='y', coord='C', ms = 30)
                    hp.projplot(_mtheta,_mphi,'o', fillstyle='right', color='b', coord='C', ms = 30) 
                    if _smlabel:
                        hp.projtext(_stheta,_sphi,'sun', size=20, color = 'k', coord='C')
                        hp.projtext(_mtheta,_mphi,'moon', size=20, color = 'k', coord='C')
                        _smlabel=False
                    
                    # title
                    _title = 'trigger localization in equatorial system\n'
                    for ii in locl:_title += '%.2f sky=%.2f deg, '%(ii,locl[ii])
                    _title+='\nObservatory: lon=%s, lat=%s, alt=%s \n'%(self.lon,self.lat,self.alt)
                    _title+='time=%.2i, sun<%.2f, airmass<%.2f\n'%(self.tdiff,_sun,_airmass)

                    # plot the available region by considering airmass, sun, moon                                      
                    theta, phi = hp.pix2ang(nside, index1)
                    radecs = SkyCoord(ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)                  
                    frame = AltAz(obstime=timenow, location=observatory)    
                    altaz = radecs.transform_to(frame)
                    sun_altaz = get_sun(timenow).transform_to(altaz)
                    _sindex = np.logical_and(sun_altaz.alt <= _sun*u.deg,\
                                           altaz.secz <= _airmass)                                       
                    '''
                    if len(_sindex)>0:
                        hp.projplot(theta[_sindex][0],phi[_sindex][0],'.', color = 'grey', coord='C', ms = 5, \
                                    alpha=.2, label='available sky')
                        hp.projplot(theta[_sindex],phi[_sindex],'.', color = 'grey', coord='C', ms = 5, alpha=.05)
                        _title+='%.2f sq deg available'%float(len(_sindex)*_areasingle)
                    else:_title+='0 sq deg available now'
                    '''

                    # make plot
                    plt.title(_title)                              
                    plt.axis('off')
                    plt.legend()
                    plt.savefig(self.fname,dpi=120,format='PNG')
                    plt.close('all')

            except Exception as e:self.response += str(e)
        else:self.response += 'input all parameters needed!'

    def tmap(self):
        self.response = '%s trigger map not found in our DB'%self.argv
        for trigger in self.tlist:
            ck = re.search(trigger,self.command,re.IGNORECASE)
            if ck:
                image_url = 'http://sngroup.oapd.inaf.it/gw/pic/%s/%s_tmap.png'%(ck.group(0),ck.group(0))
                self.response=image_url
                self.attachments = [{"title": "trigger map: %s"%ck.group(0), "image_url": image_url}]

    def visibility(self):

        self.response=''
        if not self.name:
            self.name = 'GWC1'
            self.response += 'no name input, use %s'%self.name

	# current date/time	
	now = datetime.datetime.now()	
	nday=str(now.day)
	if len(nday) == 1:nday="0"+nday	
	nmonth=str(now.month)
	if len(nmonth) == 1:nmonth="0"+nmonth	
        nyear=str(now.year)[2:]

        if self.ra and self.dec:            
            coordlist = self.name+'+'+self.radeg+'+'+self.decdeg
            # retrieving the visibility plot
            image_url = "http://catserver.ing.iac.es/staralt/index.php?action=showImage&form%5Bmode%5D=1&form%5Bday%5D="+nday+"&form%5Bmonth%5D="+nmonth+"&form%5Byear%5D="+str(now.year)+"&form%5Bobservatory%5D=Cerro+Paranal+Observatory+(Chile)&form%5Bcoordlist%5D="+coordlist+"&form%5Bcoordfile%5D=&form%5Bparamdist%5D=2&form%5Bformat%5D=gif&submit=+Retrieve+"
            self.attachments = [{"title": "visibility of %s, %s"%(self.radeg,self.decdeg), "image_url": image_url}]        
        else:self.response += 'Ra and Dec are needed!'

    def extinction(self):

        self.response=''               
        if self.ra and self.dec:            
            ra1 = self.rahms.split(':')
            dec1 = self.decdms.split(':')
            coordlist='lon='+ra1[0]+'%3A'+ra1[1]+'%3A'+ra1[2]+\
                '&lat='+dec1[0]+'%3A'+dec1[1]+'%3A'+dec1[2]
            link="https://ned.ipac.caltech.edu/cgi-bin/calc?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2000.0&"+coordlist+"&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0"
            try:
		tables = pd.read_html(link,header=0)		
		filters=['Landolt','SDSS','UKIRT']
		nt=88		
		nn=0		
		for i in range(nt):
                    if tables[1].iloc[:,0][i] in filters:                       
                        if tables[1].iloc[:,0][i] == 'Landolt' and tables[1].iloc[:,1][i] == 'V':
                            av=tables[1].iloc[:,3][i]
                        nn += 1		                			
		self.response += 'Av = %.2f'%np.round(av,2)
            except Exception as e:
                response = str(e)
        else:self.response += 'Ra and Dec are needed!'

    def asteroid(self):
	
        self.response=''       
        if self.ra and self.dec and self.time:
            if not self.radius:
                self.radius = 10
                self.response += 'no radius found, use %.2f arcsec'%self.radius
            if not self.mag:
                self.mag = 0
                self.response += 'no mag cut found, search all > %.2f'%self.mag

            try:
                # handle time	
                nyear=self.time.split('-')[0]
                nmonth=self.time.split('-')[1]
                nday=self.time.split('-')[2]

                nnday=nday.split('.')
                nhour=float('0.'+nnday[1])*24
                nnhour=str(nhour).split('.')
                nmin=float('0.'+nnhour[1])*60
                nnmin=str(nmin).split('.')
                nsec=float('0.'+nnmin[1])*60

                # calculate julian date	
                t = Time(datetime.datetime(int(nyear), int(nmonth), int(nnday[0]), \
                                           int(nnhour[0]), int(nnmin[0]), int(nsec)))
                tjd=t.jd
	
                # check page	
                pagea = requests.get('https://cgi.minorplanetcenter.net/cgi-bin/checkmp.cgi')
                if pagea.status_code is not 200:
                    self.response+='error: Page not available.'	
                else:                   
                    ra1 = self.rahms.replace(':','+')
                    dec1 = self.decdms.replace(':','+')
                    link="https://minorplanetcenter.net/cgi-bin/mpcheck.cgi?year="+nyear+\
                        "&month="+nmonth+"&day="+nday+"&ra="+ra1+"&decl="+dec1+"&which=pos&TextArea=&radius="+\
                        str(self.radius)+"&limit="+str(self.mag)+"&oc=500&sort=d&mot=h&tmot=s&pdes=u&needed=f&ps=n&type=p"                    
                    self.response += link
                    self.response += '\n'

                    htm=urllib.urlopen(link).read().decode('utf-8')				
                    n1=htm.find("<pre>")
                    n2=htm.find("</pre>")		
                    if n1 == -1:
                        self.response+='Minor planets found: 0\n'
                        self.response+='No known minor planets, brighter than V = '
                        self.response+=str(self.mag)
                        self.response+=', were found in the '
                        self.response+=str(self.radius)
                        self.response+='-arcminute region around the source on '
                        self.response+='%s jd=%s'%(str(self.time), str(tjd))
                    else:                       
                        self.response+='Minor planets found!'
                        self.response+=htm[n1+5:n2]
            except Exception as e:
                self.response+=str(e)
        else:self.response += 'Ra and Dec and Time are needed!'

    def DSS(self):

        self.response=''
	size=5 # [arcmin] the size of the retrieved and saved image
	size1=2 # [arcmin] the size of the image used for the FC
	pix=1.008 # approx pixel scale
        if not self.name:
            self.name = 'GWC1'            

        if self.ra and self.dec:                  
            raS = self.rahms.replace(':','%3A')
            decS = self.decdms.replace(':','%3A')
            link='http://archive.eso.org/dss/dss/image?ra='+raS+'&dec='+decS+\
                '&equinox=J2000&name=&x='+str(size)+'&y='+str(size)+\
                '&Sky-Survey=DSS2-red&mime-type=download-fits&statsmode=WEBFORM'            
            try:
                outf='/tmp/tmpdss.fits'
                self.fname='/tmp/tmpdss.png'
                self.response += 'DSS cutout at %s,%s'%(self.radeg,self.decdeg)
                for dd in [outf,self.fname]:
                    if os.path.exists(dd):os.remove(dd)

                wget.download(link,out=outf)
		fh=fits.open(outf)
		fim = fh[0].data
		fhe = fh[0].header

		# cut image and apply scale
		x1=int(30*(size - size1))
		x2=int(30*(size + size1))
		y1=int(30*(size - size1))
		y2=int(30*(size + size1))
	
		fim=fim[y1:y2,x1:x2]		
		fim[np.isnan(fim)] = 0.0
		transform = AsinhStretch() + PercentileInterval(99.7)
		bfim = transform(fim)
		
		with warnings.catch_warnings():	#because there are deprecated keywords in the header, no need to write it out
                    warnings.simplefilter("ignore")
                    wcs = WCS(fhe)
               
		# produce and save the FC		
		fig=plt.figure(2,figsize=(5,5))
		fig1=fig.add_subplot(111,aspect='equal')		
		plt.imshow(bfim,cmap='gray_r',origin='lower')
		c = Circle((size1*30, size1*30), pix*2, edgecolor='k', lw=2,facecolor='none')
		fig1.add_patch(c)
		c = Circle((size1*30, size1*30), pix*2, edgecolor='yellow', facecolor='none')
		fig1.add_patch(c)                
		txta=fig.text(0.18,0.8,self.name,fontsize=23,color='yellow')
		txta.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
		txtb=fig.text(0.75,0.8,'DSS',fontsize=23,color='yellow')
		txtb.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])	
		fig1.add_patch(FancyArrowPatch((size1*60-15/pix-10,20),(size1*60-10,20),arrowstyle='-',color='k',linewidth=3.5))
		fig1.add_patch(FancyArrowPatch((size1*60-15/pix-10,20),(size1*60-10,20),arrowstyle='-',color='yellow',linewidth=2.0))

		txtc=fig.text(0.75,0.16,'15"',fontsize=20,color='yellow')
		txtc.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])                
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		plt.gca().yaxis.set_major_locator(plt.NullLocator())	
		plt.savefig(self.fname,dpi=120,format='PNG')
		fig.clear()

            except Exception as e:
		self.response+=str(e)
        else:self.response += 'Ra and Dec are needed!'

    def Panstarrs(self):

        self.response=''
        size=1200 #==5 arcmin: the size [px] of the retrieved and saved image
        size1=240 #=1 arcmin: the size [px] of the central part of the image for the FC
        fil='r' #filters	
        pix=0.25
        suc=0          
        if not self.name:
            self.name = 'GWC1'         
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')           
        if self.ra and self.dec:                    
            try:
                fitsurl = geturl(self.radeg, self.decdeg, size=size, filters=fil, format="fits")
                if not fitsurl:	#no image returned, not in the field
                    self.response+='Image not available'
                else:
                    outf='/tmp/tmpps.fits'                    
                    if os.path.exists(outf):os.remove(outf)

                    fh = fits.open(fitsurl[0])
                    fim = fh[0].data
                    fhe = fh[0].header
					
                    # save fits image			
                    with warnings.catch_warnings(): #to prevent the warnings in case the script is run several times and the files are rewritten
                        fits.writeto(outf, fim, fhe,overwrite=True)
						
                    # make the cut and apply scale			
                    x1=int(size/2. - size1/2.)
                    x2=int(size/2. + size1/2.)
                    y1=int(size/2. - size1/2.)
                    y2=int(size/2. + size1/2.)
			
                    fim=fim[y1:y2,x1:x2]

                    fim[np.isnan(fim)] = 0.0
                    transform = AsinhStretch() + PercentileInterval(99.9)
                    bfim = transform(fim)

                    try:	                                        
                        self.fname='/tmp/tmpps.png'
                        self.response += 'Panstarrs cutout at %s,%s'%(self.radeg,self.decdeg)                   
                        if os.path.exists(self.fname):os.remove(self.fname)
			
                        # produce and save the FC				
                        fig=plt.figure(1, figsize=(12,6))
                        fig1=fig.add_subplot(121,aspect='equal')
				                        
                        plt.imshow(bfim,cmap='gray_r',origin='lower')
                        c = Circle((size1/2., size1/2.), 1./pix, edgecolor='k', lw=2,facecolor='none')
                        fig1.add_patch(c)
                        c = Circle((size1/2., size1/2.), 1./pix, edgecolor='yellow', facecolor='none')
                        fig1.add_patch(c)
                        txta=fig.text(0.14,0.8,self.name,fontsize=23,color='yellow')
                        txta.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
                        txtb=fig.text(0.32,0.8,'Pan-STARRS',fontsize=23,color='yellow')
                        txtb.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
				
                        fig1.add_patch(FancyArrowPatch((size1-10/0.25-10,30),(size1-10,30),arrowstyle='-',color='k',linewidth=3.5))
                        fig1.add_patch(FancyArrowPatch((size1-10/0.25-10,30),(size1-10,30),arrowstyle='-',color='yellow',linewidth=2.0))

                        txtc=fig.text(0.44,0.16,'10"',fontsize=20,color='yellow')
                        txtc.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])
				
                        plt.gca().xaxis.set_major_locator(plt.NullLocator())
                        plt.gca().yaxis.set_major_locator(plt.NullLocator())                        

                        # get colour image
                        url = geturl(self.radeg,self.decdeg,size=size1,filters='grz',\
                                     output_size=None,format='png',color=True)
                        r = requests.get(url)

                        im = Image.open(BytesIO(r.content))                            
                        figA1=fig.add_subplot(122,aspect='equal')
				
                        plt.imshow(im)
                        c = Circle((size1/2., size1/2.), 1./pix, lw=2, edgecolor='k', facecolor='none')
                        figA1.add_patch(c)
                        c = Circle((size1/2., size1/2.), 1./pix, edgecolor='yellow', facecolor='none')
                        figA1.add_patch(c)                      

                        plt.gca().xaxis.set_major_locator(plt.NullLocator())
                        plt.gca().yaxis.set_major_locator(plt.NullLocator())
                        plt.subplots_adjust(wspace=0.1)
                        
                        plt.savefig(self.fname,dpi=120,format='PNG')
                        fig.clear()
                    except Exception as e:
                        self.response+=str(e)
            except Exception as e1:
                self.response+=str(e1)
        else:self.response += 'Ra and Dec are needed!'

    def TNS(self):
        self.response=''
        if self.ra and self.dec and self.radius:              
            datestart='2001-01-01' # start date for searching in TNS server
            self.response += 'search from %s till now\n'%datestart
            now=datetime.datetime.now()
            now=now.strftime("%Y-%m-%d")
            ra1=self.rahms.replace(':','%3A')
            dec1=self.decdms.replace(':','%3A')
	
            pagea = requests.get('https://wis-tns.weizmann.ac.il/search')
            if pagea.status_code is not 200:
		self.response+='error: Page not available.'		
            
            else:
		link="https://wis-tns.weizmann.ac.il/search?&name=&name_like=0&isTNS_AT=all&public=all&unclassified_at=0&classified_sne=0&ra="+ra1+"&decl="+dec1+ "&radius="+str(self.radius)+"&coords_unit=arcsec&groupid%5B%5D=null&classifier_groupid%5B%5D=null&type%5B%5D=null&date_start%5Bdate%5D="+datestart+"&date_end%5Bdate%5D="+now+"&discovery_mag_min=&discovery_mag_max=&internal_name=&redshift_min=&redshift_max=&spectra_count=&discoverer=&classifier=&discovery_instrument%5B%5D=&classification_instrument%5B%5D=&hostname=&associated_groups%5B%5D=null&ext_catid=&num_page=50&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0"	

		#a silly piece of code: do this because the table is not easily parsable :(
		fp = urllib.urlopen(link)
		mybytes = fp.read()               
		mystr = mybytes.decode("utf8")
		fp.close()
		pos=mystr.find('class=\"count rsults\"')               
		if pos>0:			
                    pos=int(mystr.find('out of <em class="placeholder">') + len('out of <em class="placeholder">'))                    
                    mystr1=mystr[pos:]
                    pos=mystr1.find('<')
                    numobj=int(mystr1[:pos])

                    self.response+='Found %s sources in TNS at %s,%s with radius %.2f arcsec: \n'%\
                        (str(numobj),self.radeg,self.decdeg,self.radius)
                    _total = 120
                    self.response += 'TNS'.ljust(_total-6)
                    self.response += 'seperation(arcsec) \n'
                    #search for each event by class="cell-name"				
                    iii=-1
                    for m in re.finditer('class=\"cell-name\">',mystr):
                        iii += 1
                        if iii==0:	#the first occurrence is skipped
                            continue 
                        mystr1=mystr[m.end()+9:]
                        lin=mystr1[:mystr1.find('\"')]
                        nam=mystr1[mystr1.find('\"')+2:mystr1.find('</a>')]
                        mystr2=mystr1[mystr1.find('class=\"cell-ra\">')+len('class=\"cell-ra\">'):]
                        raS=mystr2[:mystr2.find('</td>')]
                        mystr3=mystr2[mystr2.find('class=\"cell-decl\">')+len('class=\"cell-decl\">'):]
                        decS=mystr3[:mystr3.find('</td>')]                        

                        raSdeg, decSdeg, raShms, decSdms = readradec(raS,decS)                          
                        c1=SkyCoord(self.radeg,self.decdeg,unit='deg')
                        c2=SkyCoord(raSdeg,decSdeg,unit='deg')
                        sep=c1.separation(c2)

                        _tns = '%s %s %s %s'%(nam,raS,decS,lin)
                        self.response += _tns.ljust(_total-len(_tns))
                        self.response += '%.2f arcsec \n'%(sep.arcsec)		
                else:		
                    self.response+=('Found 0 sources in TNS at %s,%s with radius %.2f arcsec'%\
                               (self.radeg,self.decdeg,self.radius))
        else:self.response += 'Ra and Dec and radius are needed!'

    def Skymappercat(self):  
        self.response='tbd'

    def Skymapperimage(self):  
        self.response='tbd'
      
    def simbad(self):
        self.response=''        
        if self.ra and self.dec and self.radius:                   
            pagea = requests.get('http://simbad.u-strasbg.fr/simbad/')
            if pagea.status_code is not 200:
		self.response = 'error: Page not available.'	
            else:                
		with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = Simbad.query_region(SkyCoord(float(self.radeg), float(self.decdeg), \
                            unit=(u.deg, u.deg), frame='icrs'),radius=Angle(float(self.radius),"arcsec"))                                                
                try:
                    tblA=result['MAIN_ID']                   
                    n=len(tblA)                    
                    self.response += 'Number of sources: %i (with radius=%i arcsec) \n'%\
                                     (n,float(self.radius))                    
                    raS=result['RA']
                    decS=result['DEC']                        
                    _total = 100
                    self.response += 'Simbad'.ljust(_total-6)
                    self.response += 'seperation(arcsec) \n'
                    for i in range(n):
                        _radeg, _decdeg, _rahms, _decdms = readradec(raS[i],decS[i])
                        c1=SkyCoord(self.radeg,self.decdeg,unit='deg')
                        c2=SkyCoord(_radeg,_decdeg,unit='deg')                        
                        sep=c1.separation(c2)
                        self.response += tblA[i].decode("utf-8").ljust(_total-len(tblA[i].decode("utf-8")))
                        self.response += '%.2f arcsec \n'%sep.arcsec               
                except:
                    self.response += 'Number of sources: 0 with radius=%i \n'%float(self.radius)
        else:self.response += 'Ra and Dec and radius are needed!'

    def GLADE(self):
        _limnum = 50
        self.response=''
        _filter = {}
        _raname,_decname,_magname,_distname,_columns,_cat = 'RAJ2000', 'DEJ2000', 'BMAG', 'Dist',\
                        ['RAJ2000', 'DEJ2000', 'Dist', 'BMAG','PGC', 'GWGC', 'HyperLEDA', '2MASS'],\
                        'VII/281'     
        if self.ral and self.decl and self.magl and self.distl:
            try:
                ramin,ramax=min(eval(self.ral)),max(eval(self.ral))
                _filter[_raname] = '%s..%s'%(str(ramin),str(ramax))
            except:self.response+='input ra format wrong\n'                       
            try:
                decmin,decmax=min(eval(self.decl)),max(eval(self.decl))
                _filter[_decname] = '%s..%s'%(str(decmin),str(decmax))
            except:self.response+='input dec format wrong\n'      
            try:           
                bmagmin,bmagmax=min(eval(self.magl)),max(eval(self.magl))
                _filter[_magname] = '%s..%s'%(str(bmagmin),str(bmagmax))
            except:self.response+='input mag format wrong\n'        
            try:
                distmin,distmax=min(eval(self.distl)),max(eval(self.distl))
                _filter[_distname] = '%s..%s'%(str(distmin),str(distmax))
            except:self.response+='input dist format wrong\n' 
            v = Vizier(columns=_columns,column_filters=_filter)
            v.ROW_LIMIT = -1
            try:
                catalogs = v.get_catalogs(_cat)[0]   
                self.response += "%i galaxies selected \n"%len(catalogs)

                if len(catalogs)<_limnum:       
                    for _c in catalogs.columns:
                        if catalogs[_c].dtype == 'float64':
                            catalogs[_c].format = "%6.3f"
                    self.response+=str(catalogs).replace(' -- ','-'*6+' ').replace(' --- ','-'*6+' ')
                else:self.response+='I will not show them all because they are too many (>%i)'%_limnum
            except:
                self.response+='No galaxies found!'
        else:self.response += 'Ral and Decl and magl and distl are all needed!'

    def ESO(self):
        self.response=''
        if self.ra and self.dec:            
            if not self.radius:
                self.radius = 10
                self.response += 'no radius found, use %.2f arcsec\n'%self.radius
            linkESO='http://archive.eso.org/wdb/wdb/eso/eso_archive_main/query?ra='+\
                self.rahms.replace(':','+')+\
                '&dec='+self.decdms.replace(':','+')+\
                '&deg_or_hour=hours&box=%.2f&max_rows_returned=2000'%self.radius
            if self.instr:linkESO+='&instrument=%s'%self.instr
            self.response+=linkESO        
        else:self.response += 'Ra and Dec are needed!'

    def NED(self):
        self.response=''
        if self.ra and self.dec:
            if not self.radius:
                self.radius = 10
                self.response += 'no radius found, use %.2f arcsec\n'%self.radius           
            linkNED='http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?'+\
                'in_csys=Equatorial&in_equinox=J2000.0&lon='+self.radeg+\
                '&lat='+self.decdeg+'&radius=%.2f'%self.radius+\
                '&search_type=Near+Position+Search'
            self.response+=linkNED
        else:self.response += 'Ra and Dec are all needed!'

    def contact(self):
        print self.args
        self.response = 'Thanks for your contact, I will check soon!'

    def runsex(self,img,_dir=_pstpath):        
        cdef = open(_dir+"/default/default.param").readlines()    
        lparams = []
        for r in cdef:
            if r[0] != '#' and len(r.strip())>0: lparams.append(r.split()[0])
        fwlim = np.array([2.0, 2.75, 3.25, 3.75, 4.5, 6.,8., 10.])
        gconv = ['2.5_5x5', '3.0_7x7','3.5_7x7','4.0_7x7','5.0_9x9',\
                 '7.0_17x17','9.0_17x17']
        if self.fwhm<fwlim[0]:
            usefilter = "N"
            gauss = 'gauss_'+gconv[0]+'.conv'           
        elif self.fwhm>fwlim[-1]:
            usefilter = "N"
            gauss = 'gauss_'+gconv[-1]+'.conv'           
        else:
            i = len(np.compress(fwlim<self.fwhm,fwlim))
            gauss = 'gauss_'+gconv[i-1]+'.conv'
            usefilter = "Y"

        sex_verbose = 'CRITICAL'        
        lconfig = {
            "STARNNW_NAME":_dir+"/default/default.nnw",
            "FILTER_NAME":_dir+'/default/'+gauss,
            "PIXEL_SCALE":self.pixel_scale,
            "PHOT_APERTURES":self.phot_aperture,
            "ANALYSIS_THRESH":self.analysis_thresh*self.threshold,
            "DETECT_MINAREA":self.detect_minarea*self.fwhm,
            "DETECT_THRESH":str(self.threshold),
            "BACK_SIZE":self.back_size,
            "MAG_ZEROPOINT":self.mag_zeropoint,
            "SATUR_LEVEL":self.satur_level,
            "SEEING_FWHM":self.seeing,}
        sew = sewpy.SEW(params=lparams,config=lconfig,loglevel=sex_verbose)
        try:out = sew(img+'.fits')
        except:out = sew(img)
        return out['table']

##################
def escape(word):   
    replace_with = {
        '&gt;'  : '<',
        '&lt;'  : '>',
        '&amp;' : '&',
        '&quot;': '"', 
        '&#39'  : "'"}   
    for oo in replace_with: word = word.replace(oo,replace_with[oo])
    return word

def readradec(ra,dec):

    if ':' in ra:_p = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    elif len(ra.split())>1:_p = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    else:_p = SkyCoord(ra,dec,unit=(u.deg,u.deg))   

    rasec = '{0:0{width}}'.format(float("%.2f"%_p.ra.hms[2]), width=5)   
    rahms = '%.2i:%.2i:%s'%(_p.ra.hms[0],abs(_p.ra.hms[1]),rasec)
    decdms = '%.2i:%.2i:%.2f'%(_p.dec.dms[0],abs(_p.dec.dms[1]),abs(_p.dec.dms[2]))
    radeg = '%.5f'%_p.ra.deg
    decdeg = '%.5f'%_p.dec.deg
    return radeg, decdeg, rahms, decdms

'''
This piece of code is taken from: https://ps1images.stsci.edu/ps1image.html
'''
def getimages(ra,dec,size=1200,filters="grizy"):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table

def geturl(ra, dec, size=1200, output_size=None, filters="grizy", format="jpg", color=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url

def download_dss(ra,dec,width, heaigh,clobber,verbose):            
    _fname = 'tmp.in'           
    _ff = open(_fname,'w')           
    _ff.write('tmp %s %s %s %s'%(ra.replace(':',' '),dec.replace(':',' '),width,heaigh))
    _ff.close()
    pid = subprocess.Popen(shlex.split('dss2 red -i %s'%_fname),\
                           stdout=subprocess.PIPE)
    output,error = pid.communicate()
    nn,_suffix=-99,None
    for nii,ii in enumerate(output.split('\n')):
        if len(ii)==0:continue
        if ii[0] == '-':nn=nii+1                    
        if nn>0 and nii==nn:_suffix = ii.split()[8].lower()
    if not _suffix is None:           
        os.rename('tmp%s.fits'%_suffix,'tmpdss.fits')
        return 'tmpdss.fits'
    else:return None

def zscale(image, nsamples=1000, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5, max_iterations=5):

    # Sample the image
    image = np.asarray(image)
    image = image[np.isfinite(image)]
    stride = int(max(1.0, image.size / nsamples))
    samples = image[::stride][:nsamples]
    samples.sort()
    
    npix = len(samples)
    zmin = samples[0]
    zmax = samples[-1]

    # Fit a line to the sorted array of samples
    minpix = max(min_npixels, int(npix * max_reject))
    x = np.arange(npix)
    ngoodpix = npix
    last_ngoodpix = npix + 1

    # Bad pixels mask used in k-sigma clipping
    badpix = np.zeros(npix, dtype=bool)

    # Kernel used to dilate the bad pixels mask
    ngrow = max(1, int(npix * 0.01))
    kernel = np.ones(ngrow, dtype=bool)

    for niter in range(max_iterations):
        if ngoodpix >= last_ngoodpix or ngoodpix < minpix:
            break

        fit = np.polyfit(x, samples, deg=1, w=(~badpix).astype(int))
        fitted = np.poly1d(fit)(x)
        
        # Subtract fitted line from the data array
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        threshold = krej * flat[~badpix].std()

        # Detect and reject pixels further than k*sigma from the fitted line
        badpix[(flat < - threshold) | (flat > threshold)] = True

        # Convolve with a kernel of length ngrow
        badpix = np.convolve(badpix, kernel, mode='same')

        last_ngoodpix = ngoodpix
        ngoodpix = np.sum(~badpix)

    slope, intercept = fit

    if ngoodpix >= minpix:
        if contrast > 0:
            slope = slope / contrast
        center_pixel = (npix - 1) // 2
        median = np.median(samples)
        zmin = max(zmin, median - (center_pixel - 1) * slope)
        zmax = min(zmax, median + (npix - center_pixel) * slope)
    return zmin,zmax

def col2RA(column,hdr):
    crpix1 = hdr['crpix1']
    crval1 = hdr['crval1']
    cd1 = hdr['cd1_1']
    crpix2 = hdr['crpix2']
    crval2 = hdr['crval2']
    cd2 = hdr['cd2_2']
    return crval1 + (column-crpix1)*cd1

def row2dec(row,hdr):
    crpix1 = hdr['crpix1']
    crval1 = hdr['crval1']
    cd1 = hdr['cd1_1']
    crpix2 = hdr['crpix2']
    crval2 = hdr['crval2']
    cd2 = hdr['cd2_2']
    return crval2 + (row-crpix2)*cd2

def parse_slack_output(slack_rtm_output):
    """
        The Slack Real Time Messaging API is an events firehose.
        this parsing function returns None unless a message is
        directed at the Bot, based on its ID.
    """
    output_list = slack_rtm_output
    command,channel,user = None, None,None   
    if output_list and len(output_list) > 0:                
        for output in output_list:                                      
            if output and 'text' in output and AT_BOT in output['text']:                                
                command,channel,user = output['text'].split(AT_BOT)[1].strip(), \
                                       output['channel'],output['user']
    return command,channel,user

if __name__ == "__main__":

    READ_WEBSOCKET_DELAY = 0 # 1 second delay between reading from firehose
    if slack_client.rtm_connect():
        print "Bot connected and running!"      
        while True:
            try: command, channel, user = parse_slack_output(slack_client.rtm_read())
            except Exception as e:
                slack_client.api_call("chat.postMessage", channel=channel,
                                      text=str(e), as_user=True)         
            if command and channel:   

                try:print command,':\tuser:',user,'/',_usrl[user],'\tchannel:',channel,'/',channel_list[channel]
                except:print command,':\tuser:',user,'\tchannel:',channel            

#                if not channel in channel_list.keys():
#                    slack_client.api_call("chat.postMessage", channel=channel,
#                                          text=msg1, as_user=True)  

#                elif not _usrl[user] in _susers:
#                    slack_client.api_call("chat.postMessage", channel=channel,
#                                          text=msg2, as_user=True)

#                else:                    
                if True:
                    # main code
                    command = escape(command)                   
                    try: 
                        logging.info('%s %s : %s'%(_usrl[user],\
                                               channel_list[channel],command))                    
                    except:
                        logging.info('%s %s : %s'%(user,channel,command))

                    try:
                        handle_command(command, channel)
                    except Exception as e:
                        slack_client.api_call("chat.postMessage", channel=channel,
                                              text=str(e), as_user=True)                                           
            time.sleep(READ_WEBSOCKET_DELAY)
    else:
        print("Connection failed. Invalid Slack token or bot ID?")
