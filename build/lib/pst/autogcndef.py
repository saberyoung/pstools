"""############################################################################ 
2019/1/30 Start
Extract map and infos from xml, and call main()
""" ############################################################################
from __future__ import print_function
from builtins import input
import os,sys,glob,gcn,gcn.handlers,gcn.notice_types,scipy.stats,voeventparse,astropy.time
import healpy as hp
import numpy as np
import math as mt
import pst
from pst.default import *

_log = False # record log file

#################################################
""" 
Function to call every time a GCN is received.
alert types are defined by pygcn, options below:
    GRB_COORDS=1,
    TEST_COORDS=2,
    IM_ALIVE=3,
    KILL_SOCKET=4,
    MAXBC=11,
    BRAD_COORDS=21,
    GRB_FINAL=22,
    HUNTS_SRC=24,
    ALEXIS_SRC=25,
    XTE_PCA_ALERT=26,
    XTE_PCA_SRC=27,
    XTE_ASM_ALERT=28,
    XTE_ASM_SRC=29,
    COMPTEL_SRC=30,
    IPN_RAW=31,
    IPN_SEG=32,
    SAX_WFC_ALERT=33,
    SAX_WFC_SRC=34,
    SAX_NFI_ALERT=35,
    SAX_NFI_SRC=36,
    XTE_ASM_TRANS=37,
    spare38=38,
    IPN_POS=39,
    HETE_ALERT_SRC=40,
    HETE_UPDATE_SRC=41,
    HETE_FINAL_SRC=42,
    HETE_GNDANA_SRC=43,
    HETE_TEST=44,
    GRB_CNTRPART=45,
    SWIFT_TOO_FOM=46,
    SWIFT_TOO_SC_SLEW=47,
    DOW_TOD=48,
    spare50=50,
    INTEGRAL_POINTDIR=51,
    INTEGRAL_SPIACS=52,
    INTEGRAL_WAKEUP=53,
    INTEGRAL_REFINED=54,
    INTEGRAL_OFFLINE=55,
    INTEGRAL_WEAK=56,
    AAVSO=57,
    MILAGRO_POS=58,
    KONUS_LC=59,
    SWIFT_BAT_GRB_ALERT=60,
    SWIFT_BAT_GRB_POS_ACK=61,
    SWIFT_BAT_GRB_POS_NACK=62,
    SWIFT_BAT_GRB_LC=63,
    SWIFT_BAT_SCALEDMAP=64,
    SWIFT_FOM_OBS=65,
    SWIFT_SC_SLEW=66,
    SWIFT_XRT_POSITION=67,
    SWIFT_XRT_SPECTRUM=68,
    SWIFT_XRT_IMAGE=69,
    SWIFT_XRT_LC=70,
    SWIFT_XRT_CENTROID=71,
    SWIFT_UVOT_DBURST=72,
    SWIFT_UVOT_FCHART=73,
    SWIFT_BAT_GRB_LC_PROC=76,
    SWIFT_XRT_SPECTRUM_PROC=77,
    SWIFT_XRT_IMAGE_PROC=78,
    SWIFT_UVOT_DBURST_PROC=79,
    SWIFT_UVOT_FCHART_PROC=80,
    SWIFT_UVOT_POS=81,
    SWIFT_BAT_GRB_POS_TEST=82,
    SWIFT_POINTDIR=83,
    SWIFT_BAT_TRANS=84,
    SWIFT_XRT_THRESHPIX=85,
    SWIFT_XRT_THRESHPIX_PROC=86,
    SWIFT_XRT_SPER=87,
    SWIFT_XRT_SPER_PROC=88,
    SWIFT_UVOT_POS_NACK=89,
    SWIFT_BAT_ALARM_SHORT=90,
    SWIFT_BAT_ALARM_LONG=91,
    SWIFT_UVOT_EMERGENCY=92,
    SWIFT_XRT_EMERGENCY=93,
    SWIFT_FOM_PPT_ARG_ERR=94,
    SWIFT_FOM_SAFE_POINT=95,
    SWIFT_FOM_SLEW_ABORT=96,
    SWIFT_BAT_QL_POS=97,
    SWIFT_BAT_SUB_THRESHOLD=98,
    SWIFT_BAT_SLEW_POS=99,
    AGILE_GRB_WAKEUP=100,
    AGILE_GRB_GROUND=101,
    AGILE_GRB_REFINED=102,
    SWIFT_ACTUAL_POINTDIR=103,
    AGILE_POINTDIR=107,
    AGILE_TRANS=108,
    AGILE_GRB_POS_TEST=109,
    FERMI_GBM_ALERT=110,
    FERMI_GBM_FLT_POS=111,
    FERMI_GBM_GND_POS=112,
    FERMI_GBM_LC=113,
    FERMI_GBM_GND_INTERNAL=114,
    FERMI_GBM_FIN_POS=115,
    FERMI_GBM_ALERT_INTERNAL=116,
    FERMI_GBM_FLT_INTERNAL=117,
    FERMI_GBM_TRANS=118,
    FERMI_GBM_POS_TEST=119,
    FERMI_LAT_POS_INI=120,
    FERMI_LAT_POS_UPD=121,
    FERMI_LAT_POS_DIAG=122,
    FERMI_LAT_TRANS=123,
    FERMI_LAT_POS_TEST=124,
    FERMI_LAT_MONITOR=125,
    FERMI_SC_SLEW=126,
    FERMI_LAT_GND=127,
    FERMI_LAT_OFFLINE=128,
    FERMI_POINTDIR=129,
    SIMBADNED=130,
    FERMI_GBM_SUBTHRESH=131,
    SWIFT_BAT_MONITOR=133,
    MAXI_UNKNOWN=134,
    MAXI_KNOWN=135,
    MAXI_TEST=136,
    OGLE=137,
    CBAT=138,
    MOA=139,
    SWIFT_BAT_SUBSUB=140,
    SWIFT_BAT_KNOWN_SRC=141,
    VOE_11_IM_ALIVE=142,
    VOE_20_IM_ALIVE=143,
    FERMI_SC_SLEW_INTERNAL=144,
    COINCIDENCE=145,
    FERMI_GBM_FIN_INTERNAL=146,
    SUZAKU_LC=148,
    SNEWS=149,
    LVC_PRELIMINARY=150,
    LVC_INITIAL=151,
    LVC_UPDATE=152,
    LVC_TEST=153,
    LVC_COUNTERPART=154,
    AMON_ICECUBE_COINC=157,
    AMON_ICECUBE_HESE=158,
    CALET_GBM_FLT_LC=160,
    CALET_GBM_GND_LC=161,
    AMON_ICECUBE_EHE=169
"""

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,                
    gcn.notice_types.LVC_UPDATE)   

def process_gcn(payload, root):  
   
    # read arglist
    # decide which telescopes to be activated
    arglist,optlist = pst.load_config()
    _tell = arglist['react']['telescope']
    _dir =  arglist['data']['dir']
    if _tell is None: return('!!! Error: option tel wrong ...')    

    # judge if it's new voevet
    # if not return
    _voname = '%s/%s.xml'%(_dir,os.path.basename(root.attrib['ivorn']))
    if os.path.exists(_voname):
        if eval(arglist['show']['verbose']): print ('%s exists'%_voname)
        return

    if _log: # if record
        import logging

    # for python3, transform from bytes to string
    if sys.version_info>(3,0,0):
        payload = str(payload, encoding = "utf-8")

    # if new voevent
    # red alert, store voevet
    print(' - [Alert] -: %s'%_voname)
    with open(_voname,'w') as _vo:
        if eval(arglist['show']['verbose']): 
            print ('create %s'%_voname)
        _vo.write(payload)

    # read params from configure file
    _paramslist = {}
    _paramslist['tmp'] = {}
    for tel0 in _tell.split(','):        
        arglist,optlist = pst.load_config(tel0)           
        _paramslist['arg'] = arglist
        _paramslist[tel0] = optlist
        _info = '>>> Read params for telescope:%s'%tel0
        if _log:logging.info(_info)
        if _paramslist['arg']['show']['verbose']: print(_info)        

    # record role in email list   
    if root.attrib['role'] == 'test':
        if not eval(_paramslist['arg']['react']['test']):return

    # define email, slack, phone, ...
    _paramslist['tmp']['files'] = [_voname]
    _paramslist['arg']['email']['emailcontent']='online %s alert \n'%root.attrib['role']
    _paramslist['arg']['phone']['phonecontent']='online %s alert: '%root.attrib['role']
    _paramslist['arg']['slack']['slackcontent']='online %s alert: '%root.attrib['role']

    # for auto search, trigger is forced to be activated
    _paramslist['arg']['priorization']['trigger'] = 'True'
                        
    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    _paramslist['tmp']['voevent'] = params

    # store infos
    try: mapurl = params['skymap_fits']
    except: mapurl = False
    if mapurl:   # download map
        try: _fits = pst.get_skymap(mapurl,os.path.basename(root.attrib['ivorn']),_dir)
        except: return('### Error: no ivorn in %s, TB checked'%_fits)
        try:
            (_paramslist['tmp']['tmap'], _paramslist['tmp']['distmu'], \
             _paramslist['tmp']['distsigma'], _paramslist['tmp']['distnorm']), \
                _paramslist['tmp']['header'] = hp.read_map(_fits, \
                    field=[0, 1, 2, 3],h=True, \
                    verbose=_paramslist['arg']['show']['verbose'])
        except:
            _paramslist['tmp']['tmap'], _paramslist['tmp']['header'] = \
                                hp.read_map(_fits, h=True, \
                                verbose=_paramslist['arg']['show']['verbose'])
            _paramslist['tmp']['distmu'], _paramslist['tmp']['distsigma'], \
                _paramslist['tmp']['distnorm'] = None, None, None
    else:    # generate map
        print ('### Warning: no skymap_fits found, try bulding...')
        _paramslist['tmp']['tmap'], _paramslist['tmp']['header'] = \
                pst.build_hp_map(root,_dir+'/'+\
                os.path.basename(root.attrib['ivorn'])+'.fits',\
                nside,_coord=_paramslist['arg']['show']["coord"])
        if len(_paramslist['tmp']['tmap']) > 0:
            print ('### Warning: genearted a fits,'+\
                   ' %s'%(root.attrib['ivorn']+'.fits'))
            distmu, distsigma, distnorm = None, None, None
        else:return('### Error: failed to build fits') 

    # main process    
    pst.main(_paramslist)
