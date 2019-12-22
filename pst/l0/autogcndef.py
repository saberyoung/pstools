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
alert types are defined by pygcn, check this webpage: 
https://github.com/lpsinger/pygcn/blob/master/gcn/notice_types.py
"""

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,                
    gcn.notice_types.LVC_UPDATE,
    gcn.notice_types.LVC_RETRACTION)   

def process_gcn(payload, root):  

    # read arglist
    # decide which telescopes to be activated
    _config = pst.load_config()
    if not _config: return
    arglist,optlist = _config['general']['params'],\
                      _config['telescope']['params']
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
        _config = pst.load_config(tel0)    
        if not _config: return
        arglist,optlist = _config['general']['params'],\
                          _config['telescope']['params']
        _paramslist['arg'] = arglist
        _paramslist[tel0] = optlist
        _info = '>>> Read params for telescope:%s'%tel0
        if _log:logging.info(_info)
        if _paramslist['arg']['show']['verbose']: print(_info)        

    # record role in email list   
    if root.attrib['role'] == 'test':
        if not test:return

    # define email, slack, phone, ...
    _paramslist['tmp']['files'] = [_voname]
    _paramslist['arg']['email']['emailcontent']='online %s alert \n'%root.attrib['role']
    _paramslist['arg']['phone']['phonecontent']='online %s alert: '%root.attrib['role']
    _paramslist['arg']['slack']['slackcontent']='online %s alert: \n'%root.attrib['role']

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
                    nside,_coord=coord)
        if len(_paramslist['tmp']['tmap']) > 0:
            print ('### Warning: genearted a fits,'+\
                   ' %s'%(root.attrib['ivorn']+'.fits'))
            distmu, distsigma, distnorm = None, None, None
        else:return('### Error: failed to build fits') 

    # main process    
    pst.main(_paramslist)
