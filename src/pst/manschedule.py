"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input
import pst,os,sys
import healpy as hp
import numpy as np

def man_search(_fits, _tel, _verbose, _log):  

    # read telescopes
    if _tel is None:
        arglist,optlist = pst.load_config()
        _tel = arglist['observe']['telescope']
        if _tel is None: sys.exit('Error: option tel wrong ...')    

    _info = '>>> telescope used: %s'%_tel
    if _log: 
        import logging
        logging.info(_info)
    if _verbose: print(_info)      

    # read params from configure file
    _paramslist = {}
    _paramslist['tmp'] = {}
    for tel0 in _tel.split(','):        
        arglist,optlist = pst.load_config(tel0)           
        _paramslist['arg'] = arglist
        _paramslist[tel0] = optlist
        _info = '>>> Read params for telescope:%s'%tel0
        if _log:logging.info(_info)
        if _verbose: print(_info)        

    # define email, slack, phone, ...
    _paramslist['tmp']['files'] = []
    _paramslist['tmp']['images'] = []    
    _paramslist['arg']['email']['emailcontent']='offline alert \n'
    _paramslist['arg']['phone']['phonecontent']='offline alert: '  
    _paramslist['arg']['slack']['slackcontent']='offline alert: '

    # decide prioritization method        
    _trigger,_mass,_dist,_ngal = eval(_paramslist['arg']['priorization']['trigger']),\
                                 eval(_paramslist['arg']['priorization']['mass']),\
                                 eval(_paramslist['arg']['priorization']['dist']),\
                                 eval(_paramslist['arg']['priorization']['number'])

    # if activate trigger, however no fits given
    if _trigger and _fits is None: 
        _fits = input('!!! in par file, you specify need trigger, however, not given exactly.\n\t'+\
                      'input map name:')
    # if exists
    if _fits is not None and not os.path.exists(_fits): sys.exit('### Error: %s not exists...'%_fits)

    if _trigger and _fits is not None:   # check format and get fits map if given
        _paramslist['tmp']['tmap'], _paramslist['tmp']['distmu'], \
            _paramslist['tmp']['distsigma'], _paramslist['tmp']['distnorm'], \
            _paramslist['tmp']['header'], _paramslist['tmp']['voevent'] = \
                pst.get_hp_map(_fits,_paramslist['arg']['show']['verbose'],\
                _paramslist['arg']['priorization']['nside'])
    else:
        _paramslist['tmp']['tmap'], _paramslist['tmp']['distmu'], \
            _paramslist['tmp']['distsigma'], _paramslist['tmp']['distnorm'], \
            _paramslist['tmp']['header'], _paramslist['tmp']['voevent'] = \
            np.zeros(12*(int(_paramslist['arg']['priorization']['nside']))**2),\
            None,None,None,None,None

    # store above infos
    if _trigger:  _paramslist['tmp']['search'] = 'trigger'
    elif _mass or _dist or _ngal: _paramslist['tmp']['search'] = 'galaxy'
    else: _paramslist['tmp']['search'] = 'normal'

    if _verbose: print('>'*5,_paramslist['tmp']['search'])  
    if _log: logging.info(_paramslist['tmp']['search'])

    # main program
    pst.main(_paramslist)
