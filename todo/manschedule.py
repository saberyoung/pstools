"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input
import pst,os,sys
import healpy as hp
import numpy as np
from pst.default import *

def man_search(_fits, _tel, _verbose):  

    # read telescopes
    _config = pst.load_config()   
    if not _config: sys.exit()
    arglist,optlist = _config['general']['params'],\
                      _config['telescope']['params']
    _dir =  arglist['data']['dir']
    if _tel is None:
        _tel = arglist['react']['telescope']

    _info = '>>> telescope used: %s \n'%_tel
    if _verbose: print(_info)      

    # read params from configure file
    _paramslist = {}
    _paramslist['tmp'] = {}
    for tel0 in _tel.split(','):
        _config = pst.load_config(tel0)
        _paramslist['arg'] = _config['general']['params']
        _paramslist[tel0] = _config['telescope']['params']
        _info = ' - Read params for telescope:%s'%tel0
        if _verbose: print(_info)        

    # define email, slack, phone, ...
    _paramslist['tmp']['files'] = [_fits]
    _paramslist['arg']['email']['emailcontent']='offline alert \n'
    _paramslist['arg']['phone']['phonecontent']='offline alert: '
    _paramslist['arg']['slack']['slackcontent']='*offline* alert:\n'
    
    # decide prioritization method
    _trigger,_gal,_dist = eval(_paramslist['arg']['priorization']['trigger']),\
                          int(_paramslist['arg']['priorization']['galaxy']),\
                          eval(_paramslist['arg']['priorization']['dist'])

    # if activate trigger, however no fits given
    if _trigger and _fits is None: 
        _fits = input('!!! in par file, you specify using trigger, '+\
                      'however, not given exactly.\n\t'+\
                      'input map name:')
    # if exists
    if _fits is not None and not os.path.exists(_fits): 
        sys.exit('### Error: %s not exists...'%_fits)

    if _trigger and _fits is not None:   # check format and get fits map if given
        _paramslist['tmp']['tmap'], _paramslist['tmp']['distmu'], \
            _paramslist['tmp']['distsigma'], _paramslist['tmp']['distnorm'], \
            _paramslist['tmp']['header'], _paramslist['tmp']['voevent'] = \
                pst.get_hp_map(_fits,_paramslist['arg']['show']['verbose'],\
                               coord,nside,_dir)
    else:
        _paramslist['tmp']['tmap'], _paramslist['tmp']['distmu'], \
            _paramslist['tmp']['distsigma'], _paramslist['tmp']['distnorm'], \
            _paramslist['tmp']['header'], _paramslist['tmp']['voevent'] = \
            np.zeros(12*nside**2),None,None,None,None,None

    # store above infos
    _info = ''
    if _trigger:  _info += 'trigger '
    if _gal == 1: _info += '& galaxy lums'
    if _gal == 2: _info += '& galaxy counts'
    if _dist: 
        if _trigger and _gal>0:
            _info += '& galaxy dist'
        else:
            _paramslist['arg']['priorization']['dist']='False'
    if _verbose: print('\npriorization:',_info,'\n')

    # main program
    pst.main(_paramslist)
