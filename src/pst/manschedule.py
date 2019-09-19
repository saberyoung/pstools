"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input
import pst,os,sys

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
    for tel0 in _tel.split(','):        
        arglist,optlist = pst.load_config(tel0)           
        _paramslist['arg'] = arglist
        _paramslist[tel0] = optlist
        _info = '>>> Read params for telescope:%s'%tel0
        if _log:logging.info(_info)
        if _verbose: print(_info)        

    # define email content        
    _paramslist['arg']['email']['files'] = []
    _paramslist['arg']['email']['images'] = []
    _paramslist['arg']['email']['emailsub']+='[offline]'
    _paramslist['arg']['email']['emailcontent']='offline GW search\n'
    _paramslist['arg']['phone']['phonecontent']='offline GW search\n'    

    # decide prioritization method        
    _trigger,_mass,_dist,_ngal = eval(_paramslist['arg']['priorization']['trigger']),\
                                 eval(_paramslist['arg']['priorization']['mass']),\
                                 eval(_paramslist['arg']['priorization']['dist']),\
                                 eval(_paramslist['arg']['priorization']['number'])

    if _trigger and _fits is None: _fits = input('specify an input map name:')
    if not os.path.exists(_fits): sys.exit('### Error: %s not exists...'%_fits)
    if _trigger or _mass or _dist or _ngal:  _paramslist['arg']['search'] = 'man trigger'             
    else: _paramslist['arg']['search'] = 'man normal'

    if _verbose: print('>'*5,_paramslist['arg']['search'])  
    if _log: logging.info(_paramslist['arg']['search'])
    pst.main(_fits, _paramslist)
