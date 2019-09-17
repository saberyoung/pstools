"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input

def man_search(hpxfits, tel):  

    # read params
    _opts_list = {}
    for tel0 in _tel.split(','):
        _info = 'man mode with telescope:%s'%tel0
        print('## %s'%_info)
        logging.info(_info)
        arglist,optlist = pst.configure.config_init(tel0)           
        _opts_list['arg'] = arglist
        _opts_list[tel0] = optlist

    # define email content
    _opts_list['arg']['email']['emailsub']+='[offline]'
    _opts_list['arg']['email']['emailcontent']='offline GW search\n'
    _opts_list['arg']['phone']['phonecontent']='offline GW search\n'
    _opts_list['arg']['email']['files'] = []
    _opts_list['arg']['email']['images'] = []

    # decide prioritization method        
    _trigger,_mass,_dist,_ngal = eval(_opts_list['arg']['priorization']['trigger']),\
                                 eval(_opts_list['arg']['priorization']['mass']),\
                                 eval(_opts_list['arg']['priorization']['dist']),\
                                 eval(_opts_list['arg']['priorization']['number'])
        
    # if man mode, role=test
    _opts_list['arg']['email']['role'] = 'test'

    if _trigger or _mass or _dist or _ngal:                       
        print('Trigger search')           
        logging.info('with fits:%s'%mapname)
        if _trigger and mapname is None:mapname = input('specify an input map name:')
        pst.pstdef.main(mapname, _opts_list, 'man')

    else:
        print('Normal search')
        sys.exit('TBD')
