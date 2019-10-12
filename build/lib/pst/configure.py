"""############################################################################ 
2019/1/30 Start
Configure file generator
""" ############################################################################
from __future__ import print_function
from builtins import input

import os,sys,pst,configparser,shutil
from astropy import units as u
from astropy.coordinates import SkyCoord

pstdir = pst.__path__[0]

# define dirs
_dir1, _dir2 = './', '%s/default/'%pstdir

# define config files
configfile1,configfile2 = 'pstools.par','pst_%s.par'
doc1,doc2 = _dir2+'pstools.txt',_dir2+'pst_tel.txt'

##### params option define
# must have sth
_args = ['telescope','dir','name','obstime',\
         'limmoon','filter']

# must be number, number
_argsc = ['limdist','limra','limdec','limmag','contours','showmap','ob']

# must be int number
_argsi = ['shift','nfields','repeat']

# must be float number
_argsf = ['weight','lat', 'lon', 'alt', 'fovw', \
          'fovh', 'exptime', 'rottime','uarea','timelast',\
          'limsun','limalt','dither']

# must choose from options
_opt = {'trigger':['True','False'], 'galaxy':['0','1','2'], \
        'dist':['True','False'], 'test':['True','False'], \
        'cachemode':['1','2','3','4'], 'catalog':['1','2'], \
        'verbose':['True','False'], 'showmode':['1','2','3','4','5'], \
        'activate':['True','False'],'scheduler':['T','G','A'],\
        'order':['1','2','3','4'], 'action':['1','2']}    

def load_config(_tel='VST',_verbose=False, \
                _check=False, _cfile=False):

    # defined config list for general user and specific telescopes
    # search configure file in:
    # 1. current dir
    # 2. the PSTools default dir

    arglist,optlist,configfiles = {},{},[]
    for configfile,paramfile,_class in zip([configfile1,configfile2%_tel],\
                                           [arglist,optlist],\
                                           ['general','telescope']):        
        if os.path.exists(_dir1+configfile): _configfile = _dir1+configfile            
        elif os.path.exists(_dir2+configfile): _configfile = _dir2+configfile
        else:  # if no config file, initialize one
            _configfile = config_init(_tel,_dir1,_class,_verbose)
            print('!!! Warning: no %s file of %s found, initialize %s'%\
                  (_class,_tel,_configfile))

        if _configfile:
            # show
            if _verbose: print ('\tRead %s params for %s from %s'%(_class,_tel,_configfile))

            # read params from config list
            config = configparser.ConfigParser()
            config.read(_configfile)
       
            for s in config.sections():
                paramfile[s] = {}            
                for o in config.options(s):
                    if _check:
                        if o in _args and config.get(s,o) == '':
                            sys.exit('### Error: %s-%s-%s need to be filled'%\
                                     (_configfile,s,o))
                        if o in _argsf:
                            try: float(config.get(s,o)) 
                            except:
                                sys.exit('### Error: %s-%s-%s need to be a float number!'%\
                                     (_configfile,s,o))
                        if o in _argsi:
                            try: int(config.get(s,o)) 
                            except:
                                sys.exit('### Error: %s-%s-%s need to be an int number!'%\
                                     (_configfile,s,o))
                        if o in _argsc:
                            try: 
                                paramfile[s][o] = [float(config.get(s,o).split(',')[0]),\
                                                   float(config.get(s,o).split(',')[1])]
                            except:
                                sys.exit('### Error: %s-%s-%s need to be float, float!'%\
                                         (_configfile,s,o))
                        if o in _opt:
                            if not config.get(s,o) in _opt[o]:
                                sys.exit('### Error: %s-%s-%s, wrong options!'%\
                                         (_configfile,s,o))
                    paramfile[s][o] = config.get(s,o)
            configfiles.append(configfile)
    if _cfile:return configfiles,arglist,optlist
    else: return arglist,optlist

def config_init(_tel='VST',_dir='./',_class='telescope',_verbose=False):
    # Initialize par
 
    # define config file name
    if _class=='general':
        _tplf = _dir2+'pstools_def.par'
        if os.path.exists(_tplf):
            shutil.copyfile(_tplf, _dir+configfile1)
            return _dir+configfile1
        else:
            print ('### Error: pstools.par missing, download it '+\
                   'in https://github.com/saberyoung/pstools/tree/'+\
                   'master/src/pst/default/')
            return False

    elif _class=='telescope':        
        _tplf = _dir2+'pst_tel.par'
        if os.path.exists(_tplf):
            shutil.copyfile(_tplf, _dir+configfile2%_tel)
            return _dir+configfile2%_tel
        else:
            print ('!!! Error: temptel.par missing, download it '+\
                   'in https://github.com/saberyoung/pstools/'+\
                   'tree/master/src/pst/default/')
            return False
    else: return False

def modify_config(_configfile,_indict):

    config = configparser.ConfigParser()
    config.read(_configfile)
    _dict = {}
    for s in config.sections():
        try:_dict[s]
        except:_dict[s] = {}
        for o in config.options(s):
            if o in _indict: answ = _indict[o]
            else: answ = config.get(s,o)
            _dict[s][o] = answ
    # generate
    config = configparser.ConfigParser()
    for s in _dict:
        config[s] = {}
        for o in _dict[s]:
            config[s][o] = _dict[s][o]
    if os.path.exists(_configfile):os.remove(_configfile)
    with open(_configfile, 'w') as configfile:
        config.write(configfile)

def gen_config(_configfile):

    if configfile1 in _configfile:_doc = doc1
    else: _doc=doc2
    _docs = configparser.ConfigParser()
    _docs.read(_doc)

    # input params
    config = configparser.ConfigParser()
    config.read(_configfile)
    _dict = {}
    for s in config.sections():
        try:_dict[s]
        except:_dict[s] = {}
        for o in config.options(s):
            _ok = False
            answ = config.get(s,o)
            while not _ok:
                if o in _args and answ == '':
                    print ('\n','-'*20,_docs.get(s,o),'\n')
                    answ = input('%s-%s, TB filled:\t'%(s,o))
                elif o in _argsf:
                    try: 
                        float(answ) 
                        _ok = True
                    except:
                        print ('\n','-'*20,_docs.get(s,o),'\n')
                        answ = input('%s-%s, float:\t'%(s,o))
                elif o in _argsi:
                    try: 
                        int(answ) 
                        _ok = True
                    except:
                        print ('\n','-'*20,_docs.get(s,o),'\n')
                        answ = input('%s-%s, int:\t'%(s,o))
                elif o in _argsc:
                    try: 
                        float(answ.split(',')[0])
                        float(answ.split(',')[1])
                        _ok = True
                    except:
                        print ('\n','-'*20,_docs.get(s,o),'\n')
                        answ = input('%s-%s, "float,float":\t'%(s,o))
                elif o in _opt:
                    if not answ in _opt[o]:
                        print ('\n','-'*20,_docs.get(s,o),'\n')
                        answ = input('%s-%s, options: %s:\t'%(s,o,_opt[o]))
                    else: _ok = True
                else: 
                    _ok = True
            _dict[s][o] = answ

    # generate
    config = configparser.ConfigParser()
    for s in _dict:
        config[s] = {}
        for o in _dict[s]:
            config[s][o] = _dict[s][o]
    if os.path.exists(_configfile):os.remove(_configfile)
    with open(_configfile, 'w') as configfile:
        config.write(configfile)
