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
configfile1,configfile2 = 'pstools.par','%s.par'

##### params option define
# must have sth
_args = ['telescope','dir','title','name','obstime','limmoon','filter']

# must be float,float
_argsc = ['limdist','limra','limdec','limmag','contours','showmap']

# must be float
_argsi = ['shift','nfields','limfields','repeat']

# must be float
_argsf = ['theta','phi','min','max','weight','lat', 'lon', 'alt', \
          'fovw', 'fovh', 'exptime', 'rottime','uarea','timelast',\
          'limsun','limalt','dither','limsbjr']

# must choose from options
_opt = {'trigger':['True','False'], 'galaxy':['0','1','2'], \
        'dist':['True','False'], 'test':['True','False'], \
        'cachemode':['1','2','3','4'], 'catalog':['1','2'], \
        'verbose':['True','False'], 'showmode':['1','2','3','4','5'], \
        'coord':['C','E','G'], 'norm':['hist','log','None'],\
        'activate':['True','False'],'scheduler':['T','G','A'],\
        'order':['1','2','3','4'], 'action':['1','2']}

def load_config(_tel='VST',_verbose=False, _log=False):

    # defined config list for general user and specific telescopes
    # search configure file in:
    # 1. current dir
    # 2. the PSTools default dir

    arglist,optlist = {},{}
    for configfile,paramfile,_class in zip([configfile1,configfile2%_tel],\
                                           [arglist,optlist],\
                                           ['general','telescope']):        
        if os.path.exists(_dir1+configfile): _configfile = _dir1+configfile            
        elif os.path.exists(_dir2+configfile): _configfile = _dir2+configfile
        else:  # if no config file, initialize one
            _configfile = config_init(_tel,_dir1,_class,_log,_verbose)
            if not _verbose: 
                print ('!!! Warning: no %s file of %s found, initialize a new one, as %s'%\
                       (_class,_tel,_configfile))
            sys.exit('Modify %s before running scripts'%_configfile)

        if not _configfile: return False,False

        # show
        if _verbose: print ('Read %s params for %s from %s \n'%(_class,_tel,_configfile))

        # read params from config list
        config = configparser.ConfigParser()
        config.read(_configfile)

        # logger
        if _log:
            import logging
            logging.info('%s-%s: %s'%(_tel,_class,_configfile))
       
        for s in config.sections():
            paramfile[s] = {}            
            for o in config.options(s):
                if o in _args and config.get(s,o) == '':
                    sys.exit('### Error: %s-%s need to be filled!'%(s,o))
                if o in _argsf:
                    try: float(config.get(s,o)) 
                    except:
                        sys.exit('### Error: %s-%s need to be a float number!'%(s,o))
                if o in _argsi:
                    try: int(config.get(s,o)) 
                    except:
                        sys.exit('### Error: %s-%s need to be an int number!'%(s,o))
                if o in _argsc:
                    try: 
                        paramfile[s][o] = [float(config.get(s,o).split(',')[0]),\
                                           float(config.get(s,o).split(',')[1])]
                    except:
                        sys.exit('### Error: %s-%s need to be float, float!'%(s,o))
                if o in _opt:
                    if not config.get(s,o) in _opt[o]:
                        sys.exit('### Error: %s-%s, wrong options!'%(s,o))
                paramfile[s][o] = config.get(s,o)
    return arglist,optlist

def config_init(_tel='VST',_dir='./',_class='telescope',_log=False,_verbose=False):
    # Initialize par
 
    _warning = '!!! Warning: no %s %s par file found both in %s'%(_tel,_class,_dir)
    if _verbose: print (_warning)
    if _log: 
        import logging
        logging.info(_warning)

    # define config file name
    if _class=='general':
        print ('### Error: pstools.par missing, download it '+\
               'in https://github.com/saberyoung/pstools/tree/master/src/pst/default/')
        return False
    elif _class=='telescope':        
        _tplf = _dir2+'temptel.par'
        if os.path.exists(_tplf):
            shutil.copyfile(_tplf, _dir+configfile2%_tel)
            _warning = '!!! Warning: a default file was generated (as %s)'%(_dir+configfile2%_tel)
            if _verbose: print (_warning)
            if _log: logging.info(_warning)
            return _dir+configfile2%_tel
        else:
            print ('Error: temptel.par missing, download it '+\
                   'in https://github.com/saberyoung/pstools/'+\
                   'tree/master/src/pst/default/')
            return False
    else: return False
