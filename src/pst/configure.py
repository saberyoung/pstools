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
dirlist = {'default':'%s/default/'%pstdir,\
           'local':'./'}

# define config files
parlist = {'general': 'pstools.par',\
           'telescope':'pst_%s.par'}

# define docs and default configs
filelist = {'general':{'doc':'%s/default/fmt/doc_gen.txt'%pstdir,\
                       'par':'%s/default/fmt/def_gen.par'%pstdir},\
            'telescope':{'doc':'%s/default/fmt/doc_tel.txt'%pstdir,\
                         'par':'%s/default/fmt/def_tel.par'%pstdir}}

# Initialize par file
def config_init(_tel='VST',_dir='./',_class='telescope',_verbose=False):

    # define config file name
    if _class in ['general', 'telescope']: 
        _tplf = filelist[_class]['par']
        _ofile = _dir + parlist[_class]
        if _class=='telescope': _ofile = _ofile%_tel
    else: return False

    # copy
    if os.path.exists(_tplf):
        shutil.copyfile(_tplf,_ofile)
        if _verbose: 
            print('\tinitialize %s'%_ofile)
        return _ofile
    else:
        print ('### Error: %s missing'%_tplf)
        return False    

def load_config(_tel='VST',_verbose=False,_inclass=False, _indir=False):

    # defined config list for general user and specific telescopes
    # search configure file in:
    # 1. current dir
    # 2. the PSTools default dir

    _olist = {}
    if _inclass: _cls = [_inclass]
    else: _cls = ['general','telescope']
    if _indir: _dirl = [_indir]
    else: _dirl = dirlist.values()

    for _class in _cls:

        _olist[_class]={}
        # define config name
        _cfile = parlist[_class]
        if _class == 'telescope': _cfile = _cfile%_tel

        _configfile = None
        # search dir
        for _dir in _dirl:
            if os.path.exists(_dir+_cfile):
                _configfile = _dir+_cfile

        # if no config file, initialize one in local
        if _configfile is None:
            if _verbose:
                print('!!! Warning: no %s file of %s found from %s'%(_class,_tel,_dirl))
            _configfile = config_init(_tel,dirlist['local'],_class,_verbose)
        _olist[_class]['config'] = _configfile

        # read healper from doc
        _docs = configparser.ConfigParser()
        _docs.read(filelist[_class]['doc'])

        # read params
        _paraml,_ckresult = {},True
        if _configfile:
            if _verbose: print ('\tRead %s params for %s from %s'%(_class,_tel,_configfile))

            # read params from config list
            config = configparser.ConfigParser()
            config.read(_configfile)

            for s in config.sections():
                _paraml[s]={}
                for o in config.options(s):
                    _var = config.get(s,o)
                    _paraml[s][o] = _var

                    # check format
                    _help,_fmt,_opt = _docs.get(s,o).replace('\n','').split('//')
                    if _fmt == 'arg':
                        if _var == '':
                            _ckresult=False
                            if _verbose:
                                print ('### Error %s: %s-%s, %s'%\
                                       (_configfile,s,o,_fmt))
                    elif _fmt == 'float':
                        try: 
                            float(_var) 
                        except:
                            _ckresult=False
                            if _verbose:
                                print('### Error %s: %s-%s, %s'%\
                                      (_configfile,s,o,_fmt))
                    elif _fmt == 'int':
                        try: 
                            int(_var) 
                        except:
                            _ckresult=False
                            if _verbose:
                                print('### Error: %s-%s-%s: %s'%\
                                      (_configfile,s,o,_fmt))
                    elif _fmt == 'comma':
                        try: 
                            for ii in _var.split(','):float(ii)
                        except:
                            _ckresult=False
                            if _verbose:
                                print('### Error %s: %s-%s, %s'%\
                                      (_configfile,s,o,_fmt))
                    elif _fmt == 'opt': 
                        pass
                    else:
                        if not _var in _fmt.split(','):
                            _ckresult = False
                            if _verbose:
                                print('### Error %s: %s-%s, %s'%\
                                      (_configfile,s,o,_fmt))
        _olist[_class]['params'] = _paraml
        _olist[_class]['check'] = _ckresult
    return _olist

def modify_config(_configfile,_indict):

    config = configparser.ConfigParser()
    config.read(_configfile)
    # update config
    _dict = {}
    for s in config.sections():
        _dict[s] = {}
        for o in config.options(s):
            # upchange val
            answ = config.get(s,o)
            # change var
            for _nn in range(len(_indict)):
                _indk1,_indk2 = _indict[_nn]
                if s in _indk1 and o in _indk1:
                    answ = _indk2
            _dict[s][o] = answ

    # generate config
    config = configparser.ConfigParser()
    for s in _dict:
        config[s] = {}
        for o in _dict[s]:
            config[s][o] = _dict[s][o]
    print ('Warning: rebuild %s'%_configfile)
    if os.path.exists(_configfile):os.remove(_configfile)
    with open(_configfile, 'w') as configfile:
        config.write(configfile)

def gen_config(_configfile,_class):

    # read doc
    _docs = configparser.ConfigParser()
    _docs.read(filelist[_class]['doc'])

    # input params
    config = configparser.ConfigParser()
    config.read(_configfile)
    _dict = {}
    for s in config.sections():
        try:_dict[s]
        except:_dict[s] = {}
        for o in config.options(s):
            _ok = False
            _help,_fmt,_opt = _docs.get(s,o).replace('\n','').split('//')
            if config.get(s,o) != '':answ1 = config.get(s,o)
            else:answ1 = _opt.replace(' ','')
            while not _ok:
                if _fmt == 'arg':
                    print ('\n','-'*20,'\n',_help,'\n')
                    answ = input('%s-%s, %s:\t%s\n:\t'%(s,o,_fmt,answ1))
                    if answ == '':answ=answ1
                    if answ != '': _ok = True
                elif _fmt == 'float':
                    print ('\n','-'*20,_help,'\n')
                    answ = input('%s-%s, %s:\t%s\n:\t'%(s,o,_fmt,answ1))
                    if answ == '':answ=answ1
                    try: 
                        float(answ) 
                        _ok = True
                    except: pass
                elif _fmt == 'int':
                    print ('\n','-'*20,_help,'\n')
                    answ = input('%s-%s, %s:\t%s\n:\t'%(s,o,_fmt,answ1))
                    if answ == '':answ=answ1
                    try: 
                        int(answ) 
                        _ok = True
                    except:pass
                elif _fmt == 'comma':
                    print ('\n','-'*20,_help,'\n')
                    answ = input('%s-%s, %s:\t%s\n:\t'%(s,o,_fmt,answ1))
                    if answ == '':answ=answ1
                    try: 
                        for ii in answ.split(','): float(ii)
                        _ok = True
                    except:pass
                elif _fmt == 'opt':
                    answ = input('%s-%s, %s:\t%s\n:\t'%(s,o,_fmt,answ1))
                    if answ == '':answ=answ1
                    _ok = True
                else: 
                    print ('\n','-'*20,_help,'\n')
                    answ = input('%s-%s, %s:\t%s\n:\t'%(s,o,_fmt,answ1))
                    if answ == '':answ=answ1
                    if answ in _fmt.split(','):_ok = True
            _dict[s][o] = answ

    # generate
    config = configparser.ConfigParser()
    for s in _dict:
        config[s] = {}
        for o in _dict[s]:
            config[s][o] = _dict[s][o]
    print ('Warning: rebuild %s'%_configfile)
    if os.path.exists(_configfile):os.remove(_configfile)
    with open(_configfile, 'w') as configfile:
        config.write(configfile)

def choose(_tel,_class,_inconfig):

    # init dict and cfile
    _dicti = _inconfig[_class]['params']
    _cfile = _inconfig[_class]['config']
    _dir = '%s/'%os.path.dirname(_cfile)

    # read doc
    _docs = configparser.ConfigParser()
    _docs.read(filelist[_class]['doc'])

    # choose
    done, _dict, _klist = False, _dicti, ['basic']
    print ('-'*20+'\nPSTools, revise configure for %s\n\n\n'%_cfile)
    while not done:
        _keyl = ''
        for zz in _dict:_keyl+='%s   '%zz
        answ = input('%s   [a]ll   [b]ye   [o]ut\n:\t'%_keyl)
        if answ in ['B','b','Bye','bye','BYE']: done=True
        elif answ in ['O','o','OUT','out','Out']:
            if len(_klist)>1:
                del _klist[-1]
                print ('\n','>'*5,'go to %s\n'%_klist[-1])
                _dict = _dicti
                for _key in _klist[1:]:_dict = _dict[_key]
            else: print ('\nAlready at basic level...\n')
        elif answ in _dict.keys():
            print ('\n','>'*5,'go to %s\n'%answ)
            _klist.append(answ)
            _dict = _dict[answ]
        elif answ in ['A','a','ALL','all','All']:
            for _key in _dict:print ('\n\t--> %s : %s \n'%(_key,_dict[_key]))
        else: print ('\n!!! Error: wrong input...\n')

        try:
            # dict level
            _dict.keys()
        except:
            # val level
            answ1 = input('\t%s=%s, [h]elp (will not modify) or any input to modify it?\n:\t'%\
                          (_klist[-1],_dict))
            if answ1 in ['h','H']: # help
                for s in _docs.sections():
                    for o in _docs.options(s):
                        if o == _klist[-1] and s == _klist[-2]:
                            _help,_fmt,_opt = _docs.get(s,o).\
                                              replace('\n','').split('//')
                            print ('\n','-'*10,'\n- Help: %s\n- type:%s\n'%(_help, _fmt))
                #refine
                del _klist[-1]
                _dict = _dicti
                for _key in _klist[1:]: _dict = _dict[_key]
            else:
                del _klist[0]
                pst.modify_config(_cfile,[(_klist, answ1)])
                _config = pst.load_config(_tel,_inclass=_class, _indir=_dir)
                _dict = _config[_class]['params']         
                _klist1 = ['basic']
                del _klist[-1]
                for _key in _klist:
                    _dict = _dict[_key]
                    _klist1.append(_key)
                _klist = _klist1
