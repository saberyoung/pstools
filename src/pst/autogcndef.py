"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input
import os,glob,gcn,gcn.handlers,gcn.notice_types,scipy.stats,logging,voeventparse,astropy.time
import healpy as hp
import numpy as np
import math as mt
import pst

#################################################
# Function to call every time a GCN is received.
# Run only for notices of type LVC_INITIAL or LVC_UPDATE.

@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,            # pointing and send alerts
    gcn.notice_types.LVC_INITIAL,                
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):  

    # judge if new voevet, no? return
    _voname = root.attrib['ivorn']+'.xml'
    if os.path.exists(os.path.basename(_voname)):return
    else:pass

    # if new voevent, alert, store voevet
    print('###Alert')
    logging.info('###Alert %s'%_voname)
    with open(os.path.basename(_voname),'w') as _vo:_vo.write(payload)    

    # define parameter dictionary
    _opts_list = {}

    # decide which telescopes for follow up
    # in the "current" directory: pstools.default, pst_tel.default
    for _pstdir in ['./']:#, '%s/default/'%pst.__path__[0]]:
        _pstfile = glob.glob(_pstdir + 'pst_*.default')        
        if len(_pstfile)==0:
            logging.info('No default files found in %s'%_pstdir)
            continue
        for ff in _pstfile:
            _tool,tel0 = os.path.splitext(os.path.basename(ff))[0].split('_')
            if tel0 in _opts_list:continue
            arglist,optlist = pst.configure.config_init(tel0)           
            _opts_list['arg'] = arglist
            _opts_list[tel0] = optlist       
            _info = 'Use %s in %s for %s auto search'%(os.path.basename(ff),_pstdir,tel0)
            print(_info)
            logging.info(_info)

    # if arg default file missing
    if not 'arg' in _opts_list:
        print('!!! Error: take care, trigger is coming, however, pstools.default missing in pst directory!!!')
        logging.info('!!! Error: take care, trigger is coming, however, pstools.default missing in pst directory!!!')

    # lenght of default file >=2: one for arg, at least one for telescope
    if len(_opts_list) < 2:
        print('!!! Error: take care, trigger is coming, however, no default file found in pst directory!!!')    
        logging.info('!!! Error: take care, trigger is coming, however, no default file found in pst directory!!!')

    # transfer further informations via _opts_list dict
    # -> files, content, images
    _opts_list['arg']['email']['files'] = [os.path.basename(_voname)]
    _opts_list['arg']['email']['images'] = []

    # read GW parameters from voevent, record them to _opts_list dict
    if 'LVC' in root.attrib['ivorn']:
        _opts_list['arg']['email']['emailcontent'] = '#online LVC alert:\n'
        _opts_list['arg']['phone']['phonecontent'] = 'online LVC alert, '
    else:
        try:_tname = os.path.basename(root.attrib['ivorn']).split('#')[0]
        except:
            # related to gcn.handlers.include_notice_types           
            logging.info('Check voevent %s!'%_voname)
            print('Check voevent %s!'%_voname)
            return
        _opts_list['arg']['email']['emailcontent'] = '#online %s search:\n'%_tname
        _opts_list['arg']['phone']['phonecontent'] = 'online %s alert, '%_tname

    # check role: test or obs
    _opts_list['arg']['email']['emailcontent'] += '#\t%s alert\n'%root.attrib['role']
    _opts_list['arg']['phone']['phonecontent'] += '%s: '%root.attrib['role']

    # record role in email list
    if root.attrib['role'] == 'test':
        if not eval(_opts_list['arg']['priorization']['test']):return
    _opts_list['arg']['email']['role'] = root.attrib['role']

    if not root.attrib['role'] in ['test','observation']:
        # for LVC only two roles
        # take care for the other alert type: neutrino, GRB, etc
        logging.info('New role?')
        print('New role?')
        #return

    # for the voevent without healpy map, go main2
    if any(x in _voname for x in ['LVC']):go_main = 1
    # for alert with healpy map, go main1
    else:go_main = 2

    # for LVC preliminary initial update cases:
    if go_main == 1:
        # check stage: Preliminary or Initial or Update
        if 'Preliminary' in _voname:
            _opts_list['arg']['email']['emailcontent'] += '#\tPreliminary announcement\n'
            _opts_list['arg']['phone']['phonecontent'] += 'preliminary announcement for '
        elif 'Initial' in _voname:
            _opts_list['arg']['email']['emailcontent'] = '#\tInitial announcement\n'
            _opts_list['arg']['phone']['phonecontent'] += 'initial announcement for '
        elif 'Update' in _voname:
            _opts_list['arg']['email']['emailcontent'] = '#\tUpdate announcement\n'
            _opts_list['arg']['phone']['phonecontent'] += 'update announcement for '
        else:
            # impossible to arrvie here!
            logging.info('New stage?')    
            return
    # for neutrino, GRB, etc
    else:
        print(_voname)
        logging.info(_voname)

    # Read all of the VOEvent parameters from the "What" section.
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # store voevent
    _opts_list['arg']['voevent']=params

    # main process
    if go_main == 1: 

        # download map and read mapname    
        mapurl = params['skymap_fits']
        mapname = pst.pstdef.get_skymap(mapurl,os.path.basename(root.attrib['ivorn']),_opts_list['arg']['data']['dir'])
        _opts_list['arg']['email']['files'].append(mapname)

        logging.info('Finished checking VOevent, now starting main process for %s'%mapname)        
    
    else: # for no map case
        with open(os.path.basename(_voname), 'rb') as f:

            v = voeventparse.load(f)

            try:
                _nra,_ndec,_loc = float(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Value2.C1),\
                                  float(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Value2.C2),\
                                  float(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Error2Radius)
            except:
                logging.info('Error: no ra,dec found in voevent!!!')
                return

            if _nra == 0 and _ndec == 0:
                logging.info('Error: no ra,dec reported in voevent!!!')
                return
                
            try:
                _timeobs = str(v.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Time.TimeInstant.ISOTime)
                #_importance = next(root.iterfind('.//Why')).attrib['importance']
            except:
                logging.info('Error: no timeobs found in voevent!!!')
                return

            _timeobs1 = astropy.time.Time(_timeobs, format='isot', scale='utc')
            _mjdobs = _timeobs1.mjd

            try:_object = v.Why.Inference.Name
            except:_object = 'unKnown'           

        # create healpix fits map and read mapname
        mapname = '%s%s.fits'%(_opts_list['arg']['data']['dir'],\
                          os.path.basename(root.attrib['ivorn']))
        nside = int(_opts_list['arg']['priorization']["nside"])

        # circle with ra,dec,radius=err or max(fov)
        if _loc>0: _radius = _loc
        else:
            fovl = []
            for _ntt,_tt in enumerate(_opts_list):
                if _tt == 'arg':continue
                fovl.append(float(_opts_list[_tt]['telescope']['fovw']))
                fovl.append(float(_opts_list[_tt]['telescope']['fovh']))
            _radius = max(fovl)

        gradius = _radius*2*mt.pi/360 # from deg to radians
        _pmap = np.zeros(hp.nside2npix(nside))
        _index = pst.DeclRaToIndex(_ndec,_nra,nside)
        _pmap[_index]+=1
        _pmap=hp.sphtfunc.smoothing(_pmap,fwhm=gradius)
        
        _pmap = _pmap/sum(_pmap)       
        hlist = [('CREATOR','PSTOOLS'),
                 ('OBJECT',_object),
                 ('NSIDE',nside),
                 ('MJD-OBS',_mjdobs),
                 ('DATE-OBS',_timeobs)]

        hp.write_map(mapname,_pmap, coord='C', extra_header=hlist, overwrite=True)
        _opts_list['arg']['email']['files'].append(mapname)
        logging.info('Finished checking VOevent, now starting main 2 process for ra=%.2f dec=%.2f with error=%.2f'%\
                     (_nra,_ndec,_loc))

    pst.pstdef.main(mapname, _opts_list, 'auto')
