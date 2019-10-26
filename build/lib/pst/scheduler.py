from __future__ import print_function
from builtins import input
import pprint
from astropy import units as u
from astropy.coordinates import SkyCoord
import re
import astropy.time

_skytlist = {'1':"Photometric",\
             '2':"Clear",\
             '3':"Variable, thin cirrus",\
             '4':"Variable, thick cirrus"}

def VST(_list,triggerid,enviroment,username,password,containerid,_userpriority,_airmass, _skyTransparency, _fli,_seeing,_exp,_dither,_nexp):

    # generate OB list
    _oblist = {}   
    for ii in _list:       
        if len(ii)==0:continue
        if ii[0]=='#':continue

        try:
            ra,dec,filtro,nob = float(ii.split()[0]),float(ii.split()[1]),\
                                ii.split()[2],int(ii.split()[-1].split('-')[0])
        except:
            continue

        _time = astropy.time.Time('T'.join(re.findall('\d\d\d\d-\d\d-\d\d',ii) + \
                        re.findall('\d\d:\d\d:\d\d',ii)),format='isot',out_subfmt='date_hm')

        ''' TBD !!! '''
        _diffhour = astropy.time.TimeDelta(3600, format='sec') 

        _time1 = str(_time-_diffhour)
        _time2 = str(_time+_diffhour)

        try:_oblist[nob]
        except:_oblist[nob]=[]
        _oblist[nob].append([ra,dec,filtro,_time1,_time2])

    # start
    from pst import p2api
    p = pprint.PrettyPrinter(indent=4)
    
    # login
    api = p2api.ApiConnection(enviroment,username,password)
    print('*** logged in')

    # create folder for general
    runContainerId = containerid
    folder, folderVersion = api.createFolder(runContainerId, triggerid)
    folderId0 = folder['containerId']
    print('*** created folder')
    p.pprint(folder) 

    for _nob in _oblist:
        # create folder for OBs        
        folder, folderVersion = api.createConcatenation(folderId0, 'Concatenation %i'%_nob)
        folderId1 = folder['containerId']
        print('*** created Concatenation %i'%_nob)
        p.pprint(folder)

        _npointing = 0
        for _ra,_dec,_filter,_time1,_time2 in _oblist[_nob]:
            _npointing+=1

            # create OB
            ob, obVersion = api.createOB(folderId1, 'TOO_OMEGACAM_GW_p%i_e1'%(_npointing))
            obId = ob['obId']
            print('*** created OB')
            p.pprint(ob)

            # ra dec from deg to hms
            _p = SkyCoord(ra=_ra*u.degree, dec=_dec*u.degree)                                  

            _rahms = '%.2i:%.2i:%.3f'%(_p.ra.hms[0],abs(_p.ra.hms[1]),abs(_p.ra.hms[2]))
            _dechms = '%.2i:%.2i:%.3f'%(_p.dec.dms[0],abs(_p.dec.dms[1]),abs(_p.dec.dms[2]))

            # edit OB
            ob['userPriority']                   = _userpriority
            ob['target']['name']                 = triggerid
            ob['target']['ra']                   = _rahms
            ob['target']['dec']                  = _dechms
            ob['target']['properMotionRa']       = 0
            ob['target']['properMotionDec']      = 0
            ob['constraints']['name']            = 'No name'
            ob['constraints']['airmass']         = _airmass
            ob['constraints']['skyTransparency'] = _skytlist[_skyTransparency]
            ob['constraints']['fli']             = _fli
            ob['constraints']['seeing']          = _seeing
            ob, obVersion = api.saveOB(ob, obVersion)
            print('*** saved OB changes')

            # attach Acquisition Template
            acqTpl, acqTplVersion = api.createTemplate(obId, "OMEGACAM_img_acq")
            print('*** attached Acquisition Template')
            p.pprint(acqTpl)

            if _npointing == 1:_ias = True
            else:_ias = False
            # edit Acquisition Template
            acqTpl, acqTplVersion  = api.setTemplateParams(obId, acqTpl, {               
                'TEL.ROT.OFFANGLE'    :0,
                'TEL.GS1.ALPHA'       :'0',
                'TEL.GS1.DELTA'       :'0',
                'TEL.GS1.MAG'         :12,
                'TEL.GS2.ALPHA'       :'0',
                'TEL.GS2.DELTA'       :'0',
                'TEL.GS2.MAG'         :12,
                'TEL.ADC.TYPE'        :"NONE",
                'INS.FILT.NAME'       :_filter,
                'OCS.AG.START'        :False,
                'OCS.IA.START'        :_ias
            }, acqTplVersion)
            print('*** saved Acquisition Template changes')

            # attach Science Template
            scTpl, scTplVersion = api.createTemplate(obId, "OMEGACAM_img_obs_dither")
            print('*** attached Science Template')
            p.pprint(scTpl)

            # edit Science Template
            scTpl, scTplVersion = api.setTemplateParams(obId, scTpl, {
                'DET1.WIN1.UIT1'           :_exp,
                'SEQ.NEXPO'                :_nexp,
                'TEL.TARG.OFFSETSIZEX'     :25, # fixed dither
                'TEL.TARG.OFFSETSIZEY'     :85, # fixed dither
                'TEL.TARG.DX'              :8,
                'TEL.TARG.DY'              :8,
                'TEL.TARG.PATTERN'         :"diag",
                'INS.FILT.NAME'            :_filter,
                'OCS.STRTG'                :"Mosaic"
            }, scTplVersion)
            print('*** saved Science Template changes')

            ''' no time window for too
            # define Absolute Time Constraints
            absTCs, atcVersion = api.getAbsoluteTimeConstraints(obId)
            absTCs, atcVersion = api.saveAbsoluteTimeConstraints(obId,[
                {
                    'from': _time1,
                    'to': _time2
                }
            ], atcVersion)
            print('*** added Absolute Time Constraints')
            p.pprint(absTCs)
            '''

            # Verify OB to status (D)
            response, _ = api.verifyOB(obId, True)
            if response['observable']:
                print('*** Congratulations. Your OB' , obId, ob['name'], 'is observable!')
            else:
                print('OB', obId, 'is >>not observable<<. See messages below.')
            print(' ', '\n  '.join(response['messages']))

    '''
    # submit to ESO
    runId = 6092526
    response, _ = api.submitRun(runId)
    print('*** submitted Observing Run to ESO for Review')
    p.pprint(response)
    print('*** Goodbye ***')
    '''
    return 'https://www.eso.org/p2demo/home/container/%i'%folderId0

def schmidt(_list,_scheduler,nexp,exp,dit,_id):

    _point,_autof,_exp,_exp1 = "BLOCK\n"+\
                               "GSTOP\n"+\
                               "P A%.6f D%.6f E2000.000000\n"+\
                               "PSTART\n",\
                               "CDIS\n"+\
                               'C L00 T0 R1 E1 F5 S0 D0.0 I "AutoFocus"\n'+\
                               "CBIN 3\n"+\
                               "CSTART\n"+\
                               "CFOCUS\n",\
                               "BLOCK\n"+\
                               "CDIS\n"+\
                               "%s"+\
                               "CBIN 1\n"+\
                               "GSTART\n"+\
                               "CSTART\n",\
                               'C L%02d T%i R%i E%.1f F%i S%i D%.1f I"%s"\n'

    #type of observation
    _too = {'obs':0}

    #type of filter
    _tof = {'None':0,'u':1,'B':2,'g':3,'V':4,'r':5,'i':6}

    ral,decl=[],[]
    _commandw=''
    _autofocus=True # auto once when targets are close
    _oo = open(_scheduler,'w')

    for ii in _list:
        if len(ii)==0:continue
        if ii[0]=='#':continue
        try:
            ra,dec,filtro=float(ii.split()[0]),\
                float(ii.split()[1]),ii.split()[2]     
        except:
            continue

        # pointing
        _commandw += _point%(ra,dec)

        while _autofocus:
            _commandw += _autof
            _autofocus=False

        # exposure
        _expw=''
        for _n,_f in enumerate(filtro.split(',')):
            _save = 1 # save pic
            _expw += _exp1%(_n,_too['obs'],nexp,exp,_tof[_f],_save,dit,_id)
        _expf=_exp%_expw
        
        _commandw+=_expf               
        _oo.write(_commandw)                       
        _commandw=''
    _oo.close()
