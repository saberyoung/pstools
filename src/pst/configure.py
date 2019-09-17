"""############################################################################ 
2019/1/30 Start
A testing file
""" ############################################################################
from __future__ import print_function
from builtins import input

import os,sys,pst
import configparser
from astropy import units as u
from astropy.coordinates import SkyCoord
import logging

# http://www.eso.org/~ndelmott/obs_sites.html
tel_list = {'schmidt':{'site':'Asiago','lat':45.8617,'lon':11.5283,'alt':1045,'fovw':1,'fovh':1,'exptime':45,'rottime':45,'limra':'0,360','limdec':'-90,90','ob':'3*3','stretegy':'T','resolution':0.87},
            'VST':{'site':'Cerro Paranal','lat':-24.625,'lon':-70.4033,'alt':2635,'fovw':1,'fovh':1,'exptime':40,'rottime':60,'limra':'0,360','limdec':'-90,90','ob':'3*3','stretegy':'T','resolution':0.21},            
            'CI':{'site':'Campo Imperatore','lon':13.5877,'lat':42.44194,'alt':2135,'fovw':1,'fovh':1,'exptime':60,'rottime':60,'limra':'0,360','limdec':'-90,90','ob':'3*3','stretegy':'T','resolution':0.87},
            'REM':{'site':'La Silla','lat':-29.25,'lon':-70.73,'alt':2375,'fovw':10/60.,'fovh':10/60.,'exptime':60,'rottime':60,'limra':'0,360','limdec':'-90,90','ob':'1*1','stretegy':'G','resolution':0.87},
            'sv':{'site':'Savelli','lat':39.3126,'lon':16.7784,'alt':2635,'fovw':21.1/60.,'fovh':21.1/60.,'exptime':90,'rottime':60,'limra':'0,360','limdec':'-90,90','ob':'1*1','stretegy':'G','resolution':0.87}}

def config_init(_tel='VST',_dir1='./',_dir2='%s/default/'%pst.__path__[0]):

    # defined config list for general user and specific user
    configfile1,configfile2 = 'pstools.default','pst_%s.default'%_tel
    arglist,optlist = {},{}
    for configfile,paramfile,_class in zip([configfile1,configfile2],\
                                           [arglist,optlist],\
                                           ['arglist','optlist']):        
        if os.path.exists(_dir1+configfile):_configfile = _dir1+configfile
        elif os.path.exists(_dir2+configfile):_configfile = _dir2+configfile
        else:
            # if not exists, initialize config list
            _configfileo = config_list(_tel,_dir2,_class)
            sys.exit('initialize a default file as %s'%_configfileo)

        # read from config list
        config = configparser.ConfigParser()
        config.read(_configfile)
        logging.info('%s: %s'%(_tel,_configfile))
        for s in config.sections():
            paramfile[s] = {}
            for o in config.options(s):
                paramfile[s][o] = config.get(s,o)        
    return arglist,optlist

def config_list(_tel='VST',_dir='./',_class='arglist'):

    if _tel in tel_list.keys():
        lat,lon,alt,fovw,fovh,exptime,rottime,limra,limdec,ob,strg,_res = tel_list[_tel]['lat'],\
                                                                          tel_list[_tel]['lon'],\
                                                                          tel_list[_tel]['alt'],\
                                                                          tel_list[_tel]['fovw'],\
                                                                          tel_list[_tel]['fovh'],\
                                                                          tel_list[_tel]['exptime'],\
                                                                          tel_list[_tel]['rottime'],\
                                                                          tel_list[_tel]['limra'],\
                                                                          tel_list[_tel]['limdec'],\
                                                                          tel_list[_tel]['ob'],\
                                                                          tel_list[_tel]['stretegy'],\
                                                                          tel_list[_tel]['resolution']
    else:
        print('!!! Warning: %s not in tel_list\nModify defualt file first before running codes in pst_%s.default'%(_tel,_tel))
        lat,lon,alt,fovw,fovh,exptime,rottime,limra,limdec,ob,strg,_res = 0,0,0,1,1,60,60,'0,360','-90,90','1*1','T',1.

    # filename
    if _class=='arglist':configfile = _dir+'pstools.default'
    elif _class=='optlist':configfile = _dir+'pst_%s.default'%_tel
    else:sys.exit('Error: Wrong params option!')
    _log = open(configfile,'w')

    # arglist for general users
    if _class=='arglist':
        _log.write(
        "[priorization] ;  ------------- how to generate priorization \n"
        "trigger:\tTrue\n"
        "\t\t; include trigger prob into score\n"
        "mass:\tTrue\n"
        "\t\t; include galaxy mass into score\n"
        "dist:\tFalse\n"
        "\t\t; include galaxy distance into score\n"
        "number:\tFalse\n"
        "\t\t; include galaxy number into score\n"
        "nside:\t1024\n"
        "\t\t; nside of priorization map\n"
        "\t\t; will change according to trigger healpix map if given and different from such nside\n"
        "test:\tFalse\n"
        "\t\t; Respond also to LVC test alerts?\n"
        "\n"
        "[data]  ;  ------------- define directory for data and some specific file name\n"
        "dir:\t./\n"
        "\t\t; data directory, for putting data files\n"
        "schfile:\tall.txt\n"
        "\t\t; output scheduler file format\n"            
        "\n"
        "[observe] ;------------- observation estimate\n"
        "obstime:\tnow\n"
        "\t\t; start time: [1] now;\n"
        "\t\t;             [2] time in utc, e.g.1999-01-01T00:00:00.1234\n"
        "\t\t;             [3] +/- x minutes, e.g. -30, 30 minutes before\n"
        "timelast:\t10\n"
        "\t\t; last for how much longer\n"
        "\t\t; if it's very long, the code stopped when no visible targets\n"
        "\t\t; Otherwise, sopped when timelast finished\n"
        "order:\t1\n"
        "\t\t; [1]greedy algorithm: rank with priority \n"
        "\t\t; [2]conservative algorithm: rank from west to east \n"       
        "\t\t; [3]optimal algorithm: consider slewing angle, start from westest pointing \n"
        "\t\t; [4]optimal algorithm: consider slewing angle, start from highest ranking \n"               
        "cov:\t\n"
        "\t\t; observing coverage \n"      
            "\t\t; if a float num (0-1), e.g. 0.9, given, force cover 90% trigger prob, otherwise, auto set\n"
        "limsun:\t-12\n"
        "\t\t; limitation of <limsun\n"
        "limmoon:\tauto\n"
        "\t\t; limit of moon (in deg)\n"
        "\t\t; e.g. 30 or auto\n"
        "\t\t; if set: pointings within such radius will be removed\n"
        "\t\t; if auto: do removement accoring to lunar phase\n"        
        "limsolorobject:\tjupiter,venus\n"
        "\t\t; 'mercury','venus','moon','mars','jupiter','saturn','uranus','neptune'\n"
        "limsolorobjectnum:\t2\n"
        "\t\t; pointings within such radius (deg) will be removed\n"
        "limalt:\t2.0\n"
        "\t\t; limitation on airmass\n"
        "threshold:\t\n"
        "\t\t;pixel threshold (e.g. 1e-3):\n"
        "\t\t;remain fields only when their corresponding probs is above the threshold\n"
        "\t\t;leave blanket, I will automatically do thresholding based on map size and FoV of telescope\n"
        "template:\t\n"
        "\t\t; search if available templates\n"
        "\t\t; e.g. VST,DES,PS,SDSS\n"
        "\n"       
        "[galaxies] ;  ------------- galaxy selection\n"
        "cachefile:\t\n"
        "\t\t; cache mode: store galaxies into $dir/cache, only if the number of galaxies is large\n"
        "catalog:\tGLADE\n"
        "\t\t; galaxy cat: GLADE/GWGC/2MASS \n"
        "filter:\tB\n"
        "\t\t; filter for galaxy absolute magniture \n"
        "size:\t-1\n"
        "\t\t; galaxy catalog size number, -1 for full size\n"        
        "dist:\t0,9999\n"
        "\t\t; galaxy distance range\n"
        "\t\t; if trigger (CBC GW, or?) distance available, auto select galaxy range; otherwise, applied settings here\n"
        "\t\t; in principle, one should select all galaxies, and then wait for a 3D prob, in case, some one wants to concentrate only local volume\n"
        "\t\t; Notice: GWGC is only valid within 100 Mpc, GLADE will go up to 300Mpc however, not complete\n"
        "limra:\t0,360\n"
        "\t\t; Ra range (in deg) for galaxy selection\n"
        "limdec:\t-90,90\n"
        "\t\t; Dec range (in deg) for galaxy selection\n"
        "limmag:\t-18\n"
        "\t\t; absolute limiting magnitude\n"
        "\t\t; remian galaxies only when they are brighter than this number\n"
        "limebv:\tFalse\n"
        "\t\t; ebv cut\n"
        "\t\t; if not False: remove pointings located in high ebv (e.g. $limebv<2.5) area \n"
        "radius:\t5\n"
        "\t\t; FWHM for galaxy smoothing (in arcmin)\n"
        "\t\t; False: no smoothing\n"
        "\n"         
        "[plot] ;  ------------- plotting parameters\n"
        "verbose:\tFalse\n"
        "\t\t; verbose mode: show detialed infors, plots, etc\n"
        "showmap:\ttrigger,cum\n"
        "\t\t; if verbose, which plot to show: trigger,galaxy,cum\n"
        "interactive:\tFalse\n"
        "\t\t; Once verbose, interactive mode will enable user to set plots \n"                  
        "rot_theta:\t0\n"
        "\t\t;default paramater for healpix view\n"
        "rot_phi:\t0\n"
        "\t\t;default paramater for healpix view\n"
        "ordering:\tFalse\n"
        "\t\t;default paramater for healpix view\n"
        "coord:\tC\n"
        "\t\t;default paramater for healpix view\n"
        "\t\t; default: Equatorial, don't change\n"
        "norm:\thist\n"
        "\t\t;default paramater for healpix view\n" 
        "[email] ;  ------------- email setting\n"        
        "sendemail:\tFalse\n"
        "\t\t;send email or not\n"
        "email:\t\n"       
        "\t\t;if yes, set email sender account\n"
        "emailpass:\t\n"        
        "\t\t;if yes, set email sender password\n"
        "emailsmtp:\t\n"        
        "\t\t;if yes, set email sender host and port\n"
        "emailto:\t\n"       
        "\t\t;if yes, set email receiver account\n"
        "emailsub:\tGW alerts:\n"
        "\t\t;if yes, set email subject format\n"
        "\n"
        "[phone] ;  ------------- phone message setting\n"        
        "activate:\tFalse\n"
        "\t\t;send message or not\n"
        "account:\t\n"        
        "\t\t;twilio account\n"
        "token:\t\n"        
        "\t\t;twilio token\n"
        "to:\t\n"        
        "\t\t;set message receiver\n"
        "from:\t\n"       
        "\t\t;set message sender\n"      
        "\n"
        "[database] ;-------------- define database configuration\n"
        "activate:\tFalse\n"
        "\t\t; activate database injection, configure in src/sqlconn.py\n"
        "[web] ;-------------- define web configuration\n"
        "activate:\tFalse\n"
        "\t\t; activate file transfer or not\n"
        "ssh:\t\n"
        "\t\t; if empty, local dir, else, remote dir, transfer files via ssh\n"
        "\t\t; user:password@server:port\n"             
        "dir:\t\n"        
        "\t\t; transfer files to web dir\n"
        "refresh:\tupdate.py\n"
        "\t\t; excute python file in web dir, if using html indtead of php,etc\n"
        "[wechat] ;-------------- define wechat configuration\n"
        "activate:\tFalse\n"
        "\t\t; activate or not\n"
        "friends:\t\n"
        "\t\t; send msg to whom\n"
        "groups:\t\n"
        "\t\t; send msg to which group\n"
        "[slack] ;-------------- define slack configuration\n"
        "activate:\tFalse\n"
        "\t\t; activate file transfer or not\n"
        "SLACK_BOT_TOKEN:\txoxb-418728624791-637128267495-B6ry0WaibOwhxawkkOW06pZG\n"
        "BOT_ID:\tUJR3S7VEK\n"
        "channel:\tDJBRPHUD8\n"
    )

    # optlist    
    if _class=='optlist':
        _log.write(
        "[telescope] ;-------------- define telescope parameters\n"
        "name:\t%s\n"%_tel+
        "\t\t;telescope name\n"
        "lat:\t%s\n"%lat+
        "\t\t;Latitude:\n"        
        "lon:\t%s\n"%lon+
        "\t\t;Longitude, format as above\n"
        "alt:\t%s\n"%alt+
        "\t\t;Altitude (in meter)\n"
        "fovw:\t%.2f\n"%fovw+
        "\t\t;field of view in widtw (in deg)\n"
        "fovh:\t%.2f\n"%fovh+
        "\t\t;field of view in height (in deg)\n"            
        "exptime:\t%.2f\n"%exptime+
        "\t\t; exposure time per frame (in second)\n"
        "rottime:\t%.2f\n"%rottime+
        "\t\t; estimated time for pointing reassignment (in second)\n"
        "resolution:\t%.2f\n"%_res+
        "\t\t; image resolution (arcsec/pixel)\n"
        "\n"
        "[pointings] ;-------------- define pointing parameters\n"
        "dbfile:\t\n"
        "\t\t; database mode: store pointings into dbfile.npz - ra,dec,id,score,info \n"
        "\t\t; next time, read  pointings and score \n"
        "\t\t; fill with %s.db to open, otherwise, leave blacked \n"%_tel+      
        "radecfile:\t\n"
        "\t\t; radec mode: next time, input only ra,dec, without score \n"
        "scheduler:\t%s\n"%strg+
        "\t\t; scheduler: [T]iling; [G]alaxy; [M]cmc\n"
        "limra:\t%s\n"%limra+
        "\t\t; exposure per frame\n"
        "limdec:\t%s\n"%limdec+
        "\t\t; average time for pointing reassignment\n"
        "ob:\t%s\n"%ob+
        "\t\t;generate Observing Blocks, for tiling\n"
        "nob:\t1\n"
        "\t\t;calculate visibility every nob*tf time\n"
        "\t\t;suggest: 1 for 3*3 OB mode; 10 for galaxy search\n"
        "limnob:\t100\n"
        "\t\t;maximum OB that can be done diring night\n"
        "\t\t;for telescope which has limitation on time, otherwise, leave blanket, or very large, e.g.9999\n"
        "\n"
        "[files] ;-------------- define file name\n"
        "fskip:\t\n"
        "\t\t; file of skipping pointing list\n"
        "\t\t; for bright stars or repeating fields before, etc\n"
        "vskip:\t1.\n"
        "\t\t; seperation in degree\n"        
        "schfile:\tpst_%s.txt\n"%_tel+
        "\t\t; output scheduler file format\n"
    )

        if _tel == 'VST':
            _log.write(
            "[scheduler]  ;  ------------- scheduler\n"
            "api:\tFalse\n"
            "\t\t; use api for scheduler\n"            
            "environment:\tdemo\n"
            "\t\t; VST environment\n"
            "username:\t52052\n"
            "\t\t; VST username\n"
            "password:\ttutorial\n"
            "\t\t; VST password\n"
            "containerid:\t1455705\n"
            "\t\t; VST containerid\n"
            "userpriority:\t1\n"
            "\t\t; OB priority\n"
            "airmass:\t3\n"
            "\t\t; airmass limit\n"
            "skyTransparency:\t2\n"
            "\t\t; sky transparency:\n"
            "\t\t; 1 - Photometric\n"
            "\t\t; 2 - Clear\n"
            "\t\t; 3 - Variable, thin cirrus\n"
            "\t\t; 4 - Variable, thick cirrus\n"
            "fli:\t1.0\n"
            "\t\t; fractional_lunar_illumination\n"
            "seeing:\t2.0\n"
            "\t\t; seeing limit\n"
            "dither:\t0\n"
            "\t\t; for VST, currently, use fixed dither\n"
            "repeat:\t1\n"
            "\t\t; repeat per field, if dither, at lease 1\n"
            "filter:\tr_SDSS\n"
            "\t\t; filter\n")

        elif _tel == 'schmidt':
            _log.write(
            "[scheduler]  ;  ------------- scheduler\n"          
            "dither:\t5\n"
            "\t\t; dither (e.g. 5, in pixel)\n"
            "repeat:\t1\n"
            "\t\t; repeat per field, if dither, at lease 1\n"
            "filter:\tr\n"
            "\t\t; filter\n")

        else:
            _log.write(
            "[scheduler]  ;  ------------- scheduler\n"           
            "dither:\t0\n"
            "\t\t; dither (e.g. 5, in arcsec)\n"
            "repeat:\t0\n"
            "\t\t; repeat per field, if dither, at lease 1\n"
            "filter:\tr\n"
            "\t\t; filter\n")

    _log.close()
    return configfile
