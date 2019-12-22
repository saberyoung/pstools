"""############################################################################ 
2019/1/30 Start
define functions for plotting/visulization
""" ############################################################################
from __future__ import print_function
from builtins import input
import sys
import numpy as np
import healpy as hp
import math as mt
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import cm
import astropy.coordinates
import astropy.units as u
import pst

#################################################

""" 2d trigger/healpix map: """
def mollview(pparams):

    # read parameters
    hpmap = pparams['hpmap']
    theta = float(pparams['theta'])   
    phi = float(pparams['phi'])
    fignum = int(pparams['fignum'])
    title = pparams['title']
    ordering = pparams['ordering']
    coord = pparams['coord']
    norm = str(pparams['norm'])      
    _figsize = (int(pparams['figsize'].split(',')[0]),\
                int(pparams['figsize'].split(',')[1]))    
    _min = float(pparams['min'])
    _max = float(pparams['max'])
    theta_contour = pparams['theta_contour']
    phi_contour = pparams['phi_contour']
    colors = pparams['colors']
    distinfo = pparams['distinfo']
    tellist = pparams['tellist']
    timenow = pparams['timenow']
    if timenow == 'now': timenow = astropy.time.Time.now()            
    else:
        try: 
            timenow = astropy.time.Time(\
                    timenow, scale='utc')
        except: 
            timenow = astropy.time.Time.now() + \
                    astropy.time.TimeDelta(\
                    float(timenow)*60, format='sec')   

    # define rotate scheme
    _r = hp.Rotator(deg=True, rot=[phi,theta])
    # map is going reverse with coordinates
    _r1 = hp.Rotator(deg=True, rot=[-phi, np.pi/2. - theta])

    # rotate map
    nside = hp.npix2nside(len(hpmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))) #theta, phi

    # Get theta, phi under rotated co-ordinates
    trot, prot = _r1(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hpmap, trot, prot)

    # create figure
    fig = plt.figure(fignum, figsize=_figsize)

    # Choose color map and set background to white
    cmap = cm.YlOrRd
    cmap.set_under("w")

    # Plot GW skymap in Mollweide projection
    hp.mollview(map=rot_map, fig=fignum, rot=None, coord=coord, \
                unit='', xsize=800, title=title, nest=ordering, min=_min, max=_max, \
                flip='astro', remove_dip=False, remove_mono=False, gal_cut=0, \
                format='%g', format2='%g', cbar=False, cmap=cmap, notext=False, \
                norm=norm, hold=True, margins=None, sub=None, \
                return_projected_map=False)      

    hp.graticule() # Set grid lines

    # contour plots
    if theta_contour is not None and \
       phi_contour is not None:
        for ndd,dd in enumerate(theta_contour):
            # for .5, .9, ....
            _theta_contour, _phi_contour = theta_contour[dd], phi_contour[dd]        
            for i in range(len(_theta_contour)):
                # for each contour
                if len(_theta_contour[i])==0:continue
                _theta,_phi = _r(_theta_contour[i],_phi_contour[i])
                hp.projplot(_theta,_phi,linewidth=1,c=colors[ndd])

    # coordinate
    plot_coord(_r,coord)

    # write labels
    xx,yy = -2.,1.
    plt.text(xx, .9, 'sun: $\odot$',fontsize=20,\
             ha="center", va="center", color='y')
    plt.text(xx, 1., 'moon: $\oplus$',fontsize=20,\
             ha="center", va="center", color='b')        
    for ndd in range(len(colors)):                       
        yy+=.1
        try: plt.text(xx, yy, list(distinfo.values())[ndd],\
                      ha="center", va="center", color=colors[ndd],\
                      fontsize=20)
        except:pass

    # plot the sky: sun, moon, horizon, galactic plane, ...
    plot_sky(_r,tellist,coord,timenow,1)
    return fig

def plot_coord(r,coord='C'):
    """
    show specific coordiantes in healpy plots
    """

    _tlist,_plist = [60,120,180,240,300,360],\
                    [30,60,-30,-60]

    for _t in _tlist:
        # deg to hms
        c= astropy.coordinates.SkyCoord(ra=_t*u.degree, \
                                        dec=0*u.degree, frame='icrs')

        # select some special points
        theta,phi = pst.RadecToThetaphi(_t,0) 

        # apply rotation
        theta,phi = r(theta,phi)
        
        # visualization
        hp.projtext(theta,phi, '%ih'%c.ra.hms.h, coord=coord)

    for _p in _plist:
       
        # select some special points
        theta,phi = pst.RadecToThetaphi(0,_p) 

        # apply rotation
        theta,phi = r(theta,phi)
        
        # visualization
        hp.projtext(theta,phi, '%.f$^\circ$'%_p, coord=coord)

def plot_sky(r,optlist,coord,timenow,fignum=1):

    _colorlist = ['b','g','k','y','c','m']
    plt.figure(fignum)

    # plot the horizon
    for _st in optlist:
        for _ntt,_tt in enumerate(optlist[_st]):
            # for each telescope
            _tt,_hlat,_hlon,_halt = _tt['telescope']['name'],\
                                    _tt['telescope']['lat'],\
                                    _tt['telescope']['lon'],\
                                    _tt['telescope']['alt']       
            observatory = astropy.coordinates.EarthLocation(lat=float(_hlat)*u.deg, \
                                                            lon=float(_hlon)*u.deg, \
                                                            height=float(_halt)*u.m)
            _smlabel=True
            for _timeplus in np.arange(0,360,1):
                # define observatory
                # sometimes finals2000A is needed by astropy.utils.iers.iers, however blocked
                # download from: ftp://ftp.iers.org/products/eop/rapid/standard/finals2000A.all
                newAltAzcoordiantes = astropy.coordinates.SkyCoord(alt = 0*u.deg, \
                                        az = 0*u.deg + _timeplus*u.deg, \
                                        obstime = timenow, frame = 'altaz', \
                                        location = observatory)

                # transform to theta phi
                _htheta,_hphi = pst.RadecToThetaphi(newAltAzcoordiantes.icrs.ra.deg, \
                                        newAltAzcoordiantes.icrs.dec.deg)

                _htheta,_hphi = r(_htheta,_hphi)

                # plot
                if _smlabel:
                    hp.projplot(_htheta,_hphi,'.', color = _colorlist[_ntt], \
                            coord=coord, ms = 2, label='%s horizon now'%_tt)
                    _smlabel=False
                else:
                    hp.projplot(_htheta,_hphi,'.', color = _colorlist[_ntt], \
                            coord=coord, ms = 2)

    # plot the galactic plane        
    _hral = np.arange(0,360,10)
    _hdecl = np.zeros(len(_hral))
    _smlabel=True
    for _hra,_hdec in zip(_hral,_hdecl):
        
        # from galactic coordinates to equatorial
        _hradecs = astropy.coordinates.SkyCoord(l=_hra*u.deg, \
                                                b=_hdec*u.deg, frame='galactic')                

        # transform to theta phi
        _htheta,_hphi = pst.RadecToThetaphi(_hradecs.icrs.ra.deg,\
                                            _hradecs.icrs.dec.deg)
        _htheta,_hphi = r(_htheta,_hphi)

        # plot
        if _smlabel:
            hp.projplot(_htheta,_hphi,'x', color = 'k', \
                        coord=coord, ms = 10, label='galactic plane')
            _smlabel=False
        else:
            hp.projplot(_htheta,_hphi,'x', color = 'k', \
                         coord=coord, ms = 10)

    # plot the sun
    _sra,_sdec = astropy.coordinates.get_sun(timenow).\
                 ra.deg,astropy.coordinates.get_sun(timenow).dec.deg
    _stheta,_sphi = pst.RadecToThetaphi(_sra,_sdec)    
    _stheta,_sphi = r(_stheta,_sphi)
    hp.projtext(_stheta,_sphi,'$\odot$', color = 'y', \
                coord=coord,fontsize=20)     

    # plot the moon lasting 5 days
    for _nd in np.arange(0,5,1):
        _timeplus = _nd*24
        timelater = timenow + astropy.time.TimeDelta(_timeplus*3600, \
                                                     format='sec')                           
        _mra,_mdec = astropy.coordinates.get_moon(timelater).\
                     ra.deg,astropy.coordinates.get_moon(timelater).dec.deg
        _mtheta,_mphi = pst.RadecToThetaphi(_mra,_mdec)        
        _mtheta,_mphi = r(_mtheta,_mphi)          
        hp.projtext(_mtheta,_mphi,'$\oplus$', color = 'b', \
                    coord=coord,fontsize=20)
        hp.projtext(_mtheta,_mphi-.1,'%id'%_nd, color = 'b', \
                    coord=coord,fontsize=12)
    plt.legend()

""" galaxy plots """
def pointview(pparams):
   
    ra=pparams['ra']
    dec=pparams['dec']
    rot_phi=pparams['phi']
    rot_theta=pparams['theta']
    color=pparams['color']
    coord=pparams['coord']    
    label=pparams['label']
    fignum=pparams['fignum']
    ms=4
    if len(ra)>0:pass
    else:return

    fig = plt.figure(fignum)
    hp.graticule()
    _theta,_phi = pst.RadecToThetaphi(ra,dec)
    _rot = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])

    hp.projplot(_rot(_theta[0],_phi[0]),'x', color =color, coord=coord, ms = ms, label=label)
    hp.projplot(_rot(_theta,_phi),'x', color =color, coord=coord, ms = ms)

    plt.legend()
    return fig

def distview(pparams):

    _distmin=float(pparams['distmin'])
    _distmax=float(pparams['distmax'])
    dist0=pparams['dist']
    color1=pparams['color1']
    color2=pparams['color2']    
    try:nbin=int(pparams['nbin'])
    except:
        return ('!!! Warning: nbin should be an integer')            
    _figsize = (int(pparams['figsize'].split(',')[0]),\
                int(pparams['figsize'].split(',')[1]))    
    label=pparams['label']
    fignum=pparams['fignum']
    scale=pparams['scale']

    # cut dist
    dist0 = dist0[np.logical_and(dist0>_distmin,
                                 dist0<_distmax)]
    fig = plt.figure(fignum, _figsize)
    for nii,ii in enumerate(np.arange(_distmin,_distmax,abs(_distmax-_distmin)/nbin)):                 
        cum = len(dist0[np.logical_and(dist0<ii,dist0>_distmin)])
        pl.text(ii,-2,cum,color=color1)
    pl.hist(dist0,nbin,label=label,histtype='stepfilled',color=color2)
    if scale=='log': plt.yscale('log')
    plt.xlabel('Distance (Mpc)')
    plt.ylabel('Nuber of galaxies')
    plt.legend()
    return fig

def lumsview(pparams):

    _distmin=float(pparams['distmin'])
    _distmax=float(pparams['distmax'])    
    dist0=pparams['dist']   
    mag0=pparams['mag']
    color1=pparams['color1']  
    color2=pparams['color2']      
    nbin=int(pparams['nbin'])   
    label=pparams['label']
    fignum=pparams['fignum']
    scale=pparams['scale']    
    _figsize = (int(pparams['figsize'].split(',')[0]),\
                int(pparams['figsize'].split(',')[1]))
    fig = plt.figure(fignum, _figsize)
    Lums = 10**((4.74-mag0)/2.5)
    ticks,Lcum,Lbin = [],[],[]
    for ii in np.arange(_distmin,_distmax,nbin):
        ticks.append((ii+.5)*nbin)
        Lbin.append(sum(Lums[np.logical_and(dist0<ii,dist0>ii-nbin)]))
        Lcum.append(sum(Lums[np.logical_and(dist0<ii,dist0>_distmin)]))            
    pl.plot(ticks,Lbin,drawstyle="steps",\
            label='binned luminosity',color=color1)
    pl.plot(ticks,Lcum,drawstyle="steps",\
            label='cumulative luminosity',color=color2)
    pl.fill_between(ticks,Lbin,step="pre", alpha=1,color=color1)
    pl.title(label)   
    if scale=='log': plt.yscale('log')
    plt.xlabel('Distance (Mpc)')
    plt.ylabel('$L_{sun}$')
    plt.legend()     
    return fig

""" tilings plot """
def verticeview(pparams):

    ralist=pparams['ra']
    declist=pparams['dec']
    _fovw=pparams['fovw']
    _fovh=pparams['fovh']
    color=pparams['color']
    label=pparams['label']
    rot_phi=float(pparams['phi'])
    rot_theta=float(pparams['theta'])
    coord=pparams['coord']
    fignum=pparams['fignum']
        
    if len(ralist)>0:pass
    else:return

    fig = plt.figure(fignum)
    hp.graticule()

    try:
        float(_fovw)
        float(_fovh)
        plot=1
    except:
        plot=2

    if plot==1:
        ''' fov is a number '''
        print ('fov num')
        plot_lines(ralist,declist,float(_fovw),float(_fovh),\
                   rot_phi=rot_phi,rot_theta=rot_theta,\
                   coord=coord,color=color,label=label)
    elif plot==2:
        ''' fov is a list '''
        print ('fov list')
        for nm,(ra,dec,fovw,fovh) in enumerate(zip(ralist,declist,_fovw,_fovh)):
            if nm==0:_label=label
            else:_label=None
            plot_lines([ra],[dec],float(fovw),float(fovh),\
                rot_phi=rot_phi,rot_theta=rot_theta,coord=coord,\
                color=color,label=_label)
    try:
        _rank = pparams['rank']
        theta1,phi1 = pst.RadecToThetaphi(ralist,declist) 
        theta1,phi1 = r(theta1,phi1)  
        hp.projtext(theta1,phi1, str(_rank))
    except:
        pass
    plt.legend()
    return fig

def routeview(pparams):

    ralist=pparams['ra']
    declist=pparams['dec']
    _time=pparams['time']
    color=pparams['color']
    label=pparams['label']
    rot_phi=float(pparams['phi'])
    rot_theta=float(pparams['theta'])
    fignum=pparams['fignum']
    coord=pparams['coord']
    ms = 4

    if len(ralist)>0:pass
    else:return

    fig = plt.figure(fignum)
    hp.graticule()

    # start point
    _rot = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])
    _theta,_phi = pst.RadecToThetaphi(ralist[0],declist[0])
    hp.projplot(_rot(_theta,_phi),'*', color=color, \
                coord=coord, ms = ms)

    # routes
    theta,phi = pst.RadecToThetaphi(ralist,declist)
    hp.projplot(_rot(theta,phi),color=color,\
                coord=coord) #,label=label)

    plt.legend()
    return fig

def plot_lines(ra,dec,hh,ww,rot_theta=0,rot_phi=0,\
               color='k',coord='C',label=None):

    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])
    if not label is None: _smlabel=True
    else: _smlabel=False
    for _ra,_dec in zip(ra,dec):
        v1_ra,v2_ra,v3_ra,v4_ra,v1_dec,v2_dec,v3_dec,v4_dec=pst.vertices(_ra,_dec,hh,ww)
        ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra, v1_ra], \
                                     [v1_dec, v2_dec, v4_dec, v3_dec, v1_dec])
        theta,phi = pst.RadecToThetaphi(ra_vertices,dec_vertices)
        if _smlabel:
            hp.projplot(r(theta,phi),color=color,\
                        coord=coord,label=label)
            _smlabel=False
        else:hp.projplot(r(theta,phi),coord=coord,\
                         color=color)

def cumshow(pparams):

    # arg params
    full=pparams['full']
    fignum=pparams['fignum']
    select=pparams['select']    
    _figsize = (int(pparams['figsize'].split(',')[0]),\
                int(pparams['figsize'].split(',')[1]))

    # opt params
    try: nameloc=float(pparams['nameloc'])
    except: 
        print ('### Error: nameloc should be a float')
        return
    try: number=int(pparams['number'])
    except: 
        print ('### Error: number should be an interger')
        return
    try: showname=eval(pparams['showname'])
    except: 
        print ('### Error: showname should be True/False')
        return

    # initialize plot
    fig = plt.figure(fignum, figsize=_figsize)
    ax1=pl.axes([.1,.10,.4,.3])
    ax2=pl.axes([.1,.40,.4,.55])
    ax11=pl.axes([.5,.10,.35,.3])
    ax22=pl.axes([.5,.40,.35,.55])
    ax1.set_xlim([0,number+1])
    ax2.set_xlim([0,number+1])
    ax1.set_ylim([10**(-8),.3])
    ax1.set_yticks([10**(-8),.3])
    ax2.set_ylim([.01,1])
    ax1.set_yscale('log')
    ax1.set_xlabel(' '*50+'$N_{gal}$')      
    ax1.set_ylabel('Score')
    ax2.set_ylabel('Cumulative Score')
    ax2.tick_params(labelbottom=False) 
    ax11.set_ylim([10**(-8),.3])   
    ax22.set_ylim([.01,1])   
    ax11.set_yscale('log')
    ax11.set_xscale('log')    
    ax22.set_xscale('log')
    ax22.tick_params(labelbottom=False) 
    ax22.tick_params(labelleft=False) 
    ax11.tick_params(labelleft=False) 

    # plot
    _cc = 0
    _color = ['b','g','y','c','d']
    color1 = 'r'
    color2 = 'k'
    for _nt,_tt in enumerate(full):
       
        ##### full part
        scorei = np.array(pst.decomposit(_tt['score']))
        rai = np.array(pst.decomposit(_tt['ra']))
        deci = np.array(pst.decomposit(_tt['dec']))
        try: namei = np.array(pst.decomposit(_tt['name']))
        except: namei = np.array(pst.decomposit(_tt['ra']))

        # sort score
        idx = np.argsort(np.asarray(scorei))[::-1]

        scorei = scorei[idx]
        rai = rai[idx]
        deci = deci[idx]
        namei = namei[idx]

        # show only interesting fields
        _nm = 5
        score = scorei[:number*_nm]
        ra = rai[:number*_nm]
        dec = deci[:number*_nm]
        name = namei[:number*_nm]

        ax1.plot(np.arange(len(score))+1,score,color1+'.')
        ax2.plot(np.arange(len(score))+1,\
                [sum(score[:y]) for y in range(1, len(score) + 1)],\
                 color1+'.')   
        if showname:
            for i, txt in enumerate(name):
                ax2.annotate(txt, (i+1, \
                    sum(score[:i+1])+nameloc),\
                    fontsize=6,rotation=90)

        ax11.plot(np.arange(len(score))+1,score,color1+'.')
        ax22.plot(np.arange(len(score))+1,[sum(score[:y]) \
            for y in range(1, len(score) + 1)],color1+'.') 
        ax11.set_xlim([number+1,len(score)])    
        ax22.set_xlim([number+1,len(score)])
        ax11.set_xticks([number+1,len(score)])  

        ##### for selected dots
        tellist = None
        for _weight in select:
            for _ntt in select[_weight]:
                if select[_weight][_ntt]['name'] ==\
                   _tt['telescope']['name']:
                    tellist = select[_weight][_ntt]
        if not tellist is None:
            # for one telescope
            _title = '%i pointings with %.2f %% probability'%\
                     (number, sum(score[:number])*100)
            for _ra,_dec in zip(tellist['ra'],tellist['dec']):
                s,t = np.where(rai == _ra),\
                      np.where(deci == _dec)
                ids = np.intersect1d(s, t)                
                ax1.plot(ids+1,scorei[ids],_color[_cc]+'x')
                ax2.plot(ids+1,np.cumsum(scorei)[ids],_color[_cc]+'x')
            _title += '\n%s: %s'%\
                      (tellist['name'],_color[_cc])
            _cc+=1
    plt.title(_title)
    return fig

def interactive_show(func,_pm,_opt):   
    """
    Usage: intractive(plotfunc,rot_theta,'theta',rot_phi,'phi',[ralist,declist,_fovw,_fovh])
    """                
    answ = False
    while not answ:
        _fig = func(_pm)   
        _show = '    current options:'
        for ii in _opt:
            try: _show += ' %s=%.2e '%(ii,float(_pm[ii])) 
            except: _show += ' %s=%s '%(ii,_pm[ii])
        print(_show)
        answ = True
        good_answ=False
        while not good_answ: 
            good_answ = True
            try:
                _answ = input('>>> good options (y) or quit (q) or change'+\
                          ' options (eg. %s=%.2f) [y]: '%\
                          (_opt[0],float(_pm[_opt[0]])))    
            except:
                print ('!!! Error: parameter wrong')
                _answ = ''
                good_answ = False
            if len(_answ) > 0: 
                if _answ=='y': answ= True
                elif _answ=='q':sys.exit('quit...')
                else:
                    try:
                        if _answ[:_answ.index('=')] in _opt:
                            _pm[_answ[:_answ.index('=')]] = \
                                _answ[_answ.index('=')+1:].strip() 
                            try:_fig.clf()
                            except:print ('!!! Warning: no return from inter func')
                        else: 
                            print("!!! Warning; option "+_answ+\
                                     " not available !!!" )
                            good_answ = False
                    except: 
                        print("!!! Error: wrong input: retry")
                        good_answ = False
                    answ = False
    return _pm,_fig
