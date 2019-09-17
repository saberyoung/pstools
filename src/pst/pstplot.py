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
import astropy.coordinates
import astropy.units as u
import pst

#################################################

# trigger/healpix map:
def mollview(pparams):

    hpmap=pparams['hpmap']
    rot_phi=pparams['rot_phi']
    rot_theta=pparams['rot_theta']
    fignum=pparams['fignum']
    title=pparams['title']
    ordering=pparams['ordering']
    coord=pparams['coord']
    norm=pparams['norm'] 

    fig = plt.figure(fignum, figsize=(15, 10))

#    hp.mollview(map=None, fig=None, rot=None, coord=None, unit='', xsize=800, title='Mollweide view', nest=False, min=None, max=None, flip='astro', remove_dip=False, remove_mono=False, gal_cut=0, format='%g', format2='%g', cbar=True, cmap=None, notext=False, norm=None, hold=False, margins=None, sub=None, nlocs=2, return_projected_map=False)

#rotate_map(hpmap,rot_theta,rot_phi),fig=fignum,\
#                title=title,unit='prob',nest=ordering,coord=coord,norm=norm)
    hp.mollview(hpmap,rot=(rot_theta,rot_phi),fig=fignum,\
                title=title,unit='prob',nest=ordering,coord=coord,norm=norm)

    # coordinate
    plot_coord(rot_phi,rot_theta)  
    return fig

def mollzoom(pparams):

    hpmap=pparams['hpmap']
    rot_phi=pparams['rot_phi']
    rot_theta=pparams['rot_theta']
    fignum=pparams['fignum']
    title=pparams['title']
    ordering=pparams['ordering']
    coord=pparams['coord']
    norm=pparams['norm']   
    
    # zoom
    nside = hp.get_nside(hpmap)
    dec,ra = IndexToDeclRa(nside,hpmap.argsort()[-1])

    hp.mollzoom(rotate_map(hpmap,rot_theta,rot_phi),rot=(ra,dec),xsize=2000,fig=fignum,\
                title=title,unit='prob',nest=ordering,coord=coord,norm=norm)
    # coordinate
    plot_coord(rot_phi,rot_theta) 

def contourview(pparams):

    skymap=pparams['hpmap']   
    rot_phi=pparams['rot_phi']
    rot_theta=pparams['rot_theta']
    color1='grey'
    color2='y'
    color3='r'
    color4='w'
    _label='GW loc'
    _coord=pparams['coord']
    fignum=pparams['fignum']

    fig = plt.figure(fignum, figsize=(15, 10))    

    hp.graticule()
    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])
    _ilist,hpx = pst.contour(skymap)  
    label1,label2,label3,label4=_ilist.keys()
    index1,index2,index3,index4=_ilist.values()
    nside = hp.get_nside(hpx)
    for _index,_color,_cont in zip([index4,index3,index2,index1],\
                                   [color1,color2,color3,color4],\
                                   [label4,label3,label2,label1]):
        if len(_index)==0:continue
        theta,phi = hp.pix2ang(nside,_index)
        thetar,phir = r(theta,phi)
        theta,phi = np.array(thetar),np.array(phir)
        try:
            hp.projplot(theta[0],phi[0],_color,coord=_coord,\
                        label='%s %s'%(_label,_cont))
            hp.projplot(theta,phi,_color,coord=_coord)
        except:
            hp.projplot(theta,phi,_color,coord=_coord,\
                        label='%s %s'%(_label,_cont))

    plt.title('trigger localization in equatorial system')
    plot_coord(rot_phi,rot_theta)
    plt.axis('off')
    plt.legend()
    return fig

# galaxy map
def pointview(pparams):
   
    ra=pparams['ra']
    dec=pparams['dec']
    rot_phi=pparams['rot_phi']
    rot_theta=pparams['rot_theta']
    color=pparams['color']
    coord=pparams['coord']
    ms=pparams['ms']
    label=pparams['label']
    fignum=pparams['fignum']

    if len(ra)>0:pass
    else:return

    fig = plt.figure(fignum)
    hp.graticule()
    _theta,_phi = pst.RadecToThetaphi(ra,dec)
    _rot = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])

    hp.projplot(_rot(_theta[0],_phi[0]),'x', color =color, coord=coord, ms = ms, label=label)
    hp.projplot(_rot(_theta,_phi),'x', color =color, coord=coord, ms = ms)
#    plot_coord(rot_phi,rot_theta)
    plt.legend()
    return fig

def distview(pparams):

    _distmin=float(pparams['distmin'])
    _distmax=float(pparams['distmax'])
    dist00=pparams['distfull_mat']
    dist0=pparams['dist_mat']
    color1=pparams['color1']
    color2=pparams['color2']
    color3=pparams['color3']
    _nbindist=int(pparams['nbin'])
    label1=pparams['label1']
    label2=pparams['label2']
    ypos1=float(pparams['ypos1'])
    ypos2=float(pparams['ypos2'])
    fignum=pparams['fignum']

    fig = plt.figure(fignum)
    for nii,ii in enumerate(np.arange(_distmin,_distmax,abs(_distmax-_distmin)/10.)): 
        cum = len(dist00[np.logical_and(dist00<ii,dist00>_distmin)])
        pl.text(ii,ypos1,cum,color=color1)           
        cum = len(dist0[np.logical_and(dist0<ii,dist0>_distmin)])
        pl.text(ii,ypos2,cum,color=color2)          
    pl.hist(dist00[np.logical_and(dist00<_distmax,dist00>_distmin)], \
            _nbindist,label=label1,histtype='step',color=color3)       
    pl.hist(dist0,_nbindist,label=label2,histtype='stepfilled',color=color3)         
    plt.xlabel('Distance (Mpc)')
    plt.ylabel('Nuber of galaxies')
#    pl.tight_layout()
    plt.legend()
    return fig

def lumsview(pparams):

    _distmin=float(pparams['distmin'])
    _distmax=float(pparams['distmax'])
    dist00=pparams['distfull_mat']
    dist0=pparams['dist_mat']
    mag00=pparams['magfull_mat']
    mag0=pparams['mag_mat']
    color1=pparams['color1']  
    color2=pparams['color2']  
    color3=pparams['color3']  
    _nbinlums=int(pparams['nbin'])   
    label=pparams['label']
    fignum=pparams['fignum']

    fig = plt.figure(fignum)

    L1 = 10**((4.74-mag0)/2.5)
    L2 = 10**((4.74-mag00)/2.5)
    Lb1,Lb2 = [],[]
    for ii in np.arange(_distmin,_distmax,_nbinlums):                   
        Lb1.append(sum(L1[np.logical_and(dist0<ii,dist0>_distmin)]))
        Lb2.append(sum(L2[np.logical_and(dist00<ii,dist00>_distmin)]))
    ticks,Lbin1,Lbin2,Lbin3 = [],[],[],[]
    for ii in range(len(Lb1)):
        ticks.append((ii+.5)*_nbinlums)
        if ii==0:
            Lbin1.append(Lb1[ii])
            Lbin2.append(Lb2[ii])
        else:
            Lbin1.append(Lb1[ii]-Lb1[ii-1])
            Lbin2.append(Lb2[ii]-Lb2[ii-1])  
        
    pl.fill_between(ticks,Lbin1, step="pre", alpha=1,color='grey')
    pl.plot(ticks,Lbin2, drawstyle="steps",label='binned luminosity for %s'%label,color=color1)
    pl.plot(ticks,Lbin1, drawstyle="steps",label='binned luminosity for GLADE',color=color1)
    pl.plot(ticks,Lb2,color2,label='cumulative luminosity for %s'%label)
    pl.plot(ticks,Lb1,color3,label='cumulative luminosity for GLADE')   
    plt.yscale('log')
    plt.xlabel('Distance (Mpc)')
    plt.ylabel('$L_{sun}$')
#    pl.tight_layout()
    plt.legend()     
    return fig

def dist_gauss(mu,sigma):    

    # no return

    from matplotlib import mlab
    ax = plt.subplot(121)
    ################
    x = np.linspace(0, mu + sigma, 100)
    ax.plot(x,mlab.normpdf(x, mu, sigma),label=str(mu)+'+/-'+str(sigma))
    ax.set_yscale('log')

    ax.grid(True)
    ax.legend()
    ax.set_title('Gaussian distribution')
    ax.set_xlabel('distance (Mpc)')
    ax.set_ylabel('probability')

    ##################
    ax = plt.subplot(122)
    x = np.random.normal(mu, sigma, size=100)
    # plot the cumulative histogram
    n, bins, patches = ax.hist(x, 50, normed=1, histtype='step',
                               cumulative=True, label='Empirical')

    # Add a line showing the expected distribution.
    y = mlab.normpdf(bins, mu, sigma).cumsum()
    y /= y[-1]
    
    ax.plot(bins, y, 'k--', linewidth=1.5, label='Theoretical')

    ax.grid(True)
    ax.legend()
    ax.set_title('Cumulative step histograms')
    ax.set_xlabel('distance (Mpc)')
    ax.set_ylabel('cumulative probability')

def verticeview(pparams):

    ralist=pparams['ra']
    declist=pparams['dec']
    _fovw=pparams['fovw']
    _fovh=pparams['fovh']
    color=pparams['color']
    label=pparams['label']
    rot_phi=pparams['rot_phi']
    rot_theta=pparams['rot_theta']
    fignum=pparams['fignum']
        
    if len(ralist)>0:pass
    else:return

    fig = plt.figure(fignum)
    hp.graticule()
    plot_lines(ralist,declist,_fovw,_fovh,rot_phi=rot_phi,\
               rot_theta=rot_theta,color=color,label=label)
    try:
        _rank = pparams['rank']
        theta1,phi1 = pst.RadecToThetaphi(ralist,declist) 
        theta1,phi1 = r(theta1,phi1)  
        hp.projtext(theta1,phi1, str(_rank))
    except:pass
    plt.legend()
#    plot_coord(rot_phi,rot_theta) 
    return fig

def plot_coord(rot_phi,rot_theta):
    """
    show specific coordiantes in healpy plots
    """

    # define rotate scheme
    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])

    _ralist,_declist = [0,45,90,135,180,225,270,315],\
                       [0,30,60,-30,-60]

    for _ra in _ralist:
        for _dec in _declist:

            # select some special points
            theta1,phi1 = pst.RadecToThetaphi(_ra,_dec) 

            # apply rotation
            theta1,phi1 = r(theta1,phi1)
  
            # visualization
            if _ra == 0:
                hp.projtext(theta1,phi1, str(_dec), coord='C')
            elif _dec == 0:
                hp.projtext(theta1,phi1, str(_ra), coord='C')
            else:
                hp.projtext(theta1,phi1, str(_ra)+','+str(_dec), coord='C')

def plot_lines(ra,dec,hh,ww,rot_theta=0,rot_phi=0,color='k',label='tiles'):

    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])
    _smlabel=True
    for _ra,_dec in zip(ra,dec):
        v1_ra,v2_ra,v3_ra,v4_ra,v1_dec,v2_dec,v3_dec,v4_dec=pst.vertices(_ra,_dec,hh,ww)
        ra_vertices, dec_vertices = ([v1_ra, v2_ra, v4_ra, v3_ra, v1_ra], \
                                     [v1_dec, v2_dec, v4_dec, v3_dec, v1_dec])
        theta,phi = pst.RadecToThetaphi(ra_vertices,dec_vertices)
        if _smlabel:
            hp.projplot(r(theta,phi),color,label=label)
            _smlabel=False
        else:hp.projplot(r(theta,phi),color)

def plot_sky(optlist,timenow,_colorlist):

    if True:
        # plot the horizon
        for _ntt,_tt in enumerate(optlist):
            if _tt == 'arg':continue

            # for each telescope
            _hlat,_hlon,_halt = optlist[_tt]['telescope']['lat'],\
                                optlist[_tt]['telescope']['lon'],\
                                optlist[_tt]['telescope']['alt']

            observatory = astropy.coordinates.EarthLocation(lat=float(_hlat)*u.deg, \
                                                            lon=float(_hlon)*u.deg, \
                                                            height=float(_halt)*u.m)
            _smlabel=True
            for _timeplus in np.arange(0,360,1):
                newAltAzcoordiantes = astropy.coordinates.SkyCoord(alt = 0*u.deg, \
                                                az = 0*u.deg + _timeplus*u.deg, \
                                                obstime = timenow, \
                                                frame = 'altaz', \
                                                location = observatory)
                # transform to theta phi
                _htheta,_hphi = pst.RadecToThetaphi(newAltAzcoordiantes.icrs.ra.deg, \
                                                    newAltAzcoordiantes.icrs.dec.deg)

                plt.figure(0)
                if _smlabel:
                    hp.projplot(_htheta,_hphi,'.', color = _colorlist[_ntt], \
                                coord=optlist['arg']['plot']["coord"], ms = 2, label='%s horizon now'%_tt)
                    _smlabel=False
                else:
                    hp.projplot(_htheta,_hphi,'.', color = _colorlist[_ntt], \
                                coord=optlist['arg']['plot']["coord"], ms = 2)

        # plot the galactic plane        
        _hral = np.arange(0,360,10)
        _hdecl = np.zeros(len(_hral))
        _smlabel=True
        for _hra,_hdec in zip(_hral,_hdecl):

            # from galactic coordinates to equatorial
            _hradecs = astropy.coordinates.SkyCoord(l=_hra*u.deg, b=_hdec*u.deg, frame='galactic')                

            # transform to theta phi
            _htheta,_hphi = pst.RadecToThetaphi(_hradecs.icrs.ra.deg, _hradecs.icrs.dec.deg)

            plt.figure(0)
            if _smlabel:
                hp.projplot(_htheta,_hphi,'x', color = 'k', coord=optlist['arg']['plot']["coord"], ms = 10, label='galactic plane')
                _smlabel=False
            else:hp.projplot(_htheta,_hphi,'x', color = 'k', coord=optlist['arg']['plot']["coord"], ms = 10)

        # plot the sun, moon
        _smlabel=True
        for _timeplus in np.arange(0,36,6):
            timelater = timenow + astropy.time.TimeDelta(_timeplus*3600, format='sec') 

            _sra,_sdec = astropy.coordinates.get_sun(timelater).ra.deg,astropy.coordinates.get_sun(timelater).dec.deg
            _stheta,_sphi = pst.RadecToThetaphi(_sra,_sdec)              

            _mra,_mdec = astropy.coordinates.get_moon(timelater).ra.deg,astropy.coordinates.get_moon(timelater).dec.deg
            _mtheta,_mphi = pst.RadecToThetaphi(_mra,_mdec)

            plt.figure(0)
            hp.projplot(_stheta,_sphi,'o', color = 'y', coord=optlist['arg']['plot']["coord"], ms = 10)            
            hp.projplot(_mtheta,_mphi,'o', color = 'b', coord=optlist['arg']['plot']["coord"], ms = 10)
            if _smlabel:
                hp.projtext(_stheta,_sphi,'sun', color = 'k', coord=optlist['arg']['plot']["coord"])
                hp.projtext(_mtheta,_mphi,'moon', color = 'k', coord=optlist['arg']['plot']["coord"])
                _smlabel=False
    plt.legend()

def cumshow(pparams):

    score=pparams['score']
    _number=int(pparams['number'])
    name=pparams['name']
    fignum=pparams['fignum']
    if 'color1' in pparams: color1=pparams['color1']
    else:color1 = 'r'
    if 'color2' in pparams: color2=pparams['color2']
    else:color2 = 'k'

    fig = plt.figure(fignum)
    ax1=pl.axes([.1,.10,.4,.3])
    ax2=pl.axes([.1,.40,.4,.55])
    ax11=pl.axes([.5,.10,.35,.3])
    ax22=pl.axes([.5,.40,.35,.55])
    plt.title('%i pointings with %.2f %% probability' % (_number, sum(score[:_number])*100))

    ax1.plot(np.arange(len(score))+1,score,color1+'.')
    ax2.plot(np.arange(len(score))+1,[sum(score[:y]) for y in range(1, len(score) + 1)],color1+'.')    
    ax1.set_ylim([10**(-8),.3])
    ax1.set_yticks([10**(-8),.3])
    ax2.set_ylim([.1,1])
    ax1.set_xlim([0,_number+1])
    ax2.set_xlim([0,_number+1])
    ax1.set_yscale('log')
    for i, txt in enumerate(name):
        if i == 0:            
            ax2.annotate(txt, (i+1, sum(score[:i+1])+.15),fontsize=6,rotation=90)
        else:
            ax2.annotate(txt, (i+1, sum(score[:i+1])-.05),fontsize=6,rotation=90)

    ax1.set_xlabel(' '*50+'$N_{gal}$')      
    ax1.set_ylabel('Score')
    ax2.set_ylabel('Cumulative Score')
    ax2.tick_params(labelbottom=False) 

    ax11.plot(np.arange(len(score))+1,score,'k.')
    ax22.plot(np.arange(len(score))+1,[sum(score[:y]) for y in range(1, len(score) + 1)],'r.')   
    ax11.set_ylim([10**(-8),.3])   
    ax22.set_ylim([.1,1])
    ax11.set_xlim([_number+1,len(score)])    
    ax22.set_xlim([_number+1,len(score)])
    ax11.set_xticks([_number+1,len(score)])    
    ax11.set_yscale('log')
    ax11.set_xscale('log')    
    ax22.set_xscale('log')
    ax22.tick_params(labelbottom=False) 
    ax22.tick_params(labelleft=False) 
    ax11.tick_params(labelleft=False) 
    return fig

def cumshow1(pparams):

    # highlight some dots

    # before
    _score=pparams['score']
    _idx = range(len(_score))   

    # sort
    idx=np.argsort(np.asarray(_score))[::-1]
    score=_score[idx]      

    # cut
    if len(score)>2000: score = score[:2000]

    _number=pparams['number']
    fignum=pparams['fignum']
    _highlight=pparams['highlight']

    _color=pparams['color']
    if 'color1' in pparams: color1=pparams['color1']
    else:color1 = 'r'
    if 'color2' in pparams: color2=pparams['color2']
    else:color2 = 'k'

    fig = plt.figure(fignum, figsize=(15, 15))
    ax1=pl.axes([.1,.10,.4,.3])
    ax2=pl.axes([.1,.40,.4,.55])
    ax11=pl.axes([.5,.10,.35,.3])
    ax22=pl.axes([.5,.40,.35,.55])
    plt.title('%i pointings with %.2f %% probability' % (_number, sum(score[:_number])*100))
    
    ax1.plot(np.arange(len(score))+1,score,color1+'.')
    ax2.plot(np.arange(len(score))+1,[sum(score[:y]) for y in range(1, len(score) + 1)],color1+'.')   

    ax1.set_ylim([10**(-8),.3])
    ax1.set_yticks([10**(-8),.3])
    ax2.set_ylim([.1,1])
    ax1.set_xlim([0,_number+1])
    ax2.set_xlim([0,_number+1])
    ax1.set_yscale('log')

    ''' remove the name
    for i, txt in enumerate(idx):
        if i == 0:ax2.annotate(txt, (i+1, sum(score[:i+1])+.15),fontsize=6,rotation=90)
        else:ax2.annotate(txt, (i+1, sum(score[:i+1])-.05),fontsize=6,rotation=90)
    '''

    ax1.set_xlabel(' '*50+'$N_{gal}$')      
    ax1.set_ylabel('Score')
    ax2.set_ylabel('Cumulative Score')
    ax2.tick_params(labelbottom=False) 

    ax11.plot(np.arange(len(score))+1,score,'k.')
    ax22.plot(np.arange(len(score))+1,[sum(score[:y]) for y in range(1, len(score) + 1)],'r.')   
    ax11.set_ylim([10**(-8),.3])   
    ax22.set_ylim([.1,1])
    ax11.set_xlim([_number+1,len(score)])    
    ax22.set_xlim([_number+1,len(score)])
    ax11.set_xticks([_number+1,len(score)])    
    ax11.set_yscale('log')
    ax11.set_xscale('log')    
    ax22.set_xscale('log')
    ax22.tick_params(labelbottom=False) 
    ax22.tick_params(labelleft=False) 
    ax11.tick_params(labelleft=False) 

    # for highlighted dots
    _cc=0
    for _hh in _highlight:        

        # every time not visible
        if len(_highlight[_hh])==0:
            print('Warning: No fields visible for %s!'%_hh)
            continue

        # for one telescope
        for _nh,_hlist in enumerate(_highlight[_hh]):

            # some time not visible
            if _hlist == False:continue

            _id_special = np.argwhere(_idx ==_hlist)
            id_special = np.argwhere(idx ==_hlist)
            try:
                ax1.plot(id_special+1,score[id_special],_color[_cc]+'x')
                ax2.plot(id_special+1,np.cumsum(score)[id_special],_color[_cc]+'x')
            except:
                ax1.plot(id_special+1,_score[_id_special],_color[_cc]+'|')
        _cc+=1
    return fig

def vis_plot(ralist,declist,lat,lon,height,timelist):
    '''TBD'''

    from astropy.visualization import astropy_mpl_style
    plt.style.use(astropy_mpl_style)
    import astropy.units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz

    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(ra=ralist*u.deg, dec=declist*u.deg)
    
    bear_mountain = astropy.coordinates.EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    utcoffset = -4*u.hour  # Eastern Daylight Time
    time = Time('2012-7-12 23:00:00') - utcoffset

    m33altaz = m33.transform_to(AltAz(obstime=time,location=bear_mountain))
    print("M33's Altitude = {0.alt:.2}".format(m33altaz))

    midnight = Time('2012-7-13 00:00:00') - utcoffset
    delta_midnight = np.linspace(-2, 10, 100)*u.hour
    frame_July13night = AltAz(obstime=midnight+delta_midnight,
                              location=bear_mountain)
    m33altazs_July13night = m33.transform_to(frame_July13night)

    m33airmasss_July13night = m33altazs_July13night.secz

    plt.plot(delta_midnight, m33airmasss_July13night)
    plt.xlim(-2, 10)
    plt.ylim(1, 4)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Airmass [Sec(z)]')
    plt.show()

    from astropy.coordinates import get_sun
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bear_mountain)
    sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)

    from astropy.coordinates import get_moon
    moon_July12_to_13 = get_moon(times_July12_to_13)
    moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)
    
    m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)

    plt.plot(delta_midnight, sunaltazs_July12_to_13.alt, color='r', label='Sun')
    plt.plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75]*3, ls='--', label='Moon')
    plt.scatter(delta_midnight, m33altazs_July12_to_13.alt,
                c=m33altazs_July12_to_13.az, label='M33', lw=0, s=8,
                cmap='viridis')
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_July12_to_13.alt < -0*u.deg, color='0.5', zorder=0)
    plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_July12_to_13.alt < -18*u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.xlim(-12, 12)
    plt.xticks(np.arange(13)*2 -12)
    plt.ylim(0, 90)
    plt.xlabel('Hours from EDT Midnight')
    plt.ylabel('Altitude [deg]')
    plt.show()

def plot_all(params=''):#skymap,ralist,declist,problist,rot_theta=0.,rot_phi=0.,num=100,scheme='mollview',\
#             _label='trigger',contour1=.68,contour2=.99,color1='r',color2='b',\
#             color3='k.',color4='y',_ordering=False,_coord='C',_norm = 'hist',\
#             _type='polygon',_height=2.,_width=2.):    
    input(params)
    ##  trigger
    if scheme == 'contour':
        contour_show(skymap,contour1,contour2,rot_phi,\
                     rot_theta,color1,color2,_coord,_label)
    elif scheme == 'mollview':
        hpmap = hp.read_map(skymap, h=False, verbose=False)
        hp.mollview(rotate_map(hpmap,rot_theta,rot_phi),fig=None,title=_label, \
                    unit='prob',nest=_ordering,coord=_coord,norm=_norm) 
 
    else:print('Warning: wrong scheme for trigger map')
    
    theta = 0.5 * np.pi - np.deg2rad(declist)
    phi = np.deg2rad(ralist)
    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])
    thetar,phir = r(theta,phi)
    theta,phi = np.array(thetar),np.array(phir)

    ##  pointings
    if _type == 'dots':     
        hp.graticule()         
        hp.projplot(theta[0],phi[0],color3,coord=_coord,ms=1,label='pointings')
        hp.projplot(theta[:num],phi[:num],color3,ms=1,coord=_coord)
#        plot_coord(rot_phi,rot_theta)   
        plt.legend(loc=(.15,-.2),ncol=4,borderaxespad=0.,fontsize = 8)
    elif _type == 'polygon':
        plot_lines(ralist[:num],declist[:num],_width,_height,color=color4,\
                   rot_phi=rot_phi,rot_theta=rot_theta)
#        plot_coord(rot_theta,rot_phi)
        plt.legend(loc=(.15,-.2),ncol=4,borderaxespad=0.,fontsize = 8)               

def interactive_show(func,_pm,_opt):   
    """
    Usage: intractive(plotfunc,rot_theta,'theta',rot_phi,'phi',[ralist,declist,_fovw,_fovh])
    """                
    answ = False
    while not answ:                   
        func(_pm)   
        _show = '    current options:'
        for ii in _opt:_show += ' %s=%.2f '%(ii,float(_pm[ii])) 
        print(_show)
        answ = True
        good_answ=False
        while not good_answ: 
            good_answ = True
            _answ = input('>>> good options (y) or quit (q) or change'+\
                          ' options (eg. %s=%.2f) [y]: '%\
                          (_opt[0],float(_pm[_opt[0]])))
            if len(_answ) > 0: 
                if _answ=='y': answ= True
                elif _answ=='q':sys.exit('quit...')
                else:
                    try:
                        if _answ[:_answ.index('=')] in _opt:
                            _pm[_answ[:_answ.index('=')]] = \
                                _answ[_answ.index('=')+1:].strip() 
                            pl.clf()                        
                        else: 
                            print("!!! Warning; option "+_answ+\
                                     " not available !!!" )
                            good_answ = False
                    except: 
                        print("!!! Error: wrong input: retry")
                        good_answ = False
                    answ = False

def show_scheduler():
    # TBD
    _showchoice = 2

    # showfile
    if _showfile:
        pl.ion()
        # galaxies
        if os.path.exists(_cachefile):     
            print('show galaxies from %s'%_cachefile)
            cat0 = Table(np.load(_cachefile)['cat'])                  
            ra0,dec0 = cat0['RAJ2000'],\
                       cat0['DEJ2000']

            galpixels_Range_num= np.zeros(hp.nside2npix(nside))
            pix_num_Range    = (DeclRaToIndex(dec0,ra0,nside))
            galpixels_Range_num[pix_num_Range]+=1

            pparams = {'hpmap':galpixels_Range_num,'title':'galaxies','rot_phi':rot_phi,\
                       'rot_theta':rot_theta,'fignum':1,'ordering':_ordering,\
                       'coord':_coord,'norm':_norm}
            optparams = ['rot_theta','rot_phi']
            pstplot.mollview(pparams)
        else:print('Notice: if wanna galaxy distribution, do cache mode to store a galaxy catalog:\ttmp_glade.npz')

        # pointings       
        _flist = read_filelist(_showfile)
        colorlist = ['r','g','b','k','y']
        __pl = {}
        __pl['ra'],__pl['dec'],__pl['count']=[],[],[]

        for __jj,_f in enumerate(_flist):
            print('read %s'%_f)
            if 'npz' in _f:
                # DB file                
                try:
                    __idl = np.load(_f)['id']
                    __ralob = np.load(_f)['ra']
                    __declob = np.load(_f)['dec']
                    __score = np.load(_f)['score']                
                except:sys.exit('Error:check input showfile!')
                for __ra0,__dec0 in zip(__ralob,__declob):   
                    # OB level
                    for __ra1,__dec1 in zip(__ra0,__dec0): 
                        # pointing level                      
                        __same = np.where(np.sqrt((__pl['ra']-__ra1)**2+(__pl['dec']-__dec1)**2)<_dither) 
                        if np.array(__pl['count'])[__same]:
                            # check if exists already                            
                            __pl['count'][__same[0][0]]+=1
                        else:
                            __pl['ra'].append(__ra1)
                            __pl['dec'].append(__dec1)
                            __pl['count'].append(1)    
                if _showchoice == 2:
                    pparams = {'ra':__pl['ra'],'dec':__pl['dec'],'rot_phi':rot_phi,'fignum':1,\
                               'rot_theta':rot_theta,'color':colorlist[__jj%5],\
                               'fovw':_fovw*_obfovw,'fovh':_fovh*_obfovh}
                    optparams = ['rot_theta','rot_phi']
                    pstplot.verticeview(pparams)
                    __pl = {}
                    __pl['ra'],__pl['dec'],__pl['count']=[],[],[]

            else:
                # asc file                
                for _ii in open(_f).readlines():
                    if _ii[0]=='#' or len(_ii)==0:continue                    
                    __ra1 = np.float64(_ii.split()[0])
                    __dec1 = np.float64(_ii.split()[1])
                    __same = np.where(np.sqrt((__pl['ra']-__ra1)**2+(__pl['dec']-__dec1)**2)<_dither)                    
                    if np.array(__pl['count'])[__same]:
                        # check if exists already                            
                        __pl['count'][__same[0][0]]+=1
                    else:
                        __pl['ra'].append(__ra1)
                        __pl['dec'].append(__dec1)
                        __pl['count'].append(1)   
                if _showchoice == 2:
                    pparams = {'ra':__pl['ra'],'dec':__pl['dec'],'rot_phi':rot_phi,'fignum':1,\
                               'rot_theta':rot_theta,'color':colorlist[__jj%5],\
                               'fovw':_fovw*_obfovw,'fovh':_fovh*_obfovh}
                    optparams = ['rot_theta','rot_phi']
                    pstplot.verticeview(pparams)                            
                    __pl = {}
                    __pl['ra'],__pl['dec'],__pl['count']=[],[],[]

        if _showchoice == 1:
            # plot            
            for __ii,cond in enumerate([np.where(np.array(__pl['count'])==1),
                                    np.where(np.array(__pl['count'])==2),
                                    np.where(np.array(__pl['count'])==3),
                                    np.where(np.array(__pl['count'])==4),
                                    np.where(np.array(__pl['count'])==5)]):
                __ral = np.array(__pl['ra'])[cond]
                __decl = np.array(__pl['dec'])[cond]                  
          
                pparams = {'ra':__ral,'dec':__decl,'rot_phi':rot_phi,'fignum':1,\
                           'rot_theta':rot_theta,'color':colorlist[__ii%5],'fovw':_fovw,'fovh':_fovh}
                optparams = ['rot_theta','rot_phi']
                pstplot.verticeview(pparams)
        input('check')       
        sys.exit('Quit...')
