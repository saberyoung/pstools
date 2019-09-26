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

# trigger/healpix map:
def mollview(pparams):

    # read parameters
    hpmap = pparams['hpmap']
    theta = float(pparams['theta'])   
    phi = float(pparams['phi'])
    fignum = int(pparams['fignum'])
    title = pparams['title']
    ordering = eval(pparams['ordering'])
    coord = pparams['coord']
    norm = str(pparams['norm'])    
    if len(pparams['figsize'].split(','))==2:
        _figsize = (int(pparams['figsize'].split(',')[0]),\
                    int(pparams['figsize'].split(',')[1]))
    else:sys.exit('check figsize in par file')
    _min = float(pparams['min'])
    _max = float(pparams['max'])
    theta_contour = pparams['theta_contour']
    phi_contour = pparams['phi_contour']
    label = pparams['label']
    colors = [zz for zz in pparams['colors'].split(',')]
    distinfo = pparams['distinfo']
    tellist = pparams['tellist']
    timenow = pparams['timenow']

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
    if hp.__version__ >= '1.12.9':
        hp.mollview(map=pst.rotate_map(hpmap,_r), fig=fignum, rot=None, coord=coord, \
                unit='', xsize=800, title=title, nest=ordering, min=_min, max=_max, \
                flip='astro', remove_dip=False, remove_mono=False, gal_cut=0, \
                format='%g', format2='%g', cbar=False, cmap=cmap, badcolor='gray', \
                bgcolor='white', notext=False, norm='log', hold=True, margins=None, \
                sub=None, nlocs=2, return_projected_map=False)
    else:       
        hp.mollview(map=rot_map, fig=fignum, rot=None, coord=coord, \
                unit='', xsize=800, title=title, nest=ordering, min=_min, max=_max, \
                flip='astro', remove_dip=False, remove_mono=False, gal_cut=0, \
                format='%g', format2='%g', cbar=False, cmap=cmap, notext=False, \
                norm=norm, hold=True, margins=None, sub=None, nlocs=2, \
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
                _theta,_phi = _r(_theta_contour[i],_phi_contour[i])
                hp.projplot(_theta,_phi,linewidth=1,c=colors[ndd])

    # coordinate
    plot_coord(_r,coord)

    # write labels
    if eval(label):
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
    for _ntt in range(len(optlist)):

        # for each telescope
        _tt,_hlat,_hlon,_halt = optlist[_ntt]['telescope']['name'],\
                                optlist[_ntt]['telescope']['lat'],\
                                optlist[_ntt]['telescope']['lon'],\
                                optlist[_ntt]['telescope']['alt']       
        observatory = astropy.coordinates.EarthLocation(lat=float(_hlat)*u.deg, \
                                                        lon=float(_hlon)*u.deg, \
                                                        height=float(_halt)*u.m)
        _smlabel=True
        for _timeplus in np.arange(0,360,1):
            # define observatory
            newAltAzcoordiantes = astropy.coordinates.SkyCoord(alt = 0*u.deg, \
                                                               az = 0*u.deg + _timeplus*u.deg, \
                                                               obstime = timenow, \
                                                               frame = 'altaz', \
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

# galaxy map
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
    nbin=int(pparams['nbin'])
    label=pparams['label']
    fignum=pparams['fignum']
    scale=pparams['scale']
    # cut dist
    dist0 = dist0[np.logical_and(dist0>_distmin,
                                 dist0<_distmax)]

    fig = plt.figure(fignum)
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

    fig = plt.figure(fignum)
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

    fig = plt.figure(fignum, figsize=_figsize)
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
