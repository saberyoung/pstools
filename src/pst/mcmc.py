def main(skymap,num,limra=False,limdec=False,fovh=3.,fovw=3.,fig=[1,2,3],radius=False,\
         rot_theta=0.,rot_phi=0.,verbose='2',intractive=False,ordering=False,\
         coord='C',norm = 'hist',color1='grey',color2='k',contour1=.99,contour2=.68):
    
    index1,index2,hpx = contour(skymap,.50,.99)  
    nside=hp.get_nside(hpx)
    dec,ra = IndexToDeclRa(nside,index2)
    if not limra:limra = [min(ra),max(ra)]
    if not limdec:limdec = [min(dec),max(dec)]

    # centering
#    print(np.mean(limra),np.mean(limdec))
#    rot_theta,rot_phi=0,0#np.mean(limra)#np.mean(limra),np.mean(limdec)

    if False:
        # contourview skymap show
        pparams = {'skymap':skymap,'contour1':contour1,'contour2':contour2,\
                   'rot_phi':rot_phi,'rot_theta':rot_theta,\
                   'color1':color1,'color2':color2,'label':'contour','coord':coord}
        optparams = ['rot_theta','rot_phi']
        if intractive:pstplot.intractive_show(pstplot.contourview,pparams,optparams)
        else:pstplot.contourview(pparams)  


    # monte carlo for tiling
    _log=[0.]
    _nloop=1
    for nn in [5,10,15,20]:
        print('searching in fovh/%i fovw/%i'%(nn,nn))
        shifth=fovh/nn
        shiftw=fovw/nn             
        answ = False
        showplot=True
        nloop=0
        while not answ:             
            answ = True
            good_answ=False            
            theta = random.uniform(0,2*np.pi)
            _print='\t with angle: %.2f'%theta
            _shifth,_shiftw = np.sqrt(shifth**2+shiftw**2)*mt.sin(theta),\
                              np.sqrt(shifth**2+shiftw**2)*mt.cos(theta)  
            while not good_answ: 
                good_answ = True                   
                _verbose = (verbose and showplot)                
                if _verbose:
                    # mollview skymap show
                    pparams = {'hpmap':skymap,'title':'skymap','rot_phi':rot_phi,\
                               'rot_theta':rot_theta,'fignum':fig[0],'ordering':ordering,\
                               'coord':coord,'norm':norm}
                    optparams = ['rot_theta','rot_phi']
                    pstplot.mollview(pparams)                     
                
                # generate pointings
                _ral,_decl=scheme.pointings(limdec=limdec,limra=limra,fovh=fovh,fovw=fovw,fig=fig[0],\
                                            shifth=_shifth,shiftw=_shiftw,rot_theta=rot_theta,\
                                            rot_phi=rot_phi,verbose=False,intractive=intractive)
                # cal prob for tiling list
                r2,d2,t = priorization.calprob_tile(skymap,_ral,_decl,fovh,fovw)    

                if _verbose:
                    # highlight selected
                    pparams = {'ra':r2[:num],'dec':d2[:num],'rot_phi':rot_phi,\
                               'rot_theta':rot_theta,'color':'r','fovw':fovw,'fovh':fovh}
                    optparams = ['rot_theta','rot_phi']
                    if intractive:pstplot.intractive_show(pstplot.verticeview,pparams,optparams)
                    else:pstplot.verticeview(pparams)

                if sum(t)>_log[-1]:                    
                    _log.append(sum(t))
                    good_answ=False
                    showplot=True                   

                    print(_print+'\tOK')
#                    plt.figure(fig[0])
                    plt.clf()
                    limra=[limra[0]+_shiftw,limra[1]+_shiftw]
                    limdec=[limdec[0]+_shifth,limdec[1]+_shifth]
                else:
                    nloop+=1
                    showplot=False
                    if nloop<_nloop:
                        answ=False                                                                    
                        print(_print+'\tno')
    
    _oo = open('radec.list','w')
    _ra,_dec,_score = [],[],[]
    print('Priorized pointings:')
    for ii,jj,kk in zip(r2[:num],d2[:num],t[:num]):
        print(ii,jj,kk)
        _ra.append(ii)
        _dec.append(jj)
        _score.append(kk)
        _oo.write('%.6f %.6f %.3f\n'%(ii,jj,kk))
    _oo.close()
    return _ra,_dec,_score
