def query_ebv(ra,dec,size=2,thresh=.25,verbose=False):
    """
    Take ra,dec and return the E(B-V) number
    URL query, see https://irsa.ipac.caltech.edu/applications/DUST
    """
    import urllib2
    import xmltodict
    url = "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?"
    url += "locstr=%.2f+%.2f+equ+j2000"%(ra,dec)
    if size<2:size=2
    if size>37.5:size=37.5
    url += '&regSize=%.2f'%size # size between 2.0 and 37.5    

    _file = urllib2.urlopen(url)
    data = _file.read()
    _file.close()

    _dict = xmltodict.parse(data)
    _ebv = _dict['results']['result'][0]['statistics']['meanValueSFD']

    ebvalue = float(_ebv.split('(mag)')[0])
    if ebvalue<thresh:
        if verbose:
            print("ra=%.2f dec=%.2f:\tebv=%.2f\tOK"%\
                  (ra,dec,float(_ebv.split('(mag)')[0])))
        return ebvalue
    else:
        if verbose:
            print("ra=%.2f dec=%.2f:\tebv=%.2f\tNo"%\
                  (ra,dec,float(_ebv.split('(mag)')[0])))
        return

def query_vizier(_catname,_size,_interactive,ramin,ramax,decmin,decmax,absmag,distmin,distmax,verbose=False):
    """
    Download galaxy catalog from vizier
    return calalog and its name
    """    

    # find catfile from catalog name
    catalog_list = Vizier.find_catalogs(_catname)    
    if len(catalog_list)==0:sys.exit('Error: unavalable catalogs')
        
    # specify
    if _catname=='GLADE':
        _raname,_decname,_magname,_distname,_columns,_cat = 'RAJ2000', 'DEJ2000', 'BMAG', 'Dist',\
                                                            ['PGC', 'GWGC', 'HyperLEDA', '2MASS', 'RAJ2000', 'DEJ2000', 'BMAG', 'Dist'],\
                                                            'VII/281'
    elif _catname=='GWGC':
        _name,_raname,_decname,_magname,_distname,_columns,_cat = 'Name','RAJ2000', 'DEJ2000','BMAG','Dist',\
                                                            ['Name','RAJ2000', 'DEJ2000','BMAG','Dist'],\
                                                            'VII/267'
    elif _catname=='NED':
        _name,_raname,_decname,_magname,_distname,_columns,_cat = 'MGC','RAJ2000', 'DEJ2000','Bmag','z',\
                                                            ['MGC','RAJ2000', 'DEJ2000','Bmag','z'],\
                                                            'VII/240'
    else:sys.exit('Error: outside galaxy catalogs')
 
    if _magname == 'Bmag':
        # for apparent mag, to abs mag
        #if readkeys.keys(catname)[2] == 'Bmag':mag00 = mag00-5*np.log10(dist00)-25   
        sys.exit('Need TBD for cat with apparent magnitude, now choose another cat with absolute magnitude!!!')

    # download catalog with vizier
    '''
    quite strange vizier:
    if one two conditions for dec selection,
    only the first one works!!!
    so, I set only one condition on dec here and use numpy to set for the second
    '''

    v = Vizier(columns=_columns,
               column_filters={_magname:'<'+str(absmag),\
                               _distname:'%s..%s'%(str(distmin),str(distmax)),\
                               _raname:'%s..%s'%(str(ramin),str(ramax)),\
                               _decname:'%s..%s'%(str(decmin),str(decmax))})
    v.ROW_LIMIT = _size
    catalogs = v.get_catalogs(_cat)[0]    

    if verbose: 
        print("%i galaxies selected from %s in %i sec"%\
              (len(catalogs),_catname,int(time.time()-start_time)))

    # return infos    
    if _catname=='GLADE':
        _name = []
        for ii in range(len(catalogs)):_name.append('%s:%s:%s:%s'%(catalogs['PGC'][ii], \
                                                                   catalogs['GWGC'][ii], \
                                                                   catalogs['HyperLEDA'][ii], \
                                                                   catalogs['_2MASS'][ii]))        
    else:_name = catalogs[_name]
    return _name,np.array(catalogs[_raname]),np.array(catalogs[_decname]),np.array(catalogs[_magname]),np.array(catalogs[_distname])
