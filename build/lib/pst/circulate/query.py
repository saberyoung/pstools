import numpy as np
import sys
import time

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

def query_vizier(_catalog,_size,filtro,limra,limdec,limmag,limdist,verbose=False):
    """
    Download galaxy catalog from vizier
    return calalog and its name
    """    
    start_time = time.time()

    # specify columns
    if _catalog == 1:
        if not filtro in ['B','K']:sys.exit('### Error: wrong filters')
        _catid, _catname, _columns = 'VII/281', 'GLADE', \
                ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                 'Dist', 'PGC', 'GWGC', 'HyperLEDA', '2MASS']                                                            
    elif _catalog == 2:
        if not filtro in ['B']:sys.exit('### Error: wrong filters')
        _catid, _catname, _columns = 'VII/267', 'GWGC', \
                ['RAJ2000', 'DEJ2000', '%sMAG'%filtro, \
                 'Dist', 'Name']
    else:sys.exit('### Error: wrong galaxy catalogs')
 

    # download catalog with vizier
    try: from astroquery.vizier import Vizier
    except:sys.exit('### Error with vizier')

    v = Vizier(columns=_columns, \
               column_filters={_columns[0]:'%s..%s'%(str(limra[0]),str(limra[1])),\
                               _columns[1]:'%s..%s'%(str(limdec[0]),str(limdec[1])),\
                               _columns[2]:'%s..%s'%(str(limmag[0]),str(limmag[1])),\
                               _columns[3]:'%s..%s'%(str(limdist[0]),str(limdist[1]))
                           })                               
    v.ROW_LIMIT = _size
    catalogs = v.get_catalogs(_catid)[0]

    if verbose: 
        print("%i galaxies selected from %s in %i sec"%\
              (len(catalogs),_catid,int(time.time()-start_time)))

    # return infos    
    if _catalog == 1:
        _name = []
        for ii in range(len(catalogs)):
            _name.append('%s:%s:%s:%s'%(catalogs[_columns[4]][ii], \
                                        catalogs[_columns[5]][ii], \
                                        catalogs[_columns[6]][ii], \
                                        catalogs['_%s'%_columns[7]][ii]))        
    else:_name = catalogs[_columns[4]]
    return np.array(_name),np.array(catalogs[_columns[0]]),\
        np.array(catalogs[_columns[1]]),np.array(catalogs[_columns[2]]),\
        np.array(catalogs[_columns[3]])
