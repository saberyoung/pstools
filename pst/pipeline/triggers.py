#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : pst/pipeline/triggers.py
# Author            : syang <saberyoung@gmail.com>
# Date              : 01.01.2019
# Last Modified Date: 23.12.2019
# Last Modified By  : syang <saberyoung@gmail.com>

import logging
import numpy as np
import healpy as hp
import astropy.time

__all__ = ('PstParseTriggers')

class PstParseTriggers():
    """
    * implement the following functions:
    -> __init__(self, hpmap=None, fithdr=None, xmlinf=None, defconf=None, logger = None)
    -> url/fits/xml/root/coo
    -> make_report
    -> calc_contours
    -> calc_area
    -> data_list
    -> set/get/resetall/reset
    """

    # Static version info
    version = 1.0
        
    def __init__(self, hpmap=None, fithdr=None,
                 xmlinf=None, defconf=None, logger = None):
        """
        Read trigger map from a given source
        only one needed, if more than 1, read the latest        

        Parameters:
        -----------
        hpmap       `list`           healpix fits for 3d localization
        fithdr      `dict`           healpix fits map informations
        xmlinf      `dict`           voevent xml file informations
        defconf     `dict`           default configure
        logger      `class`          logger
        """

        # ----- define logger ----- #
        if logger is None:
            logging.basicConfig(level = logging.INFO)
            self.logger = logging.getLogger(__name__)
        else:
            self.logger = logger

        # ----- set trigger properties ----- #
        self.data               =   {}
        self.data['hpmap']      =   hpmap
        self.data['xmlinf']     =   xmlinf
        self.data['fithdr']     =   fithdr

        # ----- define keys to check xml ----- #
        '''
        skymap_key  `string`   key for skymap url, in XML   
        radec_keys  `list`     key for coordinate and uncertainty, in XML
        obst_key    `string`   key for timeobs, in XML
        '''
        self.obst_key           =   'ISOTime'
        self.skymap_key         =   'skymap_fits'
        self.radec_keys         =   ['C1','C2','Error2Radius']
        
        # ----- define informations for report ----- #
        self.keys_to_checkxml = ('GraceID', 'AlertType', 'Group', 'FAR',
                                 'Terrestrial', 'HasNS', 'HasRemnant',
                                 'BNS', 'BBH', 'NSBH', 'Instruments',
                                 'EventPage')
        self.keys_to_checkfits = ('DISTMEAN', 'DISTSTD', 'DATE-OBS',
                                  'MJD-OBS', 'OBJECT', 'INSTRUME',
                                  'CREATOR')

        # ----- define default parameters ----- #
        self.run_config(defconf)

    def run_config(self, defconf):

        self.conf = {
            'wdir':      '/tmp/',   # `str`   working directory
            'savefits':  None,      # `str`   save file or not
            'nside':     512,       # `int`   healpix map default resolution
            'coord':     'C',       # `str`   healpix map defualt coordinate system
            'cls':       [.5,.9],   # `list`  confidence levels for contours
            'style':     'sms'      # `str`   report type
        }
        
        if defconf is None: return
        for _k in self.conf.keys():
            if _k in defconf.keys():
                self.conf[_k] = defconf[_k]
            else:
                self.logger.info ('### Warning: use default value for %s'%_k)

    '''
    parse sources, to obtain self.data
    Source options:
    -----------
        skymap_url  `string`         url of healpix fits
        skymap      `string`         filename of healpix fits
        xmlfile     `string`         the filename of XML
        root        `class`          root element of the XML document parsed by lxml.etree
        coolist     `list`           [ra,dec,err], `float`, `float`, `float`   
    '''
    
    def url(self, url, wdir=None, savefits=None):
        """ 
        apply parser on skymap url        
        """
        
        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        
        # parse fits_url
        if not url is None: self.fits(url)

        # savefits
        if not self.conf['savefits'] is None:
            flag = self.download_skymap(self.conf['wdir'], self.conf['savefits'], url)
            if flag == 1:
                self.logger.info ('download skymap via wget')
            elif flag == 2:
                self.logger.info ('download skymap via requests')
            elif flag == 3:
                self.logger.info ('### Warning: failed to download skymap,'+\
                                  'install wget or requests first')
            
    def fits(self, skymap):
        """ 
        apply parser on healpix fits or fits url
        """
        
        _hp = self.checkhp()
        if _hp: return

        if not skymap is None:            
            hpmap, fitsinfos = self.read_skymap(skymap)
            self.data['hpmap'] = hpmap            
            self.data['fithdr'] = fitsinfos

    def xml(self, xml, wdir=None, savefits=None, nside=None, coord=None):
        """ 
        apply parser on voenevt XML file
        """

        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        if not nside is None: self.conf['nside'] = nside
        if not coord is None: self.conf['coord'] = coord
        
        # read xml file, obtain root obj        
        root = self.parse_xml(xml)

        # parse root
        self.root(root, wdir=self.conf['wdir'],
                  savefits=self.conf['savefits'],
                  nside=self.conf['nside'],
                  coord=self.conf['coord'] )
            
    def root(self, root, wdir=None, savefits=None, nside=None, coord=None):
        """ 
        apply parser on root element of the XML document, parsed by lxml.etree
        """

        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        if not nside is None: self.conf['nside'] = nside
        if not coord is None: self.conf['coord'] = coord
        
        _hp = self.checkhp()
        if _hp: return
        
        if not root is None:
            
            self.data['xmlinf'] = {elem.attrib['name']:
                                   elem.attrib['value']
                                   for elem in root.iterfind('.//Param')}   
            
            # 1. find ra,dec,error box: for AMON, GRB, etc   
            coo = self.parse_root_coo(root, self.conf['nside'])
            if not coo is None:
                self.coo(coo, wdir=self.conf['wdir'],
                         savefits=self.conf['savefits'],
                         nside=self.conf['nside'],
                         coord=self.conf['coord'])

            # 2. find skymap url: for GW
            url = self.parse_root_skymap(root)
            if not url is None:
                self.url(url, wdir=self.conf['wdir'],
                         savefits=self.conf['savefits'])
            
    def coo(self, coo, wdir=None, savefits=None, nside=None, coord=None):
        """ 
        apply parser on ra, dec, error box
        obtain hpmap
        """

        if not wdir is None: self.conf['wdir'] = wdir
        if not savefits is None: self.conf['savefits'] = savefits
        if not nside is None: self.conf['nside'] = nside
        if not coord is None: self.conf['coord'] = coord
        
        _hp = self.checkhp()
        if _hp: return

        _ra,_dec,_loc = coo
        _radius = _loc*2*np.pi/360 # from deg to radians
        _pmap = np.zeros(hp.nside2npix(self.conf['nside']))
        _index = hp.pixelfunc.ang2pix(self.conf['nside'],
                        np.radians(-dec+90.), np.radians(_ra))
        _pmap[_index]+=1
        _pmap=hp.sphtfunc.smoothing(_pmap,fwhm=_radius)                
        self.data['hpmap'] = _pmap/sum(_pmap)

        # savefits
        if not self.conf['savefits'] is None:
            flag = self.make_hpmap(self.conf['wdir'],
                                   self.conf['savefits'],
                                   self.conf['coord'])
            if flag:
                self.logger.info ('make healpix map')            
            else:
                self.logger.info ('### Warning: failed to make healpix map')
                                  
    def checkhp(self):
        if not self.data['hpmap'] is None:
            self.logger.info ('hpmap has already been parsed')
            return True
        else:
            return False
        
    def make_hpmap(self, wdir, tempfile, coord):
        # generate map       
        
        skymap = '%s/%s' % (wdir, tempfile)
        try:
            hp.write_map(skymap, self.data['hpmap'], coord=coord,
                         extra_header=self.data['fithdr'], overwrite=True)
            return True
        except:
            return False

    def parse_root_skymap(self, root):
        """
        obtain healpix fits map      via skymap url
        """

        # skymap url
        xmlinfos = {elem.attrib['name']:
                    elem.attrib['value']
                    for elem in root.iterfind('.//Param')}        
        if self.skymap_key in xmlinfos:
            self.logger.info ('obtain hpmap via %s' % self.skymap_key)
            skymap_url = xmlinfos[self.skymap_key]
        else:
            self.logger.info ('### Warning: no %s found in voevent'%self.skymap_key)
            skymap_url = None
        return skymap_url
            
    def parse_root_coo(self, root, nside):
        """
        obtain ra,dec,loc
        """
        
        _coolist,_coo = [],0
        for _kk in self.radec_keys:
            _val = root.findtext('.//%s'%_kk)
            if not _val is None:
                try:
                    _coolist.append(float(_val))
                    _coo += 1
                except:
                    self.logger.info ('### Warning: %s not float, skip'%_kk)
            else:
                self.logger.info ('### Warning: no %s found in voevent'%_kk)
                
        if _coo == 3:
            _ra,_dec,_loc = _coolist
        else:
            return None
        
        """
        make header            
        """        
        # obs time
        timeobs = root.findtext('.//%s'%self.obst_key)

        # obs mjd        
        mjdobs = astropy.time.Time(timeobs, format='isot', scale='utc').mjd

        # read classifier
        try:  tid = root.attrib['ivorn']
        except:  tid = None

        # header
        hdr = [('CREATOR','PSTOOLS'),
               ('IVORN',tid),
               ('RA',_ra),
               ('DEC',_dec),
               ('ERROR-BOX',_loc),
               ('NSIDE',nside),
               ('MJD-OBS',mjdobs),
               ('DATE-OBS',timeobs)]
        
        self.data['fithdr'] = dict(hdr)
        return _coolist
    
    def make_report(self, style=None):
        """
        make a summary report on input source
        """
        
        if not style is None: self.conf['style'] = style
        
        _dict =  {}
        for el in self.keys_to_checkxml:
            try:
                _dict[el] = self.data['xmlinf'][el]
            except:
                self.logger.info ('### Warning: no %s found in voevent'%el)
        for el in self.keys_to_checkfits:
            try:
                _dict[el] = self.data['fithdr'][el]
            except:
                self.logger.info ('### Warning: no %s found in fits header'%el)

        _report = ''
        for el in _dict:
            if self.conf['style'] == 'sms':
                _report += '%s:%s '%(el,_dict[el])
            elif self.conf['style'] == 'slack':
                _report += '`%s`:%s\n'%(el,_dict[el])
            elif self.conf['style'] == 'email':
                _report += '#\t%s:\t%s\n'%(el,_dict[el])           

        if not self.conf['style'] in ['sms', 'email', 'slack']:
            self.logger.info ('### Warning: style %s not available'%self.conf['style'])
            self.logger.info ('### select from email, slack, sms')
        return _report

    def calc_contours(self, cls=None):
        '''
        calculate index in specific confidence levels
        '''
        
        if not cls is None: self.conf['cls'] = cls
        
        _hp = self.checkhp()
        if not _hp:
            self.logger.info ('### Warning: no hpmap found')
            return

        try: (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        except: hpx = self.data['hpmap']
        
        indexlist = {}
        sortid = hpx.argsort()[::-1]
        sorted_hpx = hpx[sortid]
        cumulative_hpx = np.cumsum(sorted_hpx)
        for _cl in self.conf['cls']:
            if len(sorted_hpx[cumulative_hpx<_cl]) > 0:
                _limit = sorted_hpx[cumulative_hpx<_cl][-1]
                indexlist[_cl] = [i for i in sortid if hpx[i] >= _limit]
        return indexlist
    
    def calc_area(self, cls=None):
        '''
        calculate uncertainty region area, unit in sq. deg
        '''

        if not cls is None: self.conf['cls'] = cls
        
        _ilist = self.calc_contours(cls=self.conf['cls'])
        if _ilist is None: return
        
        try: (hpx, hpd1, hpd2, hpd3) = self.data['hpmap']
        except: hpx = self.data['hpmap']

        nside = hp.get_nside(hpx)    
        areapix = (hp.nside2resol(nside,arcmin=True)/60.)**2
        _alist = {}
        for _cl in self.conf['cls']:
            _area = areapix*len(_ilist[_cl])
            _alist[_cl] = _area
        return _alist    

    def parse_xml(self, xmlfile):
        try:
            from lxml import etree
            self.logger.info ("running with lxml.etree")
        except ImportError:
            try:
                # Python 2.5
                import xml.etree.cElementTree as etree
                self.logger.info ("running with cElementTree on Python 2.5+")
            except ImportError:
                try:
                    # Python 2.5
                    import xml.etree.ElementTree as etree
                    self.logger.info ("running with ElementTree on Python 2.5+")
                except ImportError:
                    try:
                        # normal cElementTree install
                        import cElementTree as etree
                        self.logger.info ("running with cElementTree")
                    except ImportError:
                        try:
                            # normal ElementTree install
                            import elementtree.ElementTree as etree
                            self.logger.info ("running with ElementTree")
                        except ImportError:
                            self.logger.info ("Failed to import ElementTree from any known place")
                            return None
        tree = etree.parse(xmlfile)
        root = tree.getroot()        
        return root

    def data_list(self):
        return self.data
            
    def set(self, _k, _v):
        self.data[_k] = _v
    
    def get(self, _k):
        return self.data[_k]

    def resetall(self):
        self.__init__()

    def reset(self, _k):
        self.data[_k] = None
        
    @staticmethod
    def read_skymap(_fits):
        # get healpix fits        
        try: # read 3D trigger healpix map (for GW)            
            tmap, header = hp.read_map(_fits, field=[0, 1, 2, 3],h=True)
        except: # if failed, read 2D map            
            try:
                tmap, header = hp.read_map(_fits, h=True)                
            except:
                tmap, header = None, None
        return tmap, dict(header)

    @staticmethod
    def download_skymap(wdir, tmpfile, url):
        """
        Look up URL of sky map, download sky map, and parse FITS file.
        """

        skymap = '%s/%s' % (wdir, tmpfile)        
        try: 
            import wget            
            if os.path.exists(skymap): os.remove(skymap)
            wget.download(url, out=skymap)
            
            flag = 1
        except:
            try:
                import requests,tempfile,shutil                

                # Send HTTP request for sky map   
                response = requests.get(url, stream=True)   

                # Raise an exception unless the download succeeded (HTTP 200 OK)
                response.raise_for_status()

                # Create a temporary file to store the downloaded FITS file
                with tempfile.NamedTemporaryFile() as tmpfile:
                    # Save the FITS file to the temporary file
                    shutil.copyfileobj(response.raw, tmpfile)
                    tmpfile.flush()

                    # Uncomment to save FITS payload to file       
                    shutil.copy(tmpfile.name, skymap)
                    
                    flag = 2
            except:
                flag = 3
        return flag
