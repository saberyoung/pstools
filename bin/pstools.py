"""############################################################################
2019/1/30 Start version 0.1. Only module avaivale for call..
2019/3/27 version 0.2. add pstool.py. Init cookbook
2019/6/22 version 0.3. add slack.py. Possible for multiple telescopes, however limits their FoV composition
2019/9/17 version 0.4. use class
""" ############################################################################
from __future__ import print_function
from builtins import input
import argparse,time,sys,os,pst

start_time = time.time()
description = ">> PStools main algorithm"

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description=description,\
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)   

    # mandatory
    parser.add_argument("mode",choices=[1,2,3,4,5],type=int,\
        help='PSTools mode: [1] Serve alert; [2] Monitor alert; '+\
        '[3] Manual input map; [4] Generate configure file; '+\
        '[5] Tutorial')

    # optional
    parser.add_argument("-f", dest="fits",\
        help='handle healpix fits file')      
    parser.add_argument("-t", dest="tel",\
        help='specific telescopes, if multiple divided with `,`'+\
        '\t options:%s'%str(pst.configure.tel_list.keys()))
    parser.add_argument("-x",dest='xml',\
        help='xml file that would be served locally')
    parser.add_argument("-s", dest='server',\
        default='eApps',choices=['local','eApps',\
        'Atlantic_2','Atlantic_3','Linode'],\
        help='GCN server')     
    parser.add_argument("-p", dest='port',\
        default='PUBLIC',choices=['PUBLIC','PRIVATE'],\
        help='GCN port')
    parser.add_argument("-v", dest="verbose",\
        default=False,action="store_true",\
        help='Enable task progress report')  
    parser.add_argument("-l", dest="log",\
        default=False,action="store_true",\
        help='Recording task progress report')

    # read parameters
    args = parser.parse_args()
    _mode = args.mode
    _log = args.log
    _verbose = args.verbose
    _server, _port = args.server, args.port

    if _log: # if record
        import logging
        if os.path.exists('pstools.log'):os.remove('pstools.log')
        logging.basicConfig(filename='pstools.log', \
            level=logging.INFO,format='%(asctime)s.%(msecs)03d '+\
            '%(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')
        for key in logging.Logger.manager.loggerDict:
            logging.getLogger(key).setLevel(logging.CRITICAL)
    
    if _mode == 1:    # - serve alert
        from gcn.cmdline import serve_main
        if args.server != 'local':
            sys.exit('Only local server can be served')
        if args.xml is None:
            sys.exit('Pls specify a XML file to serve')
        _s1,_s2 = pst.server.lvc_server(_server, _port, _verbose)
        _info = 'Serve Mode:\thost:%s:port:%s'%(_s1,_s2)
        if _verbose: print('## %s'%_info)
        if _log: logging.info(_info)
        sys.exit(serve_main(args=[args.xml,'--host', '%s:%s'%(_s1,_s2)]))

    elif _mode == 2:  # - auto listen
        import gcn
        _s1,_s2 = pst.server.lvc_server(_server, _port, _verbose)
        _info = 'Monitor mode: listening host:%s:port:%s'%(_s1,_s2)
        if _verbose: print('## %s'%_info)
        if _log: logging.info(_info)
        gcn.listen(handler = pst.autogcndef.process_gcn, host=_s1, port=_s2)

    elif _mode == 3:  # - manual search        
        _tel = args.tel
        if _tel is None:sys.exit('Error: option tel needed for man mode..')
        pst.manschedule.man_search(args.fits, _tel)

    elif _mode == 4:  # - configure file
        _tel = args.tel
        if _tel is None:sys.exit('Error: option tel needed for man mode..')
        _tel = _tel.split(',')
        for _teli in _tel: pst.configure.config_init(_teli)        

    elif _mode == 5:  # - web tutorial
        if _verbose: print ('https://pstool-cookbook.readthedocs.io/en/latest/')
        _html = '%s/../../docs/_build/html/index.html'%pst.__path__[0]
        import webbrowser       
        webbrowser.open(_html, new=2)

print("-"*80)
print('>>>>> Completed in '+str(int(time.time()-start_time))+' sec\n')
