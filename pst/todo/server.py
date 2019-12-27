from __future__ import print_function
from builtins import input

# https://gcn.gsfc.nasa.gov/voevent.html
def lvc_server(_url,_priv,_verbose):
    if _url == 'local':
        if _verbose: 
            print('!!!Only public port for %s'%_url)
        _server,_port = '127.0.0.1',8099

    if _url == 'eApps':
        _server = '68.169.57.253'
        if _priv == 'PUBLIC':_port = 8099
        elif _priv == 'PRIVATE':
            answ = input('1.LVC or 2.AMON')
            if answ == '1': _port = 8092
            elif answ == '2': _port = 8096
            else:sys.exit('Error:wrong')

    if _url == 'Atlantic_2':
        _server = '209.208.78.170'
        if _priv == 'PUBLIC':_port = 8099
        elif _priv == 'PRIVATE':
            answ = input('1.LVC or 2.AMON')
            if answ == '1': _port = 8092
            elif answ == '2': _port = 8096
            else:sys.exit('Error:wrong')

    if _url == 'Atlantic_3':
        if _verbose:
            print('!!!Only public port for %s'%_url)
        _server,_port = '45.58.43.186',8099

    if _url == 'Linode':
        _server = '50.116.49.68'
        if _priv == 'PUBLIC':_port = 8099
        elif _priv == 'PRIVATE': _port = 8096

    return _server,_port
