#!/usr/bin/env python
from __future__ import print_function
from builtins import input
import sys,pst

# insert to specific table
_ff = sys.argv[1]

# include ra,dec,filter
_list           =   open(_ff).readlines() 

_id             =   sys.argv[2] # identifier
airmass         =   sys.argv[3] # airmass limit
dither          =   sys.argv[4] # for VST, currently, use fixed dither
nexp            =   sys.argv[5] # number of exposure per field,
exptime         =   sys.argv[6] # exposure time

airmass,dither,nexp,exptime = float(airmass),\
            float(dither),int(nexp),float(exptime)

# VST scheduler
if False:
    environment =  'production' # VST environment
    username    =  'saberyoung' # VST username
    password    =  'Ys19900615' # VST password
    containerid =  2229852    # VST containerid
else:
    environment =  'demo'       # VST environment
    username    =  '52052'      # VST username
    password    =  'tutorial'   # VST password
    containerid =  1455705    # VST containerid

# default VST setting
userpriority    =   1           # OB priority
skytransparency =  '3'          # sky transparency:
                                # 1 - Photometric
                                # 2 - Clear
                                # 3 - Variable, thin cirrus
                                # 4 - Variable, thick cirrus
fli             =   1.0         # fractional_lunar_illumination
seeing          =   2.0         # seeing limit

#
_weburl = pst.VST(_list,_id,environment,username,\
	password,containerid,userpriority,airmass,\
	skytransparency,fli,seeing,exptime,\
	dither,nexp)
print('VST scheduler done\tcheck VST OB in p2 webpage, %s'%_weburl)
