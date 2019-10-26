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
nexp            =   sys.argv[5] # number of exposure per field
exptime         =   sys.argv[6] # exposure time

airmass,dither,nexp,exptime = float(airmass),\
            float(dither),int(nexp),float(exptime)

_ofname         =   'schmidtOB.txt'

#
pst.schmidt(_list,_ofname,nexp,exptime,dither,_id)
print('schmidt scheduler done\tcheck schmidt OB in %s'%_ofname)
