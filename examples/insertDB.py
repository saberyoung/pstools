#!/usr/bin/env python
from __future__ import print_function
from builtins import input
import sys,pst

# define databse
cc = {}
cc['db'] = 'gw'
cc['host'] = 'localhost'
cc['user'] = 'syang'
cc['passwd'] = 'augusto90'

# insert to specific table
_dict = eval(sys.argv[1])
print ('Insert %s to DB'%_dict)
pst.insert_values(cc,'ligoevents',_dict)
