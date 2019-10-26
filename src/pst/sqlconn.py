'''
Database scripts
'''
from __future__ import print_function
from builtins import input
import sys, pymysql

####################################################
def dbConnect(cc):
   try:
      conn = pymysql.connect (host =       cc['host'],
                              user =       cc['user'],
                              passwd =     cc['passwd'],
                              db =         cc['db'],
                              autocommit = True)
   except pymysql.Error as e:
      print ("Error %d: %s" % (e.args[0], e.args[1]))
      conn = False
   return conn

def query(cc,command):

   connection = dbConnect(cc)     
   if not connection:
      sys.exit('!!! ERROR: problem with the mysql database !!!')

   lista=''   
   try:
      cursor = connection.cursor (pymysql.cursors.DictCursor)
      for i in command:          
         cursor.execute(i)
         lista = cursor.fetchall ()
         if cursor.rowcount == 0: pass
         cursor.close ()
   except pymysql.Error as e:        
      print(("Error %d: %s" % (e.args[0], e.args[1])))
   
   connection.close()
   return lista

def insert_values(cc,table,values):

   connection = dbConnect(cc)
   if not connection:
      sys.exit('!!! ERROR: problem with the mysql database !!!')

   def dictValuePad(key):
      return '%(' + str(key) + ')s'

   def insertFromDict(table, dict):
      """Take dictionary object dict and produce sql for 
      inserting it into the named table"""
      sql = 'INSERT INTO '+ table
      sql += ' ('
      sql += ', '.join(dict)
      sql += ') VALUES ('
      sql += ', '.join(map(dictValuePad, dict))
      sql += ');'
      return sql

   sql = insertFromDict(table, values)

   try:
      cursor = connection.cursor (pymysql.cursors.DictCursor)
      cursor.execute(sql, values)
      resultSet = cursor.fetchall()
      if cursor.rowcount == 0: pass
      cursor.close ()
   except pymysql.Error as e:
      print(("Error %d: %s" % (e.args[0], e.args[1])))
      sys.exit (1)
   connection.close()
