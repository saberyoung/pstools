'''
Database scripts
'''
import sys
import pymysql

####################################################
def getconnection(database):
   cc = {}
   cc["schmidt"] = {}
   cc["schmidt"]['db'] = 'gw'
   cc["schmidt"]['host'] = 'localhost'
   cc["schmidt"]['user'] = 'syang'
   cc["schmidt"]['passwd'] = 'augusto90'
   return cc[database]

def dbConnect(cc):

   try:
      conn = pymysql.connect (host =       cc['host'],
                              user =       cc['user'],
                              passwd =     cc['passwd'],
                              db =         cc['db'],
                              autocommit = True)
   except pymysql.Error as e:
      sys.exit("Error %d: %s" % (e.args[0], e.args[1]))
   return conn

def insert_values(table,values):

   try:
      cc = getconnection('schmidt')
      connection = dbConnect(cc)     
   except:
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
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except pymysql.Error as e:
      print(("Error %d: %s" % (e.args[0], e.args[1])))
      sys.exit (1)  
   connection.close()

def update_value(table,column,value,filename,filename0='filename'):  

   try:
      cc = getconnection('schmidt')
      connection = dbConnect(cc)       
   except:
      sys.exit('!!! ERROR: problem with the mysql database !!!')

   try:
      cursor = connection.cursor (pymysql.cursors.DictCursor)
      if value in [True,False,'NULL',None]:
         cursor.execute ("UPDATE "+str(table)+" set "+column+"="+str(value)+" where "+str(filename0)+"= "+"'"+str(filename)+"'"+"   ")
      else:
         cursor.execute ("UPDATE "+str(table)+" set "+column+"="+"'"+str(value)+"'"+" where "+str(filename0)+"= "+"'"+str(filename)+"'"+"   ")
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except pymysql.Error as e:
      print(("Error %d: %s" % (e.args[0], e.args[1])))
      sys.exit (1)  
   connection.close()

def delete_values(archive,filename,column):

   try:
      cc = getconnection('schmidt')
      connection = dbConnect(cc)      
   except:
      sys.exit('!!! ERROR: problem with the mysql database !!!')

   try:
      cursor = connection.cursor (pymysql.cursors.DictCursor)
      cursor.execute ("delete  from "+str(archive)+" where "+str(column)+"="+"'"+filename+"'")
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except pymysql.Error as e:
      print(("Error %d: %s" % (e.args[0], e.args[1])))
      sys.exit (1)
   connection.close()
   return resultSet

def query(command):

   try:
      cc = getconnection('schmidt')
      connection = dbConnect(cc)     
   except:
      sys.exit('!!! ERROR: problem with the mysql database !!!')

   lista=''   
   try:
      cursor = connection.cursor (pymysql.cursors.DictCursor)
      for i in command:          
         cursor.execute(i)
         lista = cursor.fetchall ()
         if cursor.rowcount == 0:
            pass
         cursor.close ()
   except pymysql.Error as e:        
      print(("Error %d: %s" % (e.args[0], e.args[1])))
   
   connection.close()
   return lista
