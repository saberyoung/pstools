'''PlateSystem
Define a platesystem, a collection fields ("plates") of a constant size
covering the entire sky.
Edwin Valentijn - 7 March 2004 version with smooth grid - more advanced than PAE

improved negative dec performance - plotting
'''
__version__ = '@(#)$Revision$'

import math

def sum(args):
   '''The sum function is ONLY built in in Python 2.3 and newer'''
   import operator
   return reduce(operator.add, args, 0)

class PlateSystem:
    '''Define a plate system

    The plate system is defined by a set of fields covering the sky.
    Each field has a size, a centre and a minimum amount of overlap with
    neighbouring fields.
    The fields are arranged in strips perpendicular to the meredian.
    This means that with increasing declination there are less fields per strip.
    Script optimizes for contnuity form one dec strip to the next;
    therfor the positions in various dec stips are dependent on each other.
    '''
    def __init__(self, fieldsize=1.0, overlap=5, efflimit=4):
        '''fieldsize - the size of the output fields in degrees
        overlap      -  % fraction of overlap between fields
        efflimit     -  % fraction minimum overlap
        '''
        self.fieldsize      = fieldsize                             # 1.0
        self.efflimit       = efflimit                              # 4
        self.overlap        = 1. +overlap/100.                      # 1.05
        self.distance       = fieldsize /  self.overlap             # 1/1.05
# make dec grid first
        numberstrips        = int(90.0 / self.distance)
                                                       # make seq of dec centers
        dec                 = range (numberstrips)
        for j in range (0, numberstrips):
            dec[j]    = j * 90.0/numberstrips                  
        self.numberstrips   = numberstrips +1                       # add polar field
        dec.append(89.9)                                            # add polar field 
        self.dec            = dec                                   # key grid of stips in de
#       print self.dec        
        self.n,self.e       = self.get_numberra_centers(0)
#
#  solution for numberoffields in RA for each dec strip
#
    def get_numberra_centers(self,log): # determine number of fields for each dec
        '''
        log   - print log to screen '''
        nold     = 400
        dold     = 0
        dec      = self.dec 
        n =  range (self.numberstrips)
        e =  range (self.numberstrips)
        for i in range (0, self.numberstrips):
           smallcircle      = 360.0 * math.cos( self.dec[i]* math.pi / 180.)
           numberperstrip   = int(round(smallcircle / self.distance))    #better int(round())
           effectiveoverlap = 100*( (self.fieldsize/(smallcircle/numberperstrip)) -1.)
           if log:
             print "--------------------------------------------------------------------"
             print "rawoverlap=", "%5.1f" %effectiveoverlap
                                            # finetune on grid 1
           while effectiveoverlap < self.efflimit :       # handles both negative and limit on overlap
                numberperstrip  +=1
                effectiveoverlap = 100*( (self.fieldsize/(smallcircle/numberperstrip)) -1.)
                if log: print "plusone overlap"
                                            # finetune on grid 2
                                            # smooth out the jumps in strips
                                            # CREATES SMOOTH GRID  but  not for pole
           if (nold-numberperstrip)<dold and dec[i]<84 :
                numberperstrip   = nold - dold
                effectiveoverlap = 100*( (self.fieldsize/(smallcircle/numberperstrip)) -1.)
                if log:
                   print "plusone smooth"
           if log:
                print "                                                         delta",(nold -numberperstrip)
                print i+1, "dec=", "%5.2f" %self.dec[i], "area=" ,"%5.2f" %smallcircle,\
                       "numberperstrip=",numberperstrip, "overlap in %=", "%5.1f" %effectiveoverlap
           if i > 0 :
                dold = nold-numberperstrip 
           nold = numberperstrip
           n[i] = numberperstrip
           e[i] = effectiveoverlap
        if log:
           print "PlateSystem: total number found=",sum(n)
        self.n = n
        self.e = e
        return n,e
#
    def nearest_dec(self,decin):
        ''' get nearest fieldcenter dec  for any dec'''
        dec  = self.dec
        decin = abs (decin)  # allow pass negative dec
        diff = range( len( dec))
        for i in range (0,len( dec)):
            diff[i] = abs( decin-dec[i])
        j = diff.index(min(diff))
        return dec[j],j
#        
    def get_ra(self,decin):                 # get seq of ra's for given decin
        ''' get seq of ra's for a single decin'''
        decplate,j = self.nearest_dec(decin)# get numberperstrip and index
        n = self.n
        ras             = range (n[j])        
        offsetra        = 360./n[j]
        for i in range (0,n[j]):
            ras[i] = i * offsetra
        return ras,decplate
#    
    def get_field_centre(self,rain,decin):        
        ''' get nearest ra fieldcenters, needs dec too'''
        if decin < 0:
            decsign = -1
        else:
            decsign = 1
        dec = abs(decin)
        if dec > 90.0:
            print" dec >90 we can't go there"
            dec = 90.
        ras, decplate = self.get_ra(dec)
        diff = range(len(ras))       
        for i in range(0, len(ras)):
            diff[i] = abs(rain-ras[i])
        j = diff.index(min(diff))
        rasplate = ras[j]
        return rasplate, decsign*decplate
#    
    def print_grid(self):
        ''' print grid properties in file'''
        print "Version with SMOOTHING of GRID " 
        print "sequence     Dec      number   overlap"
        print "         [degrees]             [%]"       
        for i in range (0, self.numberstrips):
           print "%8i" %(i+1), "%9.2f" %self.dec[i],"%9i" %self.n[i], \
                 "%9.1f" %self.e[i]
        print  "PlateSystem: total number of fields=", sum(self.n)
        return
# plot
# input parameters - range on sky grid -for plot only
    def plot_it(self,rastart=-1.,raend=360.,decstart=-1.,decend=90.):
        ''' plots the grid in range
        can be used as survey preparation tool
        selects fields within rectangle on sky
        returns ra,dec of selected fields in self
        autoscaling implemented                   ''' 

        import dislin

        rastart = float(rastart)
        decstart= float(decstart)
        raend   = float(raend)
        decend  = float(decend)
        selected_fields_ra = range(50000)
        selected_fields_dec= range(50000) 
        if abs(decend-decstart)<1.0 :
            decend  +=0.5
            decstart-=0.5
        if abs(raend-rastart)<1.0:
            rastart -=0.5
            raend   +=0.5
        symsize =   12
        symtype =    8
        xlen    = 2200
#
        #dislin.metafl ('tiff')      # outfile format ; befor disini 'post'=ps 'pdf'
        #dislin.metafl ('post')      # outfile format ; befor disini 'post'=ps 'pdf'<-----
        #dislin.metafl ('pscl')      # postscript in color; drak background
                                     # produces erroneous float bounding box
        dislin.metafl('cons')        # cons on screen (fills screen ***)
     #   dislin.metafl('xwin')         #        <------------- lines to change for ps output
        #dislin.setfil('C:\myPy\skygridpole.tiff')   # default is dislin
        #dislin.setfil('C:\myPy\skygridpole.ps')     # default is dislin
        #dislin.setfil('C:\myPy\skygridstrip.tiff')  # default is dislin
#        dislin.setfil('C:\myPy\skygridstrip2005.ps')    # default is dislin   <-----
        dislin.disini ()
        dislin.hwfont ()             # standard hardware fonts
        #?dislin.page(6400,200)         # ?? does it do anything ??
        dislin.pagera ()             # plot border
        dislin.axspos (450, 1800)    # lower left corner
        dislin.axslen (xlen, 1400)   # length height of plot
        dislin.name   ('RA degrees' ,'X')
        dislin.name   ('Dec degrees','Y')
        dislin.digits (-1, 'XY')     # -1 integer
        dislin.ticks  (3, 'X')       # number of ticks between labels
        dislin.ticks  (10, 'Y')      # number of ticks between labels
        dislin.texmod ('ON')
        dislin.titlin ('OmegaCam \it Sky Grid  1 degree', 1)
        ticksx = 20
        if raend-rastart<100: ticksx = 10     # autoscaling very premature
        if raend-rastart<10: ticksx = 1      
        ticksy = 10.
        if abs(decend-decstart) < 15: ticksy = 1.
        dislin.graf   (rastart,raend, 10*int(rastart/10), ticksx,\
                       decstart, decend, 10*int(decstart/10), ticksy) # axi
        dislin.title  ()             # print title defined in titlin - up to 4 lines
        dislin.color  ('blue')
        dislin.hsymbl (symsize)           # symbol size
#
        y    = self.dec
        ymin = range (len(y)-1)         # append negative dec s
        for ii in range(len(ymin)): ymin[ii] = -1 * y[ii+1]
        ymin.extend(y)
        ymin.sort()
        y = ymin                        # append done
        f  =  self.fieldsize
        count = 0 
        dc = 1400/abs(decend-decstart) #pix/degree for plot coord (instead of user)
        rc = xlen/abs(raend-rastart)   #pix/degree
        for i in range (0, len(y)):
            y1  = y[i]
            cosdec= math.cos( y1 * math.pi / 180.)
            x,de   = self.get_ra(y1)
            for j in range (0, len(x)):
                if x[j]> rastart and x[j] < raend and y1 >decstart and y1<decend:
                    count += 1
                    selected_fields_ra[count-1] = x[j]
                    selected_fields_dec[count-1]= y1
                    dislin.rlsymb (symtype,x[j],y1)   #  symbols in user coord
                    dislin.color ('green')
                    dislin.rectan (int(450+( (x[j]-(rastart+(0.5/cosdec))) *f*rc)),\
                                   int(1800-((0.5+y1-decstart)*f*dc)), \
                                   int(f*rc/cosdec),\
                                   int(f*dc))  # rectangles in plotcoord
                    dislin.color  ('blue')
        self.selected_field_ra  = selected_fields_ra[0:count]
        self.selected_field_dec = selected_fields_dec[0:count]
        e = self.e
        e= e[0:len(e)-1]             # skip polar field for overlap computation
        x= 'PlateSystem: total number of pointings=  '+str(sum(self.n))+\
           ' selected pointings='+str(count)+ ' average overlap=  '\
          +str(sum(e)/len(e))[0:4]+'%'
        print x
        dislin.titlin(x,3)
        dislin.titlin('optimized continuity',2)
        dislin.title  ()
        dislin.dash   ()
        dislin.xaxgit ()
        dislin.disfin ()
        return 
#    
# ways of invoking the class
#                              make
if __name__ == '__main__':
    grid = PlateSystem()
#    a = grid.get_numberra_centers(1)      # makes full solution - print log
    a = grid.get_numberra_centers(0)       #  full solution without log
#   nearest dec for free input dec
    for i in range(1,9):
        print "dec fieldccenter", 9.*i, grid.nearest_dec(9.*i)
        print "ra- dec nearest ra =", 9.*i,'dec=',9.*i ,grid.get_field_centre(9.*i,9.*i)
#
# print the key table of the plate grid
    a= grid.print_grid()
# plotit        
    a= grid.plot_it()
# survey preparation - list of coordinates:
    print "example of preparing a survey"
    b = grid.plot_it(16.5*15,17*15,15,18)
    for i in range(0,len(grid.selected_field_ra)):
	       print i+1,grid.selected_field_ra[i],grid.selected_field_dec[i]
#
