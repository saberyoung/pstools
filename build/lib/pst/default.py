"""PStools default settings"""

''' contents to show'''
# - from XML:
strxml = ['GraceID', 'AlertType', 'Group', 'FAR', \
          'Terrestrial', 'HasNS', 'HasRemnant', \
          'BNS', 'BBH', 'NSBH', \
          'Instruments', 'EventPage' ]
# - from fits header:
strhdr = ['DISTMEAN', 'DISTSTD', 'DATE-OBS', \
          'MJD-OBS']

''' healpix resolution'''
# default nside of priorization map
# will change according to trigger healpix map if any
nside = 1024 

# for computer contour lines
# decrease nside in order to be faster
nside1 = 64

''' figure setting'''
# default figure size
figsize = '15,10'

# 2D contours to show for triggers
# divided by comma, varied from 0 to 1
contours = [.5,.9]

# define colors for above contours
# number should be >= contours
color_contour = ['b','g','k','y','c','m']
                
# define colors for tiling/galaxies
color_field = ['b','g','k','y','c','m']

# default paramater for healpix view
# healpix ordering
ordering = False

#default paramater for healpix view
# theta = pi/2 - decl
# unit in deg
theta = 0

#default paramater for healpix view
# phi = ra
# unit in deg
phi = 0

#default paramater for healpix view
# coordinate system
# options: C, E, G
coord = 'C'

#default paramater for healpix view
# normalization
# options:
#  hist: histagram
#  log: logarithmic
#  None: linear
norm = 'None'

# minimum range value for 2d normalized map
minv = 0

# maximum range value for 2d normalized map
maxv = 1e-4

# healpix title
title = '2D sky map'

''' react setting'''
#when monitoring, for the role of triggers,
# if True: respond also to test alerts
# if False: respond only to observation alerts
test = False

''' galaxy construction'''

# size for querying the galaxy catalog
# 100 for 100 lines of catalog
# -1 for the full catalog
size = -1

# cache mode for galaxy and tiling generation: 
# 1. read galaxies/tilings from $cachefile (if setted and exists), see if cashed galaxies/tilings could also meet the galaxy/tiling cuts (i.e. ra, dec, mag, fov, ...), if everything works, read cached galaxies/tilings and stored them into $cachefile, otherwise, clobber $cachefile, query Vizier/generate tilings, and store to a new $cachefile; 
# 2. the same as 1, but instead not store to a new $cachefile; 
# 3. always query Vizier/generate tilings, clobber exiting $cachefile if any, and then store yields to a new $cachefile; 
# 4. the same as 3, but instead not store to a new $cachefile
cachemode = 1

# file name in which provides pre-defined galaxies. Since npz file would be used, if no .npz included in the file suffix, a .npz will be put in the end. If blanket will skip
cachefileg = 'tmp_pst_galaxies.npz'

# file name in which provides pre-defined tilings. Since npz file would be used, if not .npz included in the suffix, a .npz will be put in the end. leave blanket will skip
cachefilet = 'tmp_pst_%s.npz'

# in order to save time, set a number that will calculate visibility every a number of fields, only when trigger is not huge so that they are in the same portion of sky. If not int given, will calculate visibility for each fields, namely, nfields=1
nfields = 1

# remain a number of fields, otherwise, sometimes there're so many fields, especially for galaxy case, will slow down the process a lot
remfields = 1000
