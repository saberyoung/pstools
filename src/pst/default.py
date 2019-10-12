"""PStools default settings"""

''' contents to show'''
# - from XML:
strxml = ['GraceID', 'AlertType', 'Group', 'FAR', \
           'Terrestrial', 'HasNS', 'Instruments', 'EventPage' ]
# - from fits header:
strhdr = ['DISTMEAN', 'DISTSTD', 'DATE-OBS', \
           'MJD-OBS', 'CREATOR']

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
# Vizier query for galaxy catalog, options:
# 1. GLADE, VII/281
# 2. GWGC, VII/267
catalog = 1

# size for querying the galaxy catalog
# 100 for 100 lines of catalog
# -1 for the full catalog
size = -1

# filter for galaxy absolute magnitude
# for GLADE, there're B and K
# for GWGC, there're only B
filtro = 'B'

# distance range (in Mpc) for galaxy selection
# if trigger (CBC GW, or?) distance available,
#     auto select galaxy range
#      otherwise, applied settings here.
#    Notice: GWGC/GLADE is imcomplete beyond within 40/100 Mpc
limdist = 0,1000

# Ra range (in deg) for galaxy selection
limrag = 0,360

# Dec range (in deg) for galaxy selection
limdecg = -90,90

# absolute magnitude range (in mag) for galaxy selection
limmag = -12,-20
