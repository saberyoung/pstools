# for LVC GW only
# !!! for GRB, TBD
# from XML:
strxml = ['GraceID', 'AlertType', 'Group', 'FAR', \
           'Terrestrial', 'HasNS', 'Instruments', 'EventPage' ]
# from fits header:
strhdr = ['DISTMEAN', 'DISTSTD', 'DATE-OBS', \
           'MJD-OBS', 'CREATOR']

"""default settings"""

# default time last
deft = 24

# default nside of priorization map
# will change according to trigger healpix map if any
nside = 1024 

# for computer contour lines
# decrease nside in order to be faster
nside1 = 64

# size for querying the galaxy catalog
# 100 for 100 lines of catalog
# -1 for the full catalog
size = -1

# filter for galaxy absolute magnitude
# for GLADE, there're B and K
# for GWGC, there're only B
filtro = 'B'

# default figure size
figsize = '15,10'
                
# define colors for above contours
# number should be >= $contours
color_contour = ['b','g','k','y','c','m']
                
# define colors for tiling/galaxies
color_field = ['b','g','k','y','c','m']

# default paramater for healpix view
# healpix ordering
ordering = False
