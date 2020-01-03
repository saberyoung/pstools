from pst.pipeline.tilings import PstGetTilings
from pst.pipeline.triggers import PstParseTriggers
from pst.pipeline.galaxies import PstGetGalaxies
from pst.view.show import PstPlotter
a=PstGetTilings()
b=PstParseTriggers()
c=PstGetGalaxies()
d=PstPlotter()

#a.generate(limra=[0,360],limdec=[-20,90],fovra=3.,fovdec=3.)
b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')
#c.read()

import pylab as pl
pl.ion()
d.locshow(b,vmax=10**-5, ptype='m')
input('..')
d.contshow(b,nside=16)
input('..')
