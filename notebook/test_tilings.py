from pst.pipeline.tilings import PstGetTilings
from pst.pipeline.triggers import PstParseTriggers

a=PstGetTilings()
b=PstParseTriggers()

a.generate(limra=[0,360],limdec=[-20,90],fovra=3.,fovdec=3.)
b.url('https://gracedb.ligo.org/api/superevents/S190510g/files/bayestar.fits.gz')

c=a.calc_sun()
print (c)
