Welcome to the pstools documentation
=====================================

`Pointing Scheduler Tools (PSTools) <https://sngyang.com/pstools>`_ is 
a Python package that can keep monitoring `GCN <https://gcn.gsfc.nasa.gov/>`_ 
system via `voevent <https://gcn.gsfc.nasa.gov/voevent.html>`_, 
if any alert received handle the pixelated trigger (GW/GRB/AMON...) skymap 
on the sphere (or offline input a `Healpix <https://healpix.jpl.nasa.gov/>`_ map), 
and/or weighted by local luminous mass (via Vizier), in order to generate a 
series of geodesic rectangles as for telescope pointings, which are sorted 
in order astronomers to cover the sky localization with less pointings/time.

`pstools` provides following utilities:

* generate pointing list for specific telescopes, depending on its location (i.e. altitude, latitude, longitude) and FoV (note that depending on the FoV of telescope and the trigger localization, there're two complementary approaches, tiling and galaxy search).
* sort the pointings based on: 1) probability (which was an convolution yields of trigger probability, 3D luminous mass distribution, etc, specified by users), and visibility/executability, e.g. galactic extinction, airmass, telescope slewing time, etc.
* schedule also for multiple telescopes with various follow-up strategies.
* hire a GCN listener in order to schedule telescopes automatically, and provide machine the possibilities to interact with telescopes via API (which would make the follow-up search immediate).
* alerting users at the same time by email/SMS/slack.
* implement a Monto Carlo approach to optimise the pointing list.
* enable users to interact with machine via slack, so that the event advocator avoid wasting time for basic archive checks, e.g. TNS, NED, galactic extinction, airmass plots, etc.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   usage
   ref
   example

Author
------
Sheng Yang
http://www.sngyang.com
