# PSTools: Pointing Scheduler in responce to Triggers

**pip install**

pip install pstools --user

**source install**

source pst.bash

**Building from source**

The development version can be obtained and installed from github:

    $ git clone https://github.com/saberyoung/pst.git
    $ cd pst
    $ (sudo) python setup.py install 
    $ (or) python setup.py install --prefix=/--user

**Functions**

This tool is developed for creating pointing list for telescopes.

If you do monitor search however the unvertainty region of trigger is relatively large (compare to the FoV of your telescope), the tool implements a priorization algorithm, by considering the trigger prob, luminous mass distribution, Milky way extinction, detection effciency, galaxy distance distributio, etc, to provide you a combined map, indicating you in which area is worth to map first.
Then as specified by you (in the default file), the tool create a pointing list in order to cover more prioization, meanwhile, considering also the visibility. Considering we're more interested in galaxy, vesides the tiling stategy, the tool provide also a galaxy monitoring strategy.
An automatic gcn listener is also available so that such tool can be runned autamatically in background, keeping track to GRB or GW. More functions are still under construction.

**Usage**

pstool(.py) would initialize a default parameter list, specify your demandings and excute pstool again.
see tutorial in sngyang.com/pstools
![alt tag](GW.gif)

**Find more**
https://sngyang.com/pst

Don't hesitate to contact me (sheng.yang@inaf.it) if you found any problems in using the tool. 'Anything that can possibly go wrong, will go wrong.'

--MURPHY'S LAW on debugging
