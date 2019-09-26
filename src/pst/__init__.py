# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
PStools is a package intended to contain core functionality and some
common tools needed for performing pointings scheduler for telescopes with
Python.
"""
from .version import __version__
from . import pstdef

from .pstdef import (
    main,
    choose,    
    gwdist,
    get_skymap,  
    build_hp_map, 
    get_hp_map,
    read_filelist
)

from .query import (
    query_ebv,
    query_vizier
)

from .visibility import (        
    slew_angle,
    prob_obs_hpmap,
    prob_obs_galaxies,
    moon_phase
)

#from .mcmc import (
#    main
#)

from .autogcndef import process_gcn

from .manschedule import man_search

from .configure import (
    config_init,
    load_config
)

from .link import (
    sendemail_1,
    sendemail_2,
    wechat,
    createSSHClient,
    slack,
    phone
)

from .priorization import (
    make_hpfitsmap,
    dist_galaxydist,
    calprob_gal,
    calprob_tile
)

from .pstplot import (
    mollview,    
    pointview,
    distview,
    lumsview,
    dist_gauss,
    verticeview,
    plot_coord,
    cumshow,
    plot_lines,
    cumshow1,
    vis_plot,
    plot_all,
    interactive_show,
    show_scheduler,
    plot_sky
)

from .scheduler import (VST,schmidt)

from .scheme import (
    IndexToDeclRa,
    DeclRaToIndex,
    RadecToThetaphi,
    ThataphiToRadec,    
    vertices,
    ipix_in_box,
    compute_contours,
    radec2skycell,
    skycell2radec,
    divide_OB,
    pointings,
    galaxies,
    outcat
)
    
from .server import lvc_server
