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
    IndexToDeclRa,
    DeclRaToIndex,
    RadecToThetaphi,
    ThataphiToRadec,
    query_ebv,
    query_vizier,
    rotate_map,
    vertices,
    ipix_in_box,
    contour,
    mcmc,
    gwdist,
    trigger_validation,
    get_skymap,
    slew_angle,
    prob_obs_hpmap,
    prob_obs_galaxies,
    read_filelist,
    moon_phase   
)

from .autogcndef import process_gcn

from manschedule import man_search

from .configure import (
    config_init,
    config_list
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
    mollzoom,
    contourview,
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

from scheme import (
    radec2skycell,
    skycell2radec,
    divide_OB,
    pointings,
    galaxies,
    outcat
)
    
from .server import lvc_server
