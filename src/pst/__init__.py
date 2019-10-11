# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
PStools is a package intended to contain core functionality and some
common tools needed for performing pointings scheduler for telescopes with
Python.
"""
from .version import __version__
from . import pstdef,default

from .autogcndef import process_gcn
from .manschedule import man_search
from .server import lvc_server
from .pstdef import (
    main,
    choose,    
    gwdist,
    get_skymap,  
    build_hp_map, 
    get_hp_map,
    read_filelist,
    decomposit)
from .query import (
    query_ebv,
    query_vizier)
from .visibility import (
    moon_phase,
    slew_angle,
    sunset,
    innight)
from .configure import (
    config_init,
    load_config)
from .link import (
    sendemail,
    wechat,
    createSSHClient,
    slack,
    phone)
from .priorization import (
    make_hpfitsmap,
    dist_galaxydist,
    calprob_gal,
    calprob_tile)
from .pstplot import (
    mollview,    
    pointview,
    distview,
    lumsview,
    verticeview,
    routeview,
    plot_coord,
    cumshow,
    plot_lines,
    cumshow,
    interactive_show,
    plot_sky)
from .scheme import (
    IndexToDeclRa,
    DeclRaToIndex,
    RadecToThetaphi,
    ThataphiToRadec,    
    vertices,
    ipix_in_box,
    compute_contours,
    compute_contours_1,
    radec2skycell,
    skycell2radec,
    divide_OB,
    pointings,
    remove_fields,
    gen_pointings,
    galaxies,
    pointngsshift,
    overlapregion)
