{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This notebook demonstrates how to run a T2 analysis unit (in this case sncosmo) on a list (tar) of transient views. Transient views are summaries of all information regarding a transient exported from the core Ampel DB. Using these thus do not require an active local database. "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import re, os\n",
    "import ampel.utils.json_serialization as ampel_serialization\n",
    "import logging\n",
    "import gzip\n",
    "import importlib\n",
    "from ampel.ztf.pipeline.common.ZTFUtils import ZTFUtils\n",
    "import ampel.contrib.hu.t2.T2SNCosmo as T2SNCosmo\n",
    "import sncosmo"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Paths etc (edit as needed)\n",
    "# Transient view collection\n",
    "tvdir = '/home/jnordin/data/transientViews/'     # *** Change path *** # \n",
    "fname = tvdir+'tnscomplete_aug28.json.gz'        # *** Change path *** #\n",
    "outdir = '/home/jnordin/tmp/sncosmo/'\n",
    "# Logging    # \n",
    "logpath = '/home/jnordin/tmp/' # *** Change path *** # \n",
    "# Which model to fit\n",
    "fitmodel = 'salt2'\n",
    "#fitmodel = \"kilo_mej0.02_lant15_theta90\"\n",
    "# Limit fit to subset of sne\n",
    "fit_sne = []"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "... prob just the usual missing eof marker. Check that the number agrees\n",
      "Found 95 transient views\n"
     ]
    }
   ],
   "source": [
    "# Get a TransientView iterator\n",
    "if re.search('gz$',fname):\n",
    "    tv_iterator = ampel_serialization.load(gzip.open(fname,'rb'))\n",
    "else:    \n",
    "    tv_iterator = ampel_serialization.load(open(fname))\n",
    "tvlist = []\n",
    "try:\n",
    "    for tv in tv_iterator:\n",
    "        tvlist.append(tv)\n",
    "except EOFError as e:\n",
    "    print('... prob just the usual missing eof marker. Check that the number agrees')\n",
    "print('Found %s transient views'%(len(tvlist)))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Create a logger\n",
    "logger = logging.getLogger()\n",
    "logger.setLevel((logging.INFO))\n",
    "handler = logging.FileHandler(os.path.join(logpath, 'run_t2sncosmo.log'))\n",
    "logger.addHandler(handler)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    " What lightcurve fit is being carried out is detrmined by the run (and base) configurations which are supplied to the analysis unit."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [],
   "source": [
    "# The base config are channel parameters provided to the analysis unit at init. Usually not needed.\n",
    "base_config = {}"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "# The run config are parameters used by the analysis unit for each fit. For sncosmo these are usually model dependent\n",
    "if fitmodel=='salt2':\n",
    "    run_config = {\n",
    "            \"model\": \"salt2\",\n",
    "            \"with_upper_limits\": False,\n",
    "            \"sncosmo_kwarg\": {},\n",
    "            \"min_det\" : 8,\n",
    "            \"jd_reject_sigma\" : 5,\n",
    "            \"min_duration\" : 10,\n",
    "            \"max_duration\" : 90\n",
    "            }\n",
    "if fitmodel==\"kilo_mej0.02_lant15_theta90\":\n",
    "    run_config = {\n",
    "            \"model\": \"kilo_mej0.02_lant15_theta90\",\n",
    "            \"with_upper_limits\": False,\n",
    "            \"sncosmo_kwarg\": {\"bounds\":{'z':[0,0.05]}},\n",
    "            \"min_det\" : 8,\n",
    "            \"jd_reject_sigma\" : 5,\n",
    "            \"min_duration\" : 3,\n",
    "            \"max_duration\" : 16\n",
    "        }"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<module 'ampel.contrib.hu.t2.T2SNCosmo' from '/home/jnordin/github/ampelstatic/Ampel-contrib-HU/ampel/contrib/hu/t2/T2SNCosmo.py'>"
      ]
     },
     "execution_count": 7,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# During development\n",
    "importlib.reload(sncosmo)\n",
    "importlib.reload(T2SNCosmo)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [],
   "source": [
    "# This is where an instance of the T2 analysis module is created\n",
    "myt2 = T2SNCosmo.T2SNCosmo(logger, base_config)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Trying to fit model salt2 to ZTF17aacmbhk\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458720.6996759004 step= 2.1399106479994954 bounds= (2458637.9891088, 2458744.9846412) \n",
      "x0 0.004843591753999345 step= 0.00048435917539993456 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "563 function calls; 11 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:matplotlib.font_manager:Could not open font file /usr/share/fonts/truetype/noto/NotoColorEmoji.ttf: In FT2Font: Could not set the fontsize\n",
      "INFO:matplotlib.font_manager:generated new fontManager\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 0.16920902253544595, 't0': 2458707.665363455, 'x0': 0.0008961047180545107, 'x1': 4.999999975712109, 'c': 0.0561828496705159, 'mwebv': 0.2643604450621527, 'mwr_v': 3.1, 'z.err': 0.012667375606292969, 't0.err': 0.4818084316793829, 'x0.err': 2.9968192461634183e-05, 'x1.err': 0.0910888680801083, 'c.err': 0.04583065189430635}\n",
      "Trying to fit model salt2 to ZTF18aabilqu\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18aacmnbd\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18aadfbdc\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18aboerzk\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18abuanlw\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18abvfejf\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458723.7396990997 step= 2.0807907400000842 bounds= (2458640.9451042, 2458744.9846412) \n",
      "x0 0.0016777609224299355 step= 0.00016777609224299356 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "554 function calls; 13 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 0.033172655838799206, 't0': 2458708.525339911, 'x0': 0.0003169117908950929, 'x1': 4.999999999428496, 'c': 0.3124518569940786, 'mwebv': 0.0779294094406921, 'mwr_v': 3.1, 'z.err': 0.012604265777461956, 't0.err': 0.5792783598881215, 'x0.err': 1.85915785196425e-05, 'x1.err': 0.18140602258037664, 'c.err': 0.0824844652598451}\n",
      "Trying to fit model salt2 to ZTF18abxhnyj\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18abzbprc\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acbvkwl\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acckrsc\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18achrnju\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acrgczy\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acrvahh\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acszild\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acvsbiq\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18acxhvzb\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF18adbntlu\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aaajarm\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aablivc\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aacqcrq\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458723.7609491 step= 2.1399106479994954 bounds= (2458637.9891088, 2458744.9846412) \n",
      "x0 0.001204318862493987 step= 0.00012043188624939871 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "84 function calls; 10 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too long duration\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 9.866370092381782e-05, 't0': 2458743.7827304644, 'x0': 0.00035057377865255386, 'x1': -4.5722539407704375, 'c': -0.3080754206884615, 'mwebv': 0.08928661273855562, 'mwr_v': 3.1, 'z.err': 0.0, 't0.err': 0.0, 'x0.err': 0.0, 'x1.err': 0.0, 'c.err': 0.0}\n",
      "Trying to fit model salt2 to ZTF19aadbmum\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aagxfbk\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aagzfoj\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aahtowl\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458699.7493403 step= 2.0208303239941596 bounds= (2458643.9493403, 2458744.9908565) \n",
      "x0 0.0010457282917532267 step= 0.00010457282917532268 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "888 function calls; 4 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 5.251263868188972e-11, 't0': 2458672.9917682754, 'x0': 0.0032989697215843597, 'x1': -3.308086419540012, 'c': 0.2028433858373937, 'mwebv': 0.1691759015653947, 'mwr_v': 3.1, 'z.err': 0.002747024862644848, 't0.err': 0.00017125997692346573, 'x0.err': 0.000972957721761122, 'x1.err': 0.3467291109711741, 'c.err': 0.1499471249136155}\n",
      "Trying to fit model salt2 to ZTF19aaipyyh\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aakpjhp\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458695.7634258997 step= 2.100442361999303 bounds= (2458639.9634259, 2458744.985544) \n",
      "x0 0.002562595125873436 step= 0.0002562595125873436 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "377 function calls; 8 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too long duration\n",
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 2.0779709841711735e-09, 't0': 2458682.3283666796, 'x0': 0.001145149802954562, 'x1': 4.999993553651622, 'c': 0.02617735891443118, 'mwebv': 0.04083459344084727, 'mwr_v': 3.1, 'z.err': 0.0016673756439765045, 't0.err': 0.5387137376237661, 'x0.err': 7.446264130413178e-05, 'x1.err': 0.6580875851074506, 'c.err': 0.0518623359348791}\n",
      "Trying to fit model salt2 to ZTF19aalgghb\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19aavntzk\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abegrfy\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458686.6935764 step= 2.5807805540040136 bounds= (2458615.9488079, 2458744.9878356) \n",
      "x0 0.003224232159601188 step= 0.0003224232159601188 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "301 function calls; 22 dof.\n",
      "   model fit output:  {'z': 2.1265768856970624e-08, 't0': 2458683.0780919758, 'x0': 0.0008983138509830173, 'x1': 4.999999713229091, 'c': 0.11391920229075492, 'mwebv': 0.3316237621998447, 'mwr_v': 3.1, 'z.err': 0.006558563545757118, 't0.err': 0.5450782566331327, 'x0.err': 5.955899216099664e-05, 'x1.err': 0.12849376353351527, 'c.err': 0.055702788073484166}\n",
      "Trying to fit model salt2 to ZTF19abfinrl\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458696.6482639 step= 2.5204083320032806 bounds= (2458618.967419, 2458744.9878356) \n",
      "x0 0.002078288391443158 step= 0.00020782883914431582 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "376 function calls; 30 dof.\n",
      "   model fit output:  {'z': 1.086855024645672e-06, 't0': 2458686.72796392, 'x0': 0.0005469348376601214, 'x1': 4.999991918695873, 'c': 0.058100639699118606, 'mwebv': 0.19329911578963044, 'mwr_v': 3.1, 'z.err': 0.021090037448093416, 't0.err': 0.5309396274387836, 'x0.err': 3.155748645259229e-05, 'x1.err': 0.3252924809093001, 'c.err': 0.05466891817762254}\n",
      "Trying to fit model salt2 to ZTF19abjgcxh\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458720.6978588 step= 2.1399287040065973 bounds= (2458637.9900231, 2458744.9864583) \n",
      "x0 0.005861971324022352 step= 0.0005861971324022352 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "550 function calls; 10 dof.\n",
      "   model fit output:  {'z': 0.19999999177602804, 't0': 2458708.4017453454, 'x0': 0.0010905481133507474, 'x1': 4.9999999950335745, 'c': -0.08980156918377569, 'mwebv': 0.4332769755765021, 'mwr_v': 3.1, 'z.err': 0.013629902352551765, 't0.err': 0.49461961374618113, 'x0.err': 3.499130353051027e-05, 'x1.err': 0.08726845767938674, 'c.err': 0.04687172026451425}\n",
      "Trying to fit model salt2 to ZTF19abjgdet\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458705.6548958 step= 2.139928702004254 bounds= (2458637.9895718, 2458744.9860069) \n",
      "x0 0.0016096424620532096 step= 0.00016096424620532098 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "375 function calls; 8 dof.\n",
      "   model fit output:  {'z': 0.17376436032026274, 't0': 2458707.0286395433, 'x0': 0.00026841293218923833, 'x1': 4.999999583363747, 'c': 0.3219742689055054, 'mwebv': 0.20555845737147707, 'mwr_v': 3.1, 'z.err': 0.04417923296420847, 't0.err': 1.1642205473035574, 'x0.err': 2.34491169424353e-05, 'x1.err': 1.3331209153108086, 'c.err': 0.12315952564222798}\n",
      "Trying to fit model salt2 to ZTF19abjgdhb\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458694.2879051 step= 2.1399106479994954 bounds= (2458637.9891088, 2458744.9846412) \n",
      "x0 0.0016906957789012463 step= 0.00016906957789012463 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "491 function calls; 9 dof.\n",
      "   model fit output:  {'z': 2.3297544893274846e-07, 't0': 2458694.3839377193, 'x0': 0.00023279149577002677, 'x1': -1.7785837062499543, 'c': 1.0600796560711947, 'mwebv': 0.06850308456725536, 'mwr_v': 3.1, 'z.err': 0.004149930148303904, 't0.err': 1.1989267743192613, 'x0.err': 2.8256613579893937e-05, 'x1.err': 0.4882344331905466, 'c.err': 0.08124044390600327}\n",
      "Trying to fit model salt2 to ZTF19abjgdko\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458695.6648148 step= 2.1399287040065973 bounds= (2458637.9900231, 2458744.9864583) \n",
      "x0 0.007898481939383556 step= 0.0007898481939383556 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "316 function calls; 15 dof.\n",
      "   model fit output:  {'z': 0.07514852974927676, 't0': 2458674.158989674, 'x0': 0.005367391323588258, 'x1': 0.7640894668109377, 'c': -0.0845954806327045, 'mwebv': 0.49657183611182965, 'mwr_v': 3.1, 'z.err': 0.04751646134118555, 't0.err': 3.3239125558175147, 'x0.err': 0.0009523097685674066, 'x1.err': 1.2438117276461438, 'c.err': 0.10158500833294692}\n",
      "Trying to fit model salt2 to ZTF19abjgdnl\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458694.3865509 step= 2.1399287040065973 bounds= (2458637.9900231, 2458744.9864583) \n",
      "x0 0.0014591448699056654 step= 0.00014591448699056655 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "340 function calls; 4 dof.\n",
      "   model fit output:  {'z': 0.07654584104851805, 't0': 2458677.7525695437, 'x0': 0.0010051346574691557, 'x1': 2.814951826292962, 'c': 0.11253740685555247, 'mwebv': 0.42014164580106605, 'mwr_v': 3.1, 'z.err': 0.04644907457700123, 't0.err': 2.52477769064717, 'x0.err': 0.0004597956035077392, 'x1.err': 2.1072004585131143, 'c.err': 0.2617134914921624}\n",
      "Trying to fit model salt2 to ZTF19abjpick\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458720.7974074 step= 2.500788193996996 bounds= (2458619.9474884, 2458744.9868981) \n",
      "x0 0.0012261263268835676 step= 0.00012261263268835678 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "325 function calls; 14 dof.\n",
      "   model fit output:  {'z': 4.28187093326482e-07, 't0': 2458706.4764599116, 'x0': 0.0002491023934083394, 'x1': 4.999999795757875, 'c': 0.3504518261136442, 'mwebv': 0.12689738645339468, 'mwr_v': 3.1, 'z.err': 0.17096706821029595, 't0.err': 0.7853314646054059, 'x0.err': 2.3072472776918172e-05, 'x1.err': 0.7400726423271227, 'c.err': 0.0796972680354846}\n",
      "Trying to fit model salt2 to ZTF19abjpifs\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458715.8078356 step= 2.5807618040032687 bounds= (2458615.9492593, 2458744.9873495) \n",
      "x0 0.0009013172007907775 step= 9.013172007907775e-05 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "298 function calls; 14 dof.\n",
      "   model fit output:  {'z': 0.15628688403484284, 't0': 2458706.8034261963, 'x0': 0.0001755172687470407, 'x1': 4.999999947849309, 'c': 0.048581991350990616, 'mwebv': 0.10751019456919733, 'mwr_v': 3.1, 'z.err': 0.03401567945476743, 't0.err': 1.0492361667566001, 'x0.err': 1.1709640426643782e-05, 'x1.err': 0.5290142711929411, 'c.err': 0.07736159652887875}\n",
      "Trying to fit model salt2 to ZTF19abkfxfb\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458695.7634258997 step= 2.100442361999303 bounds= (2458639.9634259, 2458744.985544) \n",
      "x0 0.004977617157845258 step= 0.0004977617157845258 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "331 function calls; 10 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 8.172891696744955e-10, 't0': 2458682.15597512, 'x0': 0.0025656375924318934, 'x1': 3.1550889932125177, 'c': 0.07571690548258747, 'mwebv': 0.17831912572878247, 'mwr_v': 3.1, 'z.err': 0.0016868531710103386, 't0.err': 0.5096778150182217, 'x0.err': 0.0002034775263529938, 'x1.err': 0.4756240002828003, 'c.err': 0.04802138860291083}\n",
      "Trying to fit model salt2 to ZTF19abkfzrq\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abkhdnw\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458695.7689699 step= 2.1004664359986784 bounds= (2458639.9689699, 2458744.9922917) \n",
      "x0 0.0008533395894819908 step= 8.533395894819909e-05 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "614 function calls; 13 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 3.873910803253722e-09, 't0': 2458696.797275948, 'x0': 0.00018627729595311287, 'x1': 4.999999997984581, 'c': 0.5860780778303334, 'mwebv': 0.08452119136207309, 'mwr_v': 3.1, 'z.err': 0.005603486613833447, 't0.err': 0.8965492844581604, 'x0.err': 2.562893622618217e-05, 'x1.err': 0.2691064514310151, 'c.err': 0.1079436929507045}\n",
      "Trying to fit model salt2 to ZTF19abkinmd\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19ablpdid\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458710.7957754997 step= 2.500788193996996 bounds= (2458619.9474884, 2458744.9868981) \n",
      "x0 0.0016367333390315175 step= 0.00016367333390315177 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "416 function calls; 6 dof.\n",
      "   model fit output:  {'z': 0.029326579060770866, 't0': 2458707.685326536, 'x0': 0.0003104048419551772, 'x1': 3.828021103264579, 'c': 0.29754620744307125, 'mwebv': 0.12401118762885932, 'mwr_v': 3.1, 'z.err': 0.0489128862969333, 't0.err': 1.6577994695398957, 'x0.err': 4.510072548153197e-05, 'x1.err': 1.8251919423543477, 'c.err': 0.09286383297576184}\n",
      "Trying to fit model salt2 to ZTF19ablpfhz\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458710.6916088 step= 2.120842128004879 bounds= (2458638.9457292, 2458744.9878356) \n",
      "x0 0.005236669989476733 step= 0.0005236669989476734 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "432 function calls; 15 dof.\n",
      "   model fit output:  {'z': 0.07245442975690829, 't0': 2458708.4985121107, 'x0': 0.0013232890806622823, 'x1': 1.1650354459950378, 'c': 0.030554036123411787, 'mwebv': 0.21871148033120985, 'mwr_v': 3.1, 'z.err': 0.06800156405803212, 't0.err': 0.6977473618462682, 'x0.err': 6.960868793045643e-05, 'x1.err': 1.6257524434131012, 'c.err': 0.14798107102024805}\n",
      "Trying to fit model salt2 to ZTF19abnjest\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458723.7396990997 step= 2.020409027999267 bounds= (2458643.963287, 2458744.9837384) \n",
      "x0 0.0012346665894448863 step= 0.00012346665894448863 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "903 function calls; 6 dof.\n",
      "   model fit output:  {'z': 0.0898918565998713, 't0': 2458714.0641002823, 'x0': 0.0002592933119695253, 'x1': 3.5187927523168607, 'c': 0.18578403094075968, 'mwebv': 0.1727194844939883, 'mwr_v': 3.1, 'z.err': 0.04841172562039733, 't0.err': 1.4954253735486418, 'x0.err': 5.21590889106987e-05, 'x1.err': 2.37421459885848, 'c.err': 0.18571881930786427}\n",
      "Trying to fit model salt2 to ZTF19abouhwr\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458705.7417360996 step= 2.139864584002644 bounds= (2458637.9909606, 2458744.9841898) \n",
      "x0 0.000787462895886285 step= 7.87462895886285e-05 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "609 function calls; 7 dof.\n",
      "   model fit output:  {'z': 0.026549805350718093, 't0': 2458707.643713875, 'x0': 0.00016248246874853007, 'x1': 4.999999184581553, 'c': 0.2746759122475648, 'mwebv': 0.09279731990013568, 'mwr_v': 3.1, 'z.err': 0.11106733438575905, 't0.err': 1.6628953374456614, 'x0.err': 2.1697399458909048e-05, 'x1.err': 0.4766011793345504, 'c.err': 0.15927658266168065}\n",
      "Trying to fit model salt2 to ZTF19aboumcn\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458715.8037153 step= 1.9006222220044584 bounds= (2458649.9548958, 2458744.9860069) \n",
      "x0 0.0014365391422838912 step= 0.00014365391422838913 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "332 function calls; 14 dof.\n",
      "   model fit output:  {'z': 0.14821971709118584, 't0': 2458712.046880441, 'x0': 0.0003004336510293887, 'x1': 0.4390492691531378, 'c': 0.006825739956343391, 'mwebv': 0.09754500869484248, 'mwr_v': 3.1, 'z.err': 0.046785383491462135, 't0.err': 0.7989184681791812, 'x0.err': 2.6209784245776243e-05, 'x1.err': 1.4543321104107352, 'c.err': 0.0665150683687894}\n",
      "Trying to fit model salt2 to ZTF19abpaspf\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458715.7041667 step= 1.8801696759928017 bounds= (2458650.9770602, 2458744.985544) \n",
      "x0 0.0022306959110814486 step= 0.00022306959110814486 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "393 function calls; 3 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 2.267454002335967e-11, 't0': 2458690.6685129427, 'x0': 0.0006127374458170521, 'x1': 4.999998718834304, 'c': 0.1922630923543216, 'mwebv': 0.19629288137865186, 'mwr_v': 3.1, 'z.err': 3.6078584122634094e-05, 't0.err': 1.2459737190511078, 'x0.err': 0.00011304858855467545, 'x1.err': 0.37987467812269404, 'c.err': 0.13563336524668967}\n",
      "Trying to fit model salt2 to ZTF19abplekc\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abplfxs\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458723.7401968 step= 1.8397662020009011 bounds= (2458652.9976968, 2458744.9860069) \n",
      "x0 0.006957622820539432 step= 0.0006957622820539432 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "497 function calls; 4 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 0.043387756436194316, 't0': 2458722.5195794175, 'x0': 0.0016337515145041061, 'x1': 0.16834343843011723, 'c': -0.009596238483109643, 'mwebv': 0.09698286387040295, 'mwr_v': 3.1, 'z.err': 0.030503285240895874, 't0.err': 1.144694346236065, 'x0.err': 0.00011122653752125783, 'x1.err': 1.2440318055233153, 'c.err': 0.05388284383795616}\n",
      "Trying to fit model salt2 to ZTF19abpnryh\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abprrpm\n",
      "Initial parameters:\n",
      "z 0.1 step= 0.004 bounds= [0, 0.2] \n",
      "t0 2458715.8091897997 step= 1.7998224520031363 bounds= (2458654.9957755, 2458744.9868981) \n",
      "x0 0.001848736264298349 step= 0.0001848736264298349 \n",
      "x1 0.0 step= 0.2 bounds= [-5, 5] \n",
      "c 0.0 step= 0.14 bounds= [-2, 5] \n",
      "\n",
      "437 function calls; 5 dof.\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too short duration\n",
      "INFO:root:Exit T2SNcosmo, too short duration\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too short duration\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too short duration\n",
      "INFO:root:Exit T2SNcosmo, too short duration\n",
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   model fit output:  {'z': 0.15935138301212504, 't0': 2458723.6065141466, 'x0': 0.0003324517968714707, 'x1': 4.122032842196239, 'c': -0.05645695130766204, 'mwebv': 0.11082329211082145, 'mwr_v': 3.1, 'z.err': 0.062450139969240745, 't0.err': 2.734160903375596, 'x0.err': 4.210760667462541e-05, 'x1.err': 7.747696431396694, 'c.err': 0.09986391119215798}\n",
      "Trying to fit model salt2 to ZTF19abqgxla\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqgxpt\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqgxqb\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqgyxp\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqgyzq\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqgyzx\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqhbuk\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqhbxw\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqhbze\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqmjdq\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqmpsr\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqycui\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqykei\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqykst\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqykti\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqykuc\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqykyd\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqyogn\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqyohb\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqyonh\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqyoom\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqyopw\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqypak\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abqyphl\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abroguf\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abroiqa\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abrojrv\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abrorbx\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abtrxfn\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abtrzvf\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n",
      "INFO:root:Exit T2SNcosmo, too few det\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abtsluq\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abtsphy\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abtspvl\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abtsroz\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttapg\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttlwl\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttqkk\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttqpp\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttqpw\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttqxa\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttrze\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttscy\n",
      "   no fit done (insufficient data?)\n",
      "Trying to fit model salt2 to ZTF19abttsst\n",
      "   no fit done (insufficient data?)\n",
      "CPU times: user 29.3 s, sys: 4.55 s, total: 33.8 s\n",
      "Wall time: 26.3 s\n"
     ]
    }
   ],
   "source": [
    "%%time\n",
    "for tv in tvlist:\n",
    "    ztf_name = ZTFUtils.to_ztf_id(tv.tran_id)    \n",
    "    if len(fit_sne)>0 and not ztf_name in fit_sne : continue\n",
    "    print('Trying to fit model %s to %s'%(fitmodel,ztf_name))    \n",
    "    run_config['savefile'] = outdir+'%s.pdf'%(ztf_name)   \n",
    "    lc = tv.get_latest_lightcurve()\n",
    "    out = myt2.run(lc, run_config)\n",
    "    if len(out)==0:\n",
    "        print('   no fit done (insufficient data?)')\n",
    "    else:\n",
    "        print('   model fit output: ',out['fit_results'])\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
