from setuptools import setup
setup(
    name='pstools',
    version='0.5.0',    
    packages=[
        'pst.view',
        'pst.utils',
    ],
    scripts=['bin/pstools',],
    entry_points = {
        'ampel.channels' : [
            'okc = ampel.contrib.okc.channels:load_channels',
        ],        
        'ampel.pipeline.t0.units' : [
            'OkcTransientFilter = ampel.contrib.okc.t0.OkcTransientFilter:OkcTransientFilter',
        ],
      }
)
