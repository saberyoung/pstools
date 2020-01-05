from setuptools import setup
from pst. __version__ import version

setup(
    name='pstools',
    version=version,
    python_requires='>=2.7',
    packages=[
        'pst',
        'pst.pipeline',
        'pst.view',
        'pst.circulate',
        'pst.interface',
    ],
    install_requires=[
        'meander',
        'astropy',
        'numpy',
        'healpy',
        'matplotlib',
        'astroquery',
    ],
    classifiers=[                     
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',       
    ], 
    zip_safe = False
)
