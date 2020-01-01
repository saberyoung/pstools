from setuptools import setup
from pst. __version__ import version

setup(
    name='pstools',
    version=version,
    packages=[
        'pst',
        'pst.pipeline',
        'pst.view',
        'pst.circulate',
        'pst.interface',
    ],
    zip_safe = False
)
