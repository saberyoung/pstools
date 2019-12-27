from setuptools import setup
from pst. __version__ import version

setup(
    name='pstools',
    version=version,
    packages=[
        'pst.pipeline',
    ],
    zip_safe = False
)
