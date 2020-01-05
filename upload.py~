from __future__ import print_function
from builtins import input
from setuptools import find_packages, setup, Command
from shutil import rmtree
import io,sys,os

# Package meta-data.
NAME = 'pstools'
DESCRIPTION = 'pstools is a telescope pointing scheduler'
URL = 'https://sngyang.com/pstools'
URL1 = 'https://github.com/saberyoung/pst'
EMAIL = 'sheng.yang@inaf.it'
AUTHOR = 'Sheng Yang'
REQUIRES_PYTHON = '>=2.7'
VERSION = "0.0.0.1.0"

# What packages are required for this module to be executed
REQUIRED = [
    'meander','astropy','numpy','healpy','matplotlib','astroquery'
]

# What packages are optional
EXTRAS = {
    'multiple processing': ['multiprocessing'],
}

# define current directory
here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    with open(os.path.join(here, NAME, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds...')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution...')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine...')
#        os.system('twine upload dist/*')
        os.system('twine upload --skip-existing --verbose dist/*')

#        self.status('Pushing git tags...')
#        os.system('git add .')
#        os.system('git commit -m v{0}'.format(about['__version__']))
#        os.system('git tag v{0}'.format(about['__version__']))

        sys.exit()

setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,  
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    download_url = URL1,
    package_dir = {'': 'src'},
    packages=['pst'],    
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='MIT license',

    scripts=['bin/pstools','bin/slack'],
    classifiers=[                     
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],        
    cmdclass={
        'upload': UploadCommand,
    },
)
