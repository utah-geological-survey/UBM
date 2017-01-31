import sys
from setuptools import setup, find_packages
# To use:
#	   python setup.py bdist --format=wininst

from UBM import __version__, __name__, __author__

# trap someone trying to install flopy with something other
#  than python 2 or 3
if not sys.version_info[0] in [2,3]:
    print('Sorry, UBM not supported in your Python version')
    print('  Supported versions: 2 and maybe 3')
    print('  Your version of Python: {}'.format(sys.version_info[0]))
    sys.exit(1)  # return non-zero value for failure

long_description = 'A modified version of the Basin Characterization Model (BCM) to calculate recharge and runoff'

try:
    import pypandoc

    long_description = pypandoc.convert('README.md', 'rst')
except:
    pass

setup(name=__name__,
      description = 'download and calculate raster data for hydrologic model',
      long_description = long_description,
      version = __version__,
      author = __author__,
      author_email = 'paulinkenbrandt@utah.gov',
      url = 'https://github.com/inkenbrandt/UBM',
      license = 'LICENSE.txt',
      install_requires=["bs4",
                        "pymodis >= 2.0.4",
                        "wellapplication",
                        "pandas",
                        "numpy"],
      packages = find_packages(exclude=['contrib', 'docs', 'tests*']))