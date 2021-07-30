from setuptools.extension import Extension
from setuptools import setup, find_packages
from setuptools.extension import Extension
import numpy

setup (
    name='plumber',
    version='0.20',
    packages=find_packages(),
    include_package_data=True,
    include_dirs=[numpy.get_include()],
    entry_points="""
    [console_scripts]
    plumber=plumber.plumber:main
    parang_finder=plumber.parang_finder:main
    """
)

