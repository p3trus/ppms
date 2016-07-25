#  -*- coding: utf-8 -*-
#
# Slave, (c) 2012-2014, see AUTHORS.  Licensed under the GNU GPL.

import sys

from setuptools import setup, find_packages

requires = ['astropy']
try:
    import numpy
except ImportError:
    requires.append('numpy')

desc = 'Astropy Quantum Design PPMS file reader.'

setup(
    name='ppms',
    version='0.1',
    author='Marco Halder',
    author_email='marco.halder@frm2.tum.de',
    license = 'BSD',
    url='https://github.com/p3trus/ppms',
    description=desc,
    long_description=open('README.md').read(),
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
)
