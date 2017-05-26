#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:19:42 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

from setuptools import setup
import subprocess

try:
    subprocess.call(['conda','install','--file','requirements.txt'])
except Exception as e:
    pass


setup(name='gsit',
      description='utilities for the ionospheric remote sensing, ASI, GNSS, IMF, etc',
      author='Sebastijan Mrak',
      url='https://github.com/aldebaran1/gsit.git',
      install_requires=['numpy','GeoData', 'scipy'],
      packages=['gsit']

)