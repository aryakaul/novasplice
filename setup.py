#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.extension import Extension


setup(name='novasplice',
      version='0.0.1',
      description='Identify novel intronic splice sites',
      url='https://github.com/aryakaul/novasplice',
      packages=['novasplice'],
      install_requires=['maxentpy', 'pybedtools'],
      entry_points = {
        'console_scripts': ['novasplice=main:main'],
      }
)
