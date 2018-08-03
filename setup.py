#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

version_py = "novasplice/_version.py"
exec(open(version_py).read())

setup(name='novasplice',
      version=__version__,
      description='Identify novel intronic splice sites',
      url='https://github.com/aryakaul/novasplice',
      packages=['novasplice'],
      py_modules=["novasplice.novasplice"],
      install_requires=['maxentpy', 'pybedtools'],
      entry_points = {
        'console_scripts': ['novasplice=novasplice.novasplice:main'],
      }
)
