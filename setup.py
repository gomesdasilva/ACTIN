#!/usr/bin/env python

import os
import subprocess

try: from setuptools import setup
except: from distutils.core import setup

setup(name = 'actin',
      version = '1.0',
      description = 'Activity Indices Calculator',
      url = 'http://github.com/gomesdasilva/actin',
      download_url = 'https://github.com/gomesdasilva/actin/archive/1.0.tar.gz',
      author = 'Joao Gomes da Silva',
      author_email = 'Joao.Silva@astro.up.pt',
      license = 'MIT',
      keywords = ['astronomy', 'activity', 'fits', 'harps', 'harps-n'],
      packages = ['actin'],
      entry_points = {
        "console_scripts": ['actin = actin.actin:main']
        },
      include_package_data = True,
      install_requires = ['appdirs']
      )

# Need to call ACTIN after installing to create user config file
subprocess.call(["actin", "-h"], stdout=open(os.devnull, 'wb'))
