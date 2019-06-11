#!/usr/bin/env python

import os
import subprocess

try: from setuptools import setup
except: from distutils.core import setup


path = os.path.dirname(os.path.realpath(__file__))
version_file = os.path.join(path, "actin", "VERSION")

try:
    with open(version_file, 'r') as file:
        version = file.read()
except: version = "unknown"


setup(name = 'actin',
      version = version,
      description = 'Activity Indices Calculator',
      url = 'http://github.com/gomesdasilva/ACTIN',
      download_url = 'https://github.com/gomesdasilva/ACTIN/archive/v1.3.2.tar.gz',
      author = 'Joao Gomes da Silva',
      author_email = 'Joao.Silva@astro.up.pt',
      license = 'MIT',
      keywords = ['astronomy', 'activity', 'fits', 'harps', 'harps-n', 'espresso','radial velocity', 'exoplanets'],
      packages = ['actin'],
      entry_points = {
        "console_scripts": ['actin = actin.actin:main']
        },
      include_package_data = True,
      install_requires = ['appdirs']
      )


# This runs ACTIN from terminal, creates ACTIN directory with config file in Application Support and gives its location
#try: subprocess.call(["actin"])
#except: pass
