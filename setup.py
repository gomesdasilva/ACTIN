#!/usr/bin/env python

try: from setuptools import setup
except: from distutils.core import setup

setup(name = 'actin',
      version = '1.0',
      description = 'Activity Indices Calculator',
      url = 'http://github.com/gomesdasilva/actin',
      author = 'Joao Gomes da Silva',
      author_email = 'Joao.Silva@astro.up.pt',
      license = 'MIT',
      keywords = 'astronomy activity fits harps harps-n',
      packages = ['actin'],
      #data_files = {"actin": ["config_lines.txt.dist"]},
      entry_points = {
        "console_scripts": ['actin = actin.actin:main']
        },
      include_package_data = True,
      install_requires = ['appdirs']
      )
