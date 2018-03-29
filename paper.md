---
title: 'ACTIN: A tool to calculate activity indices'
tags:
- stellar activity
- activity indices
- radial velocity
- exoplanets
- extrasolar planets
authors:
- name: João Gomes da Silva
  #orcid: 0000-0000-0000-1234
  affiliation: 1 # (Multiple affiliations must be quoted)
- name: Pedro Figueira
  #orcid: 0000-0000-0000-1234
  affiliation: 2
- name: Nuno Santos
  affiliation: "1, 3"
- name: João Faria
  affiliation: 1
affiliations:
- name: Instituto de Astrofísica e Ciências do Espaço, Porto, Portugal
  index: 1
- name: ESO, Santiago, Chile
  index: 2
- name: Departamento de Física e Astronomia, Faculdade de Ciências, Universidade do Porto, Portugal
  index: 3
date: 28 March 2018
bibliography: paper.bib
---

# Summary

ACTIN is a Python written programme which astronomers can use to easily calculate user defined activity indices. These indices can be used to compare variations in the stellar atmospheres with radial velocity (RV) signals in order to infer if the RV signal is stellar in origin or comes from an exterior perturber (e.g. another star or planet) or to study stellar activity *per se*.

This programme reads data from fits files returned from the pipelines of spectrographs, and also .rdb tables. It extracts useful spectral data and calculates spectral activity indices. The output is an .rdb table, with data about activity indices, RV, CCF profile parameters, and other useful data such as Julian date, CCF noise and errors. It also outputs timeseries plots of the activity indices and plots of the spectral lines used to compute the indices.

The activity indices are computed based on the method described in [@gomesdasilva2011] but updated with options to chose different weights and normalisations. A configuration file is provided with parameters (such as line core, bandpass, bandpass function, etc) for four widely used activity indices. These parameters can be changed and new spectral lines and indices can be added.

ACTIN was recently used in research and a paper was submitted to a peer-review astrophysics journal (Delgado-Mena et al.).

ACTIN 1.0 recognises fits files from the HARPS and HARPS-N spectrographs.

The code is available and will be updated on GitHub (https://github.com/gomesdasilva/ACTIN) and PyPI (https://pypi.org/project/actin/) and can be easily installed using pip.


# References
