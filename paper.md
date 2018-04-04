---
title: 'ACTIN: A tool to calculate stellar activity indices'
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
  affiliation: "2, 1"
- name: Nuno C. Santos
  affiliation: "1, 3"
- name: João P. Faria
  affiliation: "1, 3"
affiliations:
- name: Instituto de Astrofísica e Ciências do Espaço, Universidade do Porto, CAUP, Rua das Estrelas, 4150-762 Porto, Portugal
  index: 1
- name: European Southern Observatory, Alonso de Cordova 3107, Vitacura, Santiago, Chile
  index: 2
- name: Departamento de Física e Astronomia, Faculdade de Ciências, Universidade do Porto, Rua do Campo Alegre, 4169-007 Porto, Portugal
  index: 3
#date: 28 March 2018
bibliography: paper.bib
---


# Summary

Magnetic activity in the atmospheres of stars produces a number of spectroscopic signatures that are visible in the shape and strength of spectral lines. These signatures can be used to access, among other things, the variability of the magnetic activity, or its influence on other parameters such as the measured radial velocity (RV). This latter is of utmost importance for the detection and characterization of planets orbiting other stars.

ACTIN is a Python program to calculate stellar activity indices. These indices can be used to study stellar activity *per se*, or to compare variations in the stellar atmospheres with RV signals. In turn, these can be used to infer if a given observed RV signal is of stellar origin or it stems from a true barycentric movement of the star. The usage of the program does not require python expertize and the indices can be customized through a configuration ASCII file.

The program reads input data either from .fits files returned by the pipelines of spectrographs, or from .rdb[^1] tables. It extracts automatically the spectral data required to calculate spectral activity indices. The output is an .rdb table, with the calculated stellar activity indices for each date (Julian Date), as well as the RV and Cross-Correlation Function profile parameters, if available. It also outputs timeseries plots of the activity indices and plots the spectral lines used to compute the indices.

The activity indices are computed based on the method described in [@gomesdasilva2011], updated with options to chose different weights and normalisations. A configuration file is provided with parameters (such as line core, bandpass, bandpass function, etc) for four widely used activity indices. These parameters can be changed and new spectral lines and indices can be added.

ACTIN was already used in research and a paper was submitted to a peer-review astrophysics journal [@delgadomena2018, submitted].

ACTIN 1.0 comes pre-formated for the .fits files from the HARPS and HARPS-N spectrographs, and can be easily adapted to other spectrographs such as ESPRESSO. Alternatively, since .rdb files can be read as input, ACTIN 1.0 can be easily used with any spectrum.

The code is available and will be updated on GitHub[^2] and PyPI[^3] and can be easily installed using pip.

[^1]:
rdb tables are tab separated ASCII files with headers separated from the data by minus ('-') symbols with the same length as the headers.
[^2]:
https://github.com/gomesdasilva/ACTIN.
[^3]:
https://pypi.org/project/actin/.

# References
