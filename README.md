# ACTIN
### Activity Indices Calculator

Reads fits files from HARPS and HARPS-N spectrographs and outputs user defined spectral indices.


### Quick start:

Usage:

`python actin.py files [-i indices] [-s output path] [-p output path] [-obj object_name] [-del True/False]`

where

`files` are either fits files of the formats e2ds, s1d, s1d_*_rv, or ADP, or
.rdb tables with required headers 'obj', 'date', 'bjd', 'wave', 'flux', 'error_pixel' (optional). Works for one or multiple files for one or multiple stars.

Options:

`-i` bool : Calculate indices chosen in lines.config file. Default is `True`.

`-s` bool : Save output to .rdb table in specified path. If 'False' no output table is saved (default).

`-p` bool : Save plots of the lines used to calculate the indices in the specified path. If 'same' uses the table output path in '-s'. If 'False' no plots are saved (default). if `True`.

`-obj` str : Object name to override the one from fits files in case the star has multiple names (ex. Proxima, ProximaCen, Gl551). Default is `None`.


*Any enquiries or bug reports to João Gomes da Silva, Joao.Silva(at)astro.up.pt*
