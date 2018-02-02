# ACTIN
### Activity Indices Calculator

Reads fits files from HARPS and HARPS-N spectrographs and outputs user defined spectral indices.

### Requires the following Python modules:
- numpy
- astropy


### Quick start:

Usage:

`python actin.py files [-i indices] [-s output path] [-p output path] [-obj object_name] [-tl target list] [-del True/False] [-w line flux weight]`

where

`files` are either fits files of the formats e2ds, s1d, s1d_*_rv, or ADP, or
.rdb tables with required headers `obj`, `date`, `bjd`, `wave`, `flux`, `error_pixel` (optional). Works for one or multiple files for one or multiple stars.

Options:

`-i` list : List of indices to calculate. Indices ids must match the ones in the config file `config_lines.txt`. If `False` no indices are calculated (default).

`-s` str : Save output to .rdb table in specified path. If `False` no output table is saved (default).

`-p` str : Save plots of the lines used to calculate the indices in the specified path. If `same` uses the table output path as specified in `-s`. If `False` no plots are saved (default).

`-tl` list : List of stars to select from `files`. Default is `None`.

`-obj` str : Object name to override the one from fits files in case the star has multiple names (ex. Proxima, ProximaCen, Gl551). Default is `None`.

`-w` str : Weight used to calculate the flux in the lines. `blaze` uses the blaze function as weight, flux is normalized by the sum of the weight function, `None` the flux is unweighted and normalized by number of pixels.


### Example:

`python actin.py ../fits/*/*e2ds_A.fits -i I_CaII I_Ha -s ../output -p same -del True -tl Gl273 Gl581`

This will execute ACTIN for all the subdirectories inside `../fits/` with files ending with `e2ds_A.fits`, calculate the indices `I_CaII` and `I_Ha`, output the data to `../output/star_names`, save spectra of the line regions to the same directory as data, and, before running the code, delete any output file that was previously there, in this case `Gl273_e2ds_actin.rdb` and `Gl581_e2ds_actin.rdb` files. Only fits files belonging to the stars chosen in `-tl` will be read, in this case `Gl273` and `Gl581`.

*Any enquiries or bug reports to João Gomes da Silva, Joao.Silva(at)astro.up.pt*
