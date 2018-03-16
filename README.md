# ACTIN
### Activity Indices Calculator

Reads fits files from HARPS and HARPS-N spectrographs, and rdb tables, and outputs user defined spectral indices.


### Requires the following Python modules:
- numpy
- astropy
- appdirs (included)


### Installation:

Copy the github repository to a directory of your choice and install via `python setup.py install`.

The `config_lines.txt` file is the line configuration file (instructions inside). This file is used to add line parameters to calculate any index as long as the line cores and bandpasses are inside the spectral range and spectral orders range (e2ds files) of the spectrograph. ACTIN will check this at start and give an error message if line parameters don't match the spectra.

This file is available from the directory each OS uses for storing user data (ACTIN must be run once, before the file becomes available*):

For OSX: `~/Library/Application Support/<AppName>`

For Windows: `C:\Documents and Settings\<User>\Application Data\Local Settings\<AppAuthor>\<AppName>` or possibly `C:\Documents and Settings\<User>\Application Data\<AppAuthor>\<AppName>`

For Linux: `~/.local/share/<AppName>`

*Use `actin -h` as the first run.

### Quick start:

Usage:

`actin files [-i indices] [-s output path] [-p output path] [-obj object_name] [-tl target list] [-del True/False] [-w line flux weight]`

where

`files` are either fits files of the formats e2ds, s1d, s1d_*_rv, or ADP, or
.rdb tables with required headers `obj`, `date`, `bjd`, `wave`, `flux`, `error_pixel` (optional). Works for one or multiple files for one or multiple stars.

Options:

`-i` list : List of indices to calculate. Indices ids must match the ones in the config file `config_lines.txt`. If `False` no indices are calculated (default).

`-s` str : Save output to .rdb table in specified path. If `False` no output table is saved (default).

`-p` str : Save plots of the lines used to calculate the indices in the specified path. If `same` uses the table output path as specified in `-s`. If `False` no plots are saved (default).

`-tl` list : List of stars to select from `files`. Default is `None`.

`-obj` str : Object name to override the one from fits files in case the star has multiple names (ex. Proxima, ProximaCen, Gl551). Default is `None`.

`-w` str : Function to weight the integrated flux. If `blaze` the flux is multiplied by the blaze function, if `None` the flux is not weighted (default).

`-n` str : Normalisation of the flux: if `band` the sum is normalised by the bandpass wavelength value in angstroms, if `npixels` by the number of pixels in the bandpass (default), if `weight` by the sum of the weight function inside the bandpass, if `None` the integrated flux is not normalised.


### Example:

`actin ../fits/*/*e2ds_A.fits -i I_CaII I_Ha -s ../output -p same -del True -tl Gl273 Gl581`

This will execute ACTIN for all the subdirectories inside `../fits/` with files ending with `e2ds_A.fits`, calculate the indices `I_CaII` and `I_Ha`, output the data to `../output/star_names`, save spectra of the line regions to the same directory as data, and, before running the code, delete any output file that was previously there, in this case `Gl273_e2ds_actin.rdb` and `Gl581_e2ds_actin.rdb` files. Only fits files belonging to the stars chosen in `-tl` will be read, in this case `Gl273` and `Gl581`. In this case, the flux is not weighted (-w None, default) and is normalised by the number of pixels in the passband (-n npixels, default).

*Any enquiries or bug reports to João Gomes da Silva, Joao.Silva(at)astro.up.pt*
