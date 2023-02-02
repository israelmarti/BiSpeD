<h1 align="center"> BiSpeD </h1>

>The Binary Spectral Disentangling (**BiSpeD**) is a PYTHON code for astronomy spectroscopy application.


 ## Introduction

The BiSpeD code allows us to compute different mass ratios among companions (q-values) simultaneously, according to the computing resources. The spectral features of primary component (S<sub>A</sub>) are removed form observed spectra (S<sub>obs</sub>), and for each possible $q$-value defined by the user, a spectrum associated with secondary companion is made (S<sub>B</sub>). Finally, the code compares each S<sub>B</sub> spectrum with a synthetic templates catalogue to find the best match for mass ratio $q$ and secondary effective temperature T<sub>eff</sub>. Eventually, the BiSpeD code works through specific functions (tasks) with particular aims. In the next Subsections we present the main functions of BiSpeD available for the users with mandatory and optional arguments.

## Installation

BiSpeD is distributed on PyPI as a universal wheel and is available on Linux/macOS and Windows and supports Python 3.5+.

```bash
pip install bisped
```

## Libraries

BiSpeD requires the next libraries to run:
- Numpy (https://www.numpy.org)
- Astropy (https://www.astropy.org)
- Matplotlib (https://matplotlib.org)
- PyAstronomy (https://pyastronomy.readthedocs.io)
- Specutils (https://specutils.readthedocs.io)
- SciPy (https://scipy.org)
- Numba (https://numba.pydata.org)
- Multiprocessing (https://docs.python.org)


## Main Functions

- **find2c**

> This script can detect the best possible solution of mass ratio and effective temperature of secondary companion for a spectra dataset. 

>  As the RV for a secondary companion is unknown (single-line spectroscopic binary), we can define a grid of different possible mass ratios q. The script uses the radial velocity of the primary component to estimate the RV value of the secondary component. For each mass ratio q, the spectra disentangling is applied using the function `spbina` and we obtain a mean spectrum corresponding to secondary companion. It convolves the mean spectrum with a synthetic template spectra list provided by the user. We obtained a value of cross-correlation function (CCF) for each q and each effective temperature T<sub>eff</sub>. The CCF value corresponds to the best synthetic template and this could be analyzed in a 3D density diagram., with X-axis the mass ratio and Y-axis the effective temperatures. 
> Mandatory parameters:
> - `lis`: file list of observed spectra to process. The file list is a text file containing a list of spectra images for input. We employed the IRAF terminology to specify this kind of file as the form "at file" and means the file name must be preceded by the symbol @. The easiest way to generate the file list is using the LINUX command **ls**;
> - `tmp`: path to folder containing spectra templates;
> - `vgamma`: estimated value for systemic radial velocity in km/s.
> 
> Optional parameters:
> - `spa`: string for primary component spectrum name;
> - `spb`: string for secondary component spectrum name;
> - `qmin`: minimum mass ratio for cross-correlation grid analysis;
> - `qmax`: maximum mass ratio for cross-correlation grid analysis;    
> - `deltaq`: mass increments for cross-correlation grid analysis;
> - `wreg`: spectral regions for cross-correlation analysis.
> - `nproc`: number of CPU cores to use in computing (it depends on computer resources).


- **hselect**
> Extract keyword values (`fields`) from a FITS image or file list (`img`).

- **rvbina**
