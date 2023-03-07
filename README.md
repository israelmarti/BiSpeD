<h1 align="center"> BiSpeD </h1>

The Binary Spectral Disentangling (**BiSpeD**) is a PYTHON code for astronomy applications. The code finds and extracts the spectral features of secondary companion star from a binary system observation.


 ## Introduction

The BiSpeD code allows to estimate the mass ratio among companions (*q*-value). The spectral features of primary component (*S<sub>A</sub>*) are removed from observed spectra (*S<sub>obs</sub>*), and for each possible *q*-values defined by the user, a residual spectrum associated with secondary companion is estimated (*S<sub>B</sub>*). Finally, the code compares each *S<sub>B</sub>* spectrum with a synthetic templates catalogue to find the best match for mass ratio *q* and secondary effective temperature *T<sub>eff</sub>*. Eventually, the BiSpeD code works through specific functions (tasks) with them aims. In *Main Functions* section, the main functions of BiSpeD with mandatory and optional arguments are presented.

## Installation

BiSpeD is distributed on PyPI as a universal wheel and is available on Linux/macOS and Windows and supports Python 3.5+.

```bash
pip install bisped
```

## Libraries

BiSpeD requires the the next libraries to run:
- Astropy_ >= 5.0.4
- Matplotlib_ >= 3.7.0
- Multiprocessing_ >= 2.6.2.1
- Numba_ >= 0.56.4
- Numpy_  >= 1.23.5
- Progress_ >= 1.5
- PyAstronomy_ >= 0.18.1
- SciPy_ >= 1.10.1
- Specutils_ >= 1.3.1


## Main Functions

- **find2c**

> This task can detect the best possible solution of mass ratio and effective temperature of secondary companion for a spectra dataset.
> As the radial velocity (*RV*) for a secondary companion is unknown (single-line spectroscopic binary), a grid of different possible mass ratios (*q*) must be defined. The script uses the *RV* of the primary component to estimate the *RV*-value of the secondary component. For each mass ratio *q*, the spectra disentangling is applied using the task `spbina` and a mean spectrum corresponding to secondary companion is obtanied. `find2c` convolves the mean spectrum with a synthetic template spectra list provided by the user. Finally, different values of cross-correlation function (CCF) for each *q* and each effective temperature *T<sub>eff</sub>* are obtained. Each CCF value corresponds to the best synthetic template and this could be analyzed in a 3D density diagram, where X-axis is the mass ratio and Y-axis the effective temperature. 
> Mandatory parameters:
> - `lis`: file list of observed spectra to process. The file list is a text file containing a list of spectra images for input (in FITS format). The IRAF terminology to specify this kind of file as the form "at file" is applied and it means the file name must be preceded by the symbol @. The easiest way to generate the file list is using the LINUX command **ls**;
> - `tmp`: full path to the folder containing spectra templates;
> - `vgamma`: estimated value for systemic radial velocity in km/s.
> Optional parameters:
> - `spa`: name of primary component mean spectrum;
> - `spb`: name of secondary component mean spectrum;
> - `qmin`: minimum mass ratio for cross-correlation grid;
> - `qmax`: maximum mass ratio for cross-correlation grid;    
> - `deltaq`: mass increments for cross-correlation grid;
> - `wreg`: spectral regions for cross-correlation analysis.
> - `nproc`: number of CPU cores to use in computing (it depends of computer resources).

- **hselect**
> Extract keyword values (`fields`) from a FITS image or file list (`img`). In case of file list, the file name must be preceded by the symbol @.

- **rvbina**
> This task computes the cross-correlation of radial velocities among two spectra using the Fast Fourier Transform method (FFT) to perform the convolution. The output is the full discrete linear convolution of the inputs showed in a dispersion grid defined by spectral region parameter ´wreg´. It requires the previously estimation of primary and secondary mean spectra. First, the primary RV is determined by cross-correlation among observed spectrum and template; then the mean primary spectra features are removed from observed spectra and the secondary RV determination is carry out. The cross-correlation result is fitted for the most significative peak with a Gaussian fitting.
> Mandatory parameters:
> - `lis`: file list of observed spectra to process. The file name must be preceded by the symbol @.
> Optional parameters:
> - `spa`: name of primary component mean spectrum;
> - `spb`: name of secondary component mean spectrum;
> - `ta`: spectrum template for comparison with primary component (in FITS format with or without extension);
> - `tb`: spectrum template for comparison with secondary component (in FITS format with or without extension);
> - `wreg`: spectral regions for cross-correlation analysis;
> - `aconv`: absorb convergence factor;
> - `keyjd`: header keyword that contains Julian Date;
> - `fitcont`: continuum subtraction of spectra dataset prior to correlation analysis;
> - `interac`: process cross-correlation interactively.

- **rvextract**
> Analyze the convergence of iterative **rvbina** task for a spectra file list (mandatory parameter `lis`). The optional parameters are `output` to write a file with RV values and `graphic` to show the RV convergence as a function of iteration number.

- **setrvs**
> Set radial velocities for each spectrum from data-set using cross-correlation with templates defined by the user. This task can be applied for SB1 and SB2 binary system; in case of double-line stars, the parameter `tb` must be assigned to a template according to possible secondary spectral type. The only mandatory parameter is spectra file list `lis`. 
> Optional parameters:
> - `ta`: spectrum template to compare with primary component (in FITS format with or without extension);
> - `tb`: spectrum template to compare with secondary component (in FITS format with or without extension);
> - `wreg`: spectral regions for cross-correlation analysis.
> - `keyjd`: header keyword that contains Julian Date;
> - `fitcont`: continuum subtraction of spectra dataset prior to correlation analysis;
> - `interac`: process cross-correlation interactively.

- **spbina**
> Compute spectral disentangling. The only mandatory parameter is spectra file list `lis`.
> Optional parameters:
> - `spa`: name of primary component mean spectrum;
> - `spb`: name of secondary component mean spectrum;
> - `nit`: number of iterations;
> - `frat`: estimated flux ratio among components;
> - `reject`: reject pixels using a sigma clipping algorithm;
> - `q`: mass ratio among components (if the RV of secondary component is unknown);
> - `vgamma`: estimated value for systemic radial velocity in km/s;
> - `obspha`: calculate spectra for all phases;
> - `showtit`: enable user interface (internal definition).

- **splot**
> Plot and show spectrum (must be in FITS format). The mandatory parameter is the file name `file`.
> Optional parameters:
> - `xmin`: lower wavelength limit for graphic;
> - `xmax`: upper wavelength limit for graphic;
> - `ymin`: lower flux limit for graphic;
> - `ymax`: upper flux limit for graphic;
> - `scale`: flux scale factor;
> - `markpix`: mark dispersion pixel values;  
> - `newfig`: open spectrum in a new window;
> - `color`: spectrum graphic color (see **matplotlib.pyplot** for colors availables).

- **uniform**
> Explore a spectra data-set (if interactive mode is active, given by the optional parameter `interac=True`) and select the bests spectra to process. Finally, this task scales each spectrum from selected data-set to a mean continuum factor to make them morphologically consistent. The mandatory parameter is `lis` for the file list of observed spectra to process; and optional parameter is interactive mode `interac`. 

- **vgrid**
> When the orbital period for primary RV values could not be fitting and the systemic velocity is unknown, the task `vgrid` can be applied to estimate the best systemic radial velocity of binary system through a systemic velocities grid around the most probable value. 
>Mandatory parameters:
> - `lis`: file list of observed spectra to process;
> - `tmp`: full path to the folder containing spectra templates;
> Optional parameters:
> - `svmin`: lower systemic radial velocity for grid in km/s;
> - `svmax`: upper systemic radial velocity for grid in km/s;
> - `step`: radial velocity step for meshing grid in km/s;
> - `qmin`: minimum mass ratio for cross-correlation grid;
> - `qmax`: maximum mass ratio for cross-correlation grid;     
> - `deltaq`: mass increments for cross-correlation grid;
> - `wreg`: spectral regions for cross-correlation analysis.
> - `nproc`: number of CPU cores to use in computing (it depends on computer resources).

- **vexplore**
> Show and explore results for systemic radial velocities grid analysis (only for **vgrid** task). The mandatory parameter `obj` is the folder output name obtained from **vgrid** task. The different cross-correlation analysis results for the systemic velocities grid can be checked with a interactive 3-D graphic.



.. _Numpy: https://www.numpy.org
.. _Astropy: https://www.astropy.org
.. _Matplotlib: https://matplotlib.org
.. _PyAstronomy: https://pyastronomy.readthedocs.io
.. _Specutils: https://specutils.readthedocs.io
.. _SciPy: https://scipy.org
.. _Numba: https://numba.pydata.org
.. _Multiprocessing: https://pypi.org/project/multiprocessing
.. _Progress: https://pypi.org/project/progress
