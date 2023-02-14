<h1 align="center"> BiSpeD </h1>

The Binary Spectral Disentangling (**BiSpeD**) is a PYTHON code for astronomy applications. The code can find and extract the spectral features of secondary companion from a binary system observation.


 ## Introduction

The BiSpeD code allows to estimate the mass ratio among companions (q-value). The spectral features of primary component (S<sub>A</sub>) are removed form observed spectra (S<sub>obs</sub>), and for each possible q-values defined by the user, a residual spectrum associated with secondary companion is made (S<sub>B</sub>). Finally, the code compares each S<sub>B</sub> spectrum with a synthetic templates catalogue to find the best match for mass ratio q and secondary effective temperature T<sub>eff</sub>. Eventually, the BiSpeD code works through specific functions (tasks) with particular aims. In the next topics, the main functions of BiSpeD with mandatory and optional arguments are presented.

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
> As the RV for a secondary companion is unknown (single-line spectroscopic binary), we can define a grid of different possible mass ratios q. The script uses the radial velocity of the primary component to estimate the RV value of the secondary component. For each mass ratio q, the spectra disentangling is applied using the function `spbina` and we obtain a mean spectrum corresponding to secondary companion. It convolves the mean spectrum with a synthetic template spectra list provided by the user. We obtained a value of cross-correlation function (CCF) for each q and each effective temperature T<sub>eff</sub>. The CCF value corresponds to the best synthetic template and this could be analyzed in a 3D density diagram., with X-axis the mass ratio and Y-axis the effective temperatures. 
> Mandatory parameters:
> - `lis`: file list of observed spectra to process. The file list is a text file containing a list of spectra images for input. We employed the IRAF terminology to specify this kind of file as the form "at file" and means the file name must be preceded by the symbol @. The easiest way to generate the file list is using the LINUX command **ls**;
> - `tmp`: path to folder containing spectra templates;
> - `vgamma`: estimated value for systemic radial velocity in km/s.
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
> This function computes the cross-correlation of radial velocities among two spectra using the Fast Fourier Transform method (FFT) to perform the convolution. The output is the full discrete linear convolution of the inputs with a dispersion grid defined by spectral region parameter ´wreg´. This function requires the previously estimation of primary and secondary mean spectra. First, the primary RV is determined by cross-correlation among observed spectrum and template; then the mean primary spectra behaviours are removed from observed spectra and the secondary RV determination is carry out. The cross-correlation result is fitting in the peak with a Gaussian function.
> Mandatory parameters:
> - `lis`: file list of observed spectra to process.
> Optional parameters:
> - `spa`: string for primary component spectrum name;
> - `spb`: string for secondary component spectrum name;
> - `ta`: spectrum template to compare with primary component;
> - `tb`: spectrum template to compare with secondary component;
> - `wreg`: spectral regions for cross-correlation analysis;
> - `aconv`: absorb convergence factor;
> - `keyjd`: header keyword that contains Julian Date;
> - `fitcont`: continuum subtract the spectra prior to correlation;
> - `interac`: process cross-correlation interactively.

- **rvextract**
> Analyze the convergence of iterative RV assignment of a spectra file list (mandatory parameter `lis`). The optional parameters are `output` to write a file with RV values and `graphic` to show the convergence as a function of iteration number.

- **setrvs**
> Set radial velocities for each spectrum from data-set using cross-correlation with templates defined by the user. This function can be employed with SB1 and SB2 binary system: in case of double-line stars, the parameter `tb` must be assigned to a template according to secondary possible spectral type. The only mandatory parameter is spectra file list `lis`. 
> Optional parameters:
> - `ta`: spectrum template to compare with primary component;
> - `tb`: spectrum template to compare with secondary component;
> - `wreg`: spectral regions for cross-correlation analysis.
> - `keyjd`: header keyword that contains Julian Date;
> - `fitcont`: continuum subtract the spectra prior to correlation;
> - `interac`: process cross-correlation interactively.

- **spbina**
> Compute spectral disentangling. The only mandatory parameter is spectra file list `lis`.
> Optional parameters:
> - `spa`: string for primary component spectrum name;
> - `spb`: string for secondary component spectrum name;
> - `nit`: number of iterations;
> - `frat`: estimated flux ratio among components;
> - `reject`: reject pixels using a sigma clipping algorithm;
> - `q`: mass ratio among components (if the RV of secondary component is unknown);
> - `vgamma`: estimated value for systemic radial velocity in km/s;
> - `obspha`: calculate spectra for all phases;
> - `showtit`: enable user interface (internal definition).

- **splot**
> Plot and show spectrum in FITS format. The mandatory parameter is the file name `file`.
> Optional parameters:
> - `xmin`: lower wavelength limit of initial graphic;
> - `xmax`: upper wavelength limit of initial graphic;
> - `ymin`: lower flux limit of initial graphic;
> - `ymax`: upper flux limit of initial graphic;
> - `scale`: flux scale factor;
> - `markpix`: mark dispersion pixel values;  
> - `newfig`: open spectrum in a new window;
> - `color`: graphic color.

- **uniform**
> Explore a spectra data-set (if interactive mode is active, given by the optional parameter `interac=True`) and select the best spectra to process. Finally, this script scale each spectrum from selected data-set to a mean continuum factor to make them morphologically consistent. The mandatory parameter is `lis`: file list of observed spectra to process; and optional parameter is interactive mode `interac`. 

- **vgrid**
> When the object has a higher Doppler semi-amplitude but the orbital period could no be fitting and the systemic velocity could not be estimated, we employed a script called `vgrid` that apply the method for a systemic velocities grid.
>Mandatory parameters:
> - `lis`: file list of observed spectra to process;
> - `tmp`: path to folder containing spectra templates.
> Optional parameters:
> - `svmin`: lower systemic radial velocity for grid km/s;
> - `svmax`: upper systemic radial velocity for grid km/s;
> - `step`: radial velocity step for meshing grid km/s;
> - `qmin`: minimum mass ratio for cross-correlation grid analysis;
> - `qmax`: maximum mass ratio for cross-correlation grid analysis;    
> - `deltaq`: mass increments for cross-correlation grid analysis;
> - `wreg`: spectral regions for cross-correlation analysis.
> - `nproc`: number of CPU cores to use in computing (it depends on computer resources).

- **vexplore**
> Show and explore results for systemic radial velocities grid analysis (only for **vgrid** mode). The mandatory parameter is the folder output name `obj`. The user can check the different cross-correlation analysis results for the systemic velocities grid.
