<h1 align="center"> BiSpeD </h1>

The Binary Spectral Disentangling (**BiSpeD**) is a PYTHON library for astronomy applications. The package includes tasks to manipulate and processing spectral observations of binary stars; the main goal of this development is to finds and extracts the spectral features of secondary companion star.


 ## Introduction

The BiSpeD code allows to estimate the mass ratio among companions (*q*-value). The spectral features of primary component (*S<sub>A</sub>*) are removed from observed spectra (*S<sub>obs</sub>*), and for each possible *q*-values defined by the user, a residual spectrum associated with secondary companion is estimated (*S<sub>B</sub>*). Finally, the code compares each *S<sub>B</sub>* spectrum with a synthetic templates catalogue to find the best match for mass ratio *q* and secondary effective temperature *T<sub>eff</sub>*. Eventually, the BiSpeD code works through specific functions (tasks) with them aims. In *Main Functions* section, the main functions of BiSpeD with mandatory and optional arguments are presented.

## Installation

BiSpeD is distributed on PyPI as a universal wheel and is available on Linux/macOS and Windows and supports Python 3.5+.

```bash
pip install bisped
```
For optimized interactive usefull run from [IPython](https://ipython.org/install.html).

## Dependencies

BiSpeD requires the the next dependencies to run:
- [Astropy](https://www.astropy.org) >= 5.0.4
- [Matplotlib](https://matplotlib.org) >= 3.7.0
- [Multiprocessing](https://pypi.org/project/multiprocessing) >= 2.6.2.1
- [Numba](https://numba.pydata.org) >= 0.56.4
- [Numpy](https://www.numpy.org)  >= 1.23.5
- [Progress](https://pypi.org/project/progress) >= 1.5
- [PyAstronomy](https://pyastronomy.readthedocs.io) >= 0.18.1
- [SciPy](https://scipy.org) >= 1.10.1
- [Specutils](https://specutils.readthedocs.io) >= 1.3.1

## Main Functions

- **find2c**

> This task can detect the best possible solution of mass ratio and effective temperature of the secondary companion in a single-lined spectroscopic binary, for a given spectra dataset.
> As the radial velocity (*RV*) for a secondary companion is unknown (single-line spectroscopic binary), a grid of different possible mass ratios (*q*) must be defined. The script uses the *RV* of the primary component to estimate the *RV*-value of the secondary component. For each mass ratio *q* the script uses de RV of the primary, and a spectral disentangling technique is applied using the task `spbina` and a mean spectrum corresponding to secondary companion is obtanied. `find2c` convolves the mean spectrum with a synthetic template spectra list provided by the user. The reconstructed secondary spectra are correlated against synthetic templates of different temperatures. The cross-correlation function (CCF) maximum is used to evaluate the different values of cross-correlation for each *q* and each effective temperature *T<sub>eff</sub>* are obtained. Each CCF value corresponds to the best synthetic template and this could be analyzed in a 3D density diagram, where X-axis is the mass ratio and Y-axis the effective temperature. 
> Mandatory parameters:
> - `lis`: file list of observed spectra to process. The file list is a text file containing a list of spectra images for input (in FITS format). The IRAF terminology to specify this kind of file as the form "at file" is applied and it means the file name must be preceded by the symbol @. The easiest way to generate the file list is using the LINUX command **ls** (string);
> - `tmp`: full path to the folder containing spectra templates (string);
> - `vgamma`: estimated value for systemic radial velocity in km/s (float).
> Optional parameters:
> - `spa`: output name of primary component mean spectrum (string);
> - `spb`: output name of secondary component mean spectrum (string);
> - `qmin`: minimum mass ratio for cross-correlation grid (float);
> - `qmax`: maximum mass ratio for cross-correlation grid (float);    
> - `deltaq`: mass increments for cross-correlation grid (float);
> - `wreg`: spectral regions for cross-correlation analysis (string). The selected region is specified among "-" and the different regiones joined with ' " , ";
> - `nproc`: number of CPU cores to use in computing; it depends of computer resources (integer).
> 
> Example:
```python3
find2c('@lista', '/home/user/templates', wreg='4000-4320,4360-4850,4875-5900', nproc=8)
```

- **hselect**
> Extract keyword values (`fields`, string viariable type) from a FITS image or file list (`img`, string viariable type). In case of file list, the file name must be preceded by the symbol @.
> 
> Example:
```python3
hselect('@lista', 'object')
```

- **rvbina**
> This task computes the cross-correlation of radial velocities between two spectra using the Fast Fourier Transform method (FFT) to perform the convolution. The output is the full discrete linear convolution of the inputs showed in a dispersion grid defined on the the spectral region of ´wreg´. It requires the previously estimation of primary and secondary mean spectra. First, the primary RV is determined by cross-correlation among observed spectrum and an estimated template for primary component; then the mean spectroscopic features are removed from observed spectra and the secondary RV determination is carry out. The cross-correlation function is fitted for the most significative peak with a Gaussian fitting.
> Mandatory parameters:
> - `lis` (string): file list of observed spectra to process. The file name must be preceded by the symbol @.
> Optional parameters:
> - `spa`: name of primary component mean spectrum (string);
> - `spb`: name of secondary component mean spectrum (string);
> - `ta`: spectrum template, in FITS format with or without extension, for comparison with primary component (string);
> - `tb`: spectrum template, in FITS format with or without extension, for comparison with secondary component (string);
> - `wreg`: spectral regions for cross-correlation analysis (string). The selected region is specified among "-" and the different regiones joined with ' " , ";
> - `aconv`: damping convergence factor (float);
> - `keyjd`: header keyword that contains Julian Date (string);
> - `fitcont`: continuum subtraction of spectra dataset prior to correlation analysis (boolean);
> - `interac`: process cross-correlation interactively (boolean).
> 
> Example:
```python3
rvbina('@lista', spa='sp_A', spb='sp_B', wreg='3550-4700,4850-5680', keyjd='HJD', fitcont=True)
```

- **rvextract**
> Analyze the convergence of iterative **rvbina** task for a spectra file list (mandatory parameter `lis`, string viariable type). The optional parameters are `output` (string) to write a file with RV values and `graphic` (boolean) to show the RV convergence as a function of iteration number.
> 
> Example:
```python3
rvextract('@lista', output='file_RVs.txt', graph=True)
```

- **setrvs**
> Measurement of radial velocities for each spectrum from data-set using cross-correlation with templates defined by the user. This task can be applied for SB1 and SB2 binary system; in case of double-line stars, the template for secondary companion star must also be provided. The only mandatory parameter is spectra file list `lis` (string). 
> Optional parameters:
> - `ta`: spectrum template, in FITS format with or without extension, for comparison with primary component (string);
> - `tb`: spectrum template, in FITS format with or without extension, for comparison with secondary component (string);
> - `wreg`: spectral regions for cross-correlation analysis (string). The selected region is specified among '-' and the different regiones joined with ','; example: 4000-4090,4110-4320,4360-4850,4875-5290,5350-5900;
> - `keyjd`: header keyword that contains Julian Date (string);
> - `fitcont`: continuum subtraction of spectra dataset prior to correlation analysis (boolean);
> - `interac`: process cross-correlation interactively (boolean).
> 
> Example:
```python3
setrvs('@lista', ta='/home/user/templates/04800-4.50.fits', interac=True)
```

- **spbina**
> Compute spectral disentangling. The only mandatory parameter is spectra file list `lis` (string).
> Optional parameters:
> - `spa`: name of primary component mean spectrum (string);
> - `spb`: name of secondary component mean spectrum (string);
> - `nit`: number of iterations (integer);
> - `frat`: estimated flux ratio among components (float);
> - `reject`: reject pixels using a sigma clipping algorithm (boolean);
> - `q`: mass ratio among components (if the RV of secondary component is unknown, float variable type);
> - `vgamma`: estimated value for systemic radial velocity in km/s (float);
> - `obspha`: calculate spectra for all phases (boolean);
> - `showtit`: enable user interface (boolean).
> 
> Example for binary type SB2:
```python3
setrvs('@lista', nit=10, frat=0.67, interac=True)
```
> 
> Example for binary type SB1 (radial velocity for secondary companion unknown with determined mass ratio `q`):
```python3
setrvs('@lista', nit=10, frat=0.67, q=0.81, vgamma=2.1, interac=True)
```

- **splot**
> Plot and show spectrum (must be in FITS format). The mandatory parameter is the file spectrum name `file` (string).
> Optional parameters:
> - `xmin`: lower wavelength limit for graphic (float);
> - `xmax`: upper wavelength limit for graphic (float);
> - `ymin`: lower flux limit for graphic (float);
> - `ymax`: upper flux limit for graphic (float);
> - `scale`: flux scale factor (float);
> - `markpix`: mark flux pixel values (boolean);  
> - `newfig`: open spectrum in a new window (boolean);
> - `color`: spectrum graphic color (see **matplotlib.pyplot** for colors availables, string variable type).
> 
> Example for show three spectra in same figure:
```python3
splot('2023-01-05_2350.fits', xmin=3215, xmax=5500)
splot('2022-10-14_0220.fits', xmin=3215, xmax=5500, newfig=False, color='blue')
splot('2022-08-22_0310.fits', xmin=3215, xmax=5500, newfig=False, color='green')
```

- **uniform**
> Explore a spectra data-set (if interactive mode is active, given by the optional parameter `interac=True`) and select the bests spectra to process. Finally, this task scales each spectrum from selected data-set to a mean continuum factor to make them morphologically consistent. The mandatory parameter is `lis` (string) for the file list of observed spectra to process; and optional parameter is interactive mode `interac` (boolean). 
> 
> Example for interactive mode useful:
```python3
uniform('@lista',interac=True)
```

- **vexplore**
> Show and explore results for systemic radial velocities grid analysis (only for previously **vgrid** task running). The mandatory parameter `folder` (string) is the folder output name obtained from **vgrid** task. The different cross-correlation analysis results for the systemic velocities grid can be checked with a interactive 3-D graphic.
>
> Example:
```python3
vexplore('output_00/')
```

- **vgrid**
> When the orbital period for primary RV values cannot be derived from the primary RVs, and the systemic velocity is unknown, the task `vgrid` can be applied to estimate the best systemic radial velocity of binary system through a systemic velocities grid around the most probable value. 
>Mandatory parameters:
> - `lis`: file list of observed spectra to process (string);
> - `tmp`: full path to the folder containing spectra templates (string);
> Optional parameters:
> - `svmin`: lower systemic radial velocity for grid in km/s (float);
> - `svmax`: upper systemic radial velocity for grid in km/s (float);
> - `step`: radial velocity step for meshing grid in km/s (float);
> - `qmin`: minimum mass ratio for cross-correlation grid (float);
> - `qmax`: maximum mass ratio for cross-correlation grid (float);     
> - `deltaq`: mass increments for cross-correlation grid (float);
> - `wreg`: spectral regions for cross-correlation analysis (string). The selected region is specified among '-' and the different regiones joined with ','; example: 4000-4090,4110-4320,4360-4850,4875-5290,5350-5900;
> - `nproc`: number of CPU cores to use in computing; it depends of computer resources (integer).
>
> Example:
```python3
vgrid('@lista', '/home/user/templates', svmin=4.6, svmax=8.1, step=0.2, qmin=0.1, qmax=0.45, nproc=8)
```
