# PREFLOW manual
<!-- #### September 4, 2021 -->
<!-- #### Tenyo Kawamura (Kavli IPMU, University of Tokyo) -->

## Introduction
**PREFLOW** is an X-ray timing model for black hole binaries (BHBs), which calculates both power spectrum and cross spectrum.
It takes into account the two key processes expected to happen in the accretion flow, i.e., 
propagation of mass accretion rate fluctuations and X-ray reflection.

Despite being a timing model, **PREFLOW** can work on the spectral fitting package XSPEC.
Since the model is currently written in Python, users can use the **PREFLOW** model on PyXSPEC, not on standard XSPEC.
This manual describes the installation of **PREFLOW** to PyXSPEC, its use in PyXSPEC, and model parameters.
To perform fitting, timing data also need to be read by PyXSPEC.
Thus, the manual contains an explanation of how to make PyXSPEC read timing data.

## About PREFLOW
The propagation of mass accretion rate fluctuations and the X-ray reflection to the accretion flow 
are expected to be responsible for the fast X-ray variability observed in BHBs.
Taking account of these processes, the **PREFLOW** model is designed 
not only to derive the accretion flow geometry and viscous time-scales 
but also to measure energy spectra for various components.
Note that the **PREFLOW** does not model QPOs.
Here, the information which should help users to try the model is briefly summarized.
See [Kawamura et al. 2021](https://arxiv.org/abs/2107.12517) (submitted to MNRAS) 
for the details of modeling the variability.

The **PREFLOW** model makes the following assumptions:
- The truncated disk / hot inner flow geometry.
- Mass accretion rate fluctuations propagates inwards 
	between the inner region of the truncated disk, **rout**--**rds**, 
	and the entire hot flow, **rds**--**rin**.
- The hot inner flow is inhomogeneous and described by two Comptonization components.
	They are named the soft Comptonization and hard Comptonization 
	for the outer (**rds**-**rsh**) and inner (**rsh**-**rin**) components.
	
In summary, the variable accretion flow consists of three spectral components.

| Region              | Spectral component        |
| ----                | ----                      |
| **rout**--**rds**   | (Variable) disk           |
| **rds**--**rsh**    | Soft Comptonization       |
| **rsh**--**rin**    | Hard Comptonization       |

The **PREFLOW** model also assumes 
the soft Comptonization and hard Comptonization illuminate the accretion flow, 
resulting in the reflection.
The reflection spectra are referred to as the soft reflection and hard reflection, respectively.
The impulse response is simply assumed to be a top-hat function with two parameters, 
**tref** and **dtref**.

## Installation to PyXSPEC
This model installation instruction assumes PyXSPEC is already installed in users' computer
(see [Build and Install PyXspec](https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/buildinstall.html) 
for the installation of PyXSPEC).

1. Download all the files in [the **PREFLOW** repository](https://github.com/tenyokawamura/PREFLOW) 
	 by clicking the `Code` button and following instructions.
	 The following files should be contained:
	- preflow.py
	- preflow_h.py
	- fourier_h.py
	- lmodel_preflow.py
	- test_preflow.py
	- quick_plot.py
	- README.py

2. Put all these files in a clean directory on users' computer.

3. Go to the directory.

		>>> cd /path/to/the/directory/

4. Launch an interactive PyXSPEC session:

		>>> python
		>>> import xspec

5. Import the source codes of the **PREFLOW** model:

		>>> from preflow import *
		>>> from lmodel_preflow import *

6. Add the **PREFLOW** model to PyXSPEC models:

		>>> xspec.AllModels.addPyMod(preflow, lmodel_preflow_info, 'add')

Then users can use the model **PREFLOW** in the same way as other additive models.
But keep in mind that **PREFLOW** is a *timing* model, not a spectral model, 
so using **PREFLOW** together with spectral models to define a single model should be meaningless.

## How to use the model on PyXSPEC
In the use of **PREFLOW** on PyXSPEC, 
users need to regard energy (keV) as frequency (Hz) for X-axis and 
regard flux (photons/cm^2/s/keV) as power spectrum ((rms/mean)^2/Hz), 
cross (real, imaginary, amplitude ) spectrum ((rms/mean)^2/Hz), 
phase lag (rad), or time lag (s) for Y-axis, depending on which quantity is calculated.

After the installation of the **PREFLOW** model, 
users can plot a power spectrum with following commands:

		>>> xspec.AllData.dummyrsp(0.1, 50., 1000, scaleType='lin')
		>>> xspec.Model('preflow')
		>>> xspec.Plot.device='/xw' # change a plot device if necessary.
		>>> xspec.Plot('emo')
	
or equivalently,

		>>> python -i quick_plot.py`

The power spectrum multiplied by the frequency in units of (rms/mean)^2 is plotted with the 'emo' keyword.

The timing quantity calculated is switched with the parameter **quant** (see 'Parameters' for further details ).
For example, the following commands plot a time lag spectrum.

		>>> xspec.Model('preflow', setPars={33:6})
		>>> xspec.Plot('mo')

Since the time lag spectrum can be negative, it may be more useful to plot the time lag spectrum on a linear scale.

		>>> xspec.Plot.iplot('mo')
		>>> log y off
		>>> r y -0.1

## Parameters
| Number | Parameter     | Units      | Definition                                                                              | Comments |
| ----   | ----          | ----       | ----                                                                                    | ----     |
| 1      | **Mass**      | Solar mass | BH mass.                                                                                |          |
| 2      | **rin**       | Rg         | Inner radius of the hard Comptonization.                                                |          |
| 3      | **rsh**       | Rg         | Transition radius between the hard and soft Comptonization.                             |          |
| 4      | **rds**       | Rg         | Transition radius between the soft Comptonization and variable disk.                    |          |
| 5      | **rout**      | Rg         | Outer radius of the variable disk.                                                      |          |
| 6      | **Nring**     | -          | Number of rings separating the variable flow.                                           |          |
| 7      | **tref**      | sec        | Start time of the reflection impulse response.                                          |          |
| 8      | **dtref**     | sec        | Time width of the reflection impulse response.                                          |          |
| 9      | **Fvar**      | -          | Fractional variability of the mass accretion rate in radial decade.                     |          |
| 10     | **Bdisk**     | -          | Coefficient of the viscous frequency in the variable disk.                              |          |
| 11     | **mdisk**     | -          | Radial index of the viscous frequency in the variable disk.                             |          |
| 12     | **Bflow**     | -          | Coefficient of the viscous frequency in the hot flow.                                   |          |
| 13     | **mflow**     | -          | Radial index of the viscous frequency in the hot flow.                                  |          |
| 14     | **stress**    | -          | Selection of the boundary condition of the emissivity. 1: stressed, 2: stress-free.     | It should be frozen.         |
| 15     | **gamma**     | -          | Radial index of the emissivity.                                                         |          |
| 16     | **Emin**      | keV        | Lower bound of the energy band.                                                         | Display purpose only. It should be always frozen.         |
| 17     | **Emax**      | keV        | Upper bound of the energy band.                                                         | Display purpose only. It should be always frozen.         |
| 18     | **Fdisk**     | -          | Fraction of the variable disk count rate in the energy band.                            | Detector effective area should be taken into account.         |
| 19     | **Fscp**      | -          | Fraction of the soft Comptonization count rate in the energy band.                      |          |
| 20     | **Fhcp**      | -          | Fraction of the hard Comptonization count rate in the energy band.                      |          |
| 21     | **Fsref**     | -          | Fraction of the soft reflection count rate in the energy band.                          |          |
| 22     | **Fhref**     | -          | Fraction of the hard reflection count rate in the energy band.                          |          |
| 23     | **Eminr**     | keV        | Lower bound of the reference band.                                                      | Display purpose only. It should be always frozen.         |
| 24     | **Emaxr**     | keV        | Upper bound of the reference band.                                                      | Display purpose only. It should be always frozen.         |
| 25     | **Fdiskr**    | -          | Fraction of the variable disk count rate in the reference band.                         |          |
| 26     | **Fscpr**     | -          | Fraction of the soft Comptonization count rate in the reference band.                   |          |
| 27     | **Fhcpr**     | -          | Fraction of the hard Comptonization count rate in the reference band.                   |          |
| 28     | **Fsrefr**    | -          | Fraction of the soft reflection count rate in the reference band.                       |          |
| 29     | **Fhrefr**    | -          | Fraction of the hard reflection count rate in the reference band.                       |          |
| 30     | **Eminrr**    | keV        | Lower bound of the reference band for the reflection.                                   | Display purpose only. It should be always frozen.         |
| 31     | **Emaxrr**    | keV        | Upper bound of the reference band for the reflection.                                   | Display purpose only. It should be always frozen.         |
| 32     | **Fscprr**    | -          | Fraction of the soft Comptonization count rate in the reference band for the reflection.|          |
| 33     | **quant**     | -          | Selection of timing quantity calculated.                                                | It should be always frozen. See 'Additional comments' for further details.          |
| 34     | **print**     | -          | Selection of display. 1: details of calculations are displayed, 2: not displayed.       | It should be always frozen and set to 2 in fitting to avoid displaying too much information.        |
| 35     | **norm**      | -          | Normalization.                                                                          | It should be frozen. See 'Additional comments' for further details.         |

### Additional comments
- Truncated disk model is employed as the accretion flow geometry.
	The geometric parameters must satisfy **rin**<=**rsh**<=**rds**<=**rout**.
- For both the energy band and reference band, the sum of the fraction of count rate for each spectral component 
	must not be larger than 1 (**Fdisk(r)**+**Fscp(r)**+**Fhcp(r)**+**Fsref(r)**+**Fhref(r)**<=1). 
	The rest of the fraction is assumed to come from the constant, invariable component 
	(**Fconst(r)**=1-(**Fdisk(r)**+**Fscp(r)**+**Fhcp(r)**+**Fsref(r)**+**Fhref(r)**)).
- Count rate should be convolved with a detector effective area (not in units of \*\*/cm^2) 
	since the variability is calculated for the number of photons *detected* (not emitted from a source).
- Reference band for the reflection is the energy band, 
	which is taken into account to calculate the variability of X-ray flux illuminating the disk.
	All the illuminating X-rays are assumed to come from the soft and hard Comptonization.
	For the reference band of reflection, 
	the sum of fraction of the count rate for the soft and hard Comptonization corresponds to 1, 
	which means **F_scprr**<=1.
- Impulse response of the reflection is simply described with a top-hat function.
	**tref** is the time at the rising edge, **dtref** is the duration of the hat.
- Photon energy parameters (**Emin**, **Emax**, **Eminr**, **Emaxr**, **Eminrr**, **Emaxrr**) 
	are just used to help users to know the information.
	The information on energy bands should be helpful 
	when users try to fit timing data for different energy bands simultaneously.
- Normalization parameter **norm** provided by XSPEC should be frozen since 
	- **Fvar** is responsible for the amplitude for 
		power spectrum, cross real spectrum, cross imaginary spectrum and cross amplitude spectrum, 
	- Many parameters, such as viscous frequency prescription and accretion flow geometry, 
		determines the amplitude for
		phase lag spectrum and time lag spectrum.
- Parameter **quant** is used for selection of timing quantity calculated:
	| Value  | Quantity calculated       | Units           |
	| ----   | ----                      | ----            |
	| 1      | Power spectrum            | (rms/mean)^2/Hz |
	| 2      | Cross real spectrum       | (rms/mean)^2/Hz |
	| 3      | Cross imaginary spectrum  | (rms/mean)^2/Hz |
	| 4      | Cross amplitude spectrum  | (rms/mean)^2/Hz |
	| 5      | Phase lag spectrum        | radian          |
	| 6      | Time lag spectrum         | sec             |
- The PREFLOW model has a number of parameters, which makes fitting difficult.
	Users are recommended to freeze the fractions of count rate for the reference band and the reference band for the reflection 
	by using reliable spectral decomposition or simply assuming spectra.
- Depending on the energy bands used, a fraction of the count rate may not be independent of another fraction of the count rate.
	Users should keep this caution in mind in performing fit.
- If the timing quantity to be calculated is a power spectrum (**quant**=1), 
	the parameters on the reference band are not used.

## Timing data
To fit timing data with the **PREFLOW** model on PyXSPEC, the timing data need to be read by PyXSPEC.
Thus, users need to prepare the data file (.pha) and response files (.rmf and .arf) so that PyXSPEC can read the timing data.

### PHA file
Under construction.

### RMF file
Under construction.

### ARF file
Under construction.


## Additional notes
### Computation efficiency
The model calculation efficiency is determined by 
the number of rings separating the accretion flow, **Nring**, and 
the frequency range of timing quantities.

**Nring** determines the resolution of the model calculation.
If **Nring** is too small, oscillation is seen at high frequencies in timing quantities.
Users are recommended to set **Nring** to ~20 to make a rough estimate and to ~50 to perform the fitting.

When PyXSPEC has timing data, 
the frequency range, over which timing quantities are calculated, 
is determined by the frequency range of the timing data.
Unfortunately, ignoring some frequencies by using the `ignore` method 
does *not* change the frequency range of the model calculation.
Thus, it may make fitting faster to remove unnecessary data from the timing data file.

