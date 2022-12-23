# PREFLOW manual
<!-- #### December 21, 2022 -->
<!-- #### Tenyo Kawamura (Kavli IPMU, University of Tokyo) -->

## Introduction
**PREFLOW** is an X-ray timing model for black hole binaries, which calculates both power spectra and cross spectra.
The model is based on the propagation of mass accretion rate fluctuations but accounts for reverberation and spectral pivoting.
By explicitly connecting time-averaged emission with rapid variability, **PREFLOW** can be used to a spectral-timing analysis.

Despite being a timing model, **PREFLOW** can work on the spectral fitting package **XSPEC**.
Since the model is currently written in Python, users can use the **PREFLOW** model on **PyXSPEC**, not on standard **XSPEC**.
This manual describes the installation of **PREFLOW** to **PyXSPEC**, its use in **PyXSPEC**, and model parameters.
In fitting, timing data also need to be read by **PyXSPEC**.
Thus, the manual instructs how to format timing data to be loaded on **PyXSPEC**.

## About PREFLOW
The propagation of mass accretion rate fluctuations is expected to be an underlying process driving rapid X-ray variability.
Fundamentally, this process describes the variability of the mass accretion rate, but what we can observe are X-ray photons.
In the conversion from the mass accretion rate into X-ray radiation, spectral pivoting and reverberation become relevant.
Given that these processes are not independent of the propagating fluctuations process, correlations between these processes affect the fast variability.
These connections have motivated us to account for these three processes altogether in **PREFLOW**.
Here, we give an instructive description intended for the practical use of the model. 
The physical description is found at [Kawamura et al. 2022b](https://arxiv.org/abs/2209.14492) (See also [Kawamura et al. 2022a](https://arxiv.org/abs/2107.12517) for the previous version of the model).

The **PREFLOW** model makes the following assumptions:
- The truncated disk / hot inner flow geometry.
- Mass accretion rate fluctuations can propagate inwards 
	across the inner region of the truncated disk, **rout**-**rds**, 
	and the hot inner flow, **rds**-**rin**.
- The hot inner flow can be inhomogeneous and described by two Comptonization components (at maximum).
	They are named soft Comptonization and hard Comptonization 
	for the outer (**rds**-**rsh**) and inner (**rsh**-**rin**) components.
	
In summary, the variable accretion flow consists of three spectral components at maximum.

| Region              | Spectral component        |
| ----                | ----                      |
| **rout**-**rds**    | (Variable) disk           |
| **rds**-**rsh**     | Soft Comptonization       |
| **rsh**-**rin**     | Hard Comptonization       |

The number of spectral regions can be reduced by setting model parameters appropriately.

The **PREFLOW** model also assumes the soft Comptonization and hard Comptonization illuminate the accretion disc, resulting in reverberation.
The reflection spectra are referred to as soft reflection and hard reflection, respectively.
The impulse response is simply assumed to be a top-hat function with two parameters, **t0** and **dt0**.

The contribution of each spectral region to flux can be given in several ways, which makes **PREFLOW** have some flavors.
The simplest way is to specify it directly.
In this case, the model named **preflow** is used for timing fits.
We can also compute it by explicitly assuming spectral models.
**preflows** adopts **cutoffpl** and **relxill** for Comptonization and reflection, respectively, while **preflowsCp** uses **nthcomp** and **relxillCp**.
Both of them model the disk emission with **diskbb** in common.
**PREFLOWS** and **PREFLOWSCP** calculate energy spectra as well as power spectra and cross spectram, thus can be used for spectral-timing fits.
We summarize the features of these flavors below:

| Model               | Disk spectrum    | Comptonization spectrum | Reflection spectrum |
| ----                | ----             | ----                    | ----                |
| **preflow**         | -                | -                       | -                  |
| **preflows**        | **diskbb**       | **cutoffpl**            | **relxill**         |
| **preflowsCp**      | **diskbb**       | **nthcomp**             | **relxillCp**       |


## Installation to **PyXSPEC**
This model installation instruction assumes **PyXSPEC** is already installed in users' computer
(see [Build and Install PyXspec](https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/buildinstall.html) 
for the installation of **PyXSPEC**).
The **relxill** package needs to be installed for the use of **preflows** and **preflowsCp**, though it does not for the use of **preflow**
(see [Download and Installation](http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/) for the installation of **relxill**).
Set the **RELXILL_TABLE_PATH** environment variable in users', following the installation guide of **relxill**.

1. Download all the files in [the **PREFLOW** repository](https://github.com/tenyokawamura/PREFLOW) 
	 by clicking the `Code` button and following the instructions.
	 The following files should be contained:
	- preflow.py / preflows.py / preflowscp.py
	- preflow_h.py
	- fourier_h.py
	- lmodel_preflow.py
	- test_preflow.py / test_preflows.py / test_preflowscp.py
	- quick_plot_preflow.py / quick_plot_preflows.py / quick_plot_preflowscp.py 
	- README.py

2. Put all these files in a clean directory on users' computer.

3. Go to the directory.

		>>> cd /path/to/the/directory/

4. Launch an interactive **PyXSPEC** session:

		>>> python
		>>> import xspec

5. Import the source codes of the **preflow** model:

		>>> from preflow import *
		>>> from lmodel_preflow import *

6. Add the **preflow** model to **PyXSPEC** models:

		>>> xspec.AllModels.addPyMod(preflow, lmodel_preflow_info, 'add')

Then users can use the model **preflow** in the same way as other additive models.

**preflows** and **preflowsCp** can be activated in the same manner:

		>>> python
		>>> import xspec
		>>> from preflows import *
		>>> from lmodel_preflow import *
		>>> xspec.AllModels.addPyMod(preflows, lmodel_preflows_info, 'add')

for **preflows**, 

		>>> python
		>>> import xspec
		>>> from preflowscp import *
		>>> from lmodel_preflow import *
		>>> xspec.AllModels.addPyMod(preflowscp, lmodel_preflowscp_info, 'add')

for **preflowscp**.

## How to use the model on **PyXSPEC**
In the use of **PREFLOW** on **PyXSPEC**, 
users need to regard energy (keV) as frequency (Hz) for X-axis and 
regard flux (photons/cm^2/s/keV) as power spectrum ((rms/mean)^2/Hz), 
cross (real, imaginary, amplitude ) spectrum ((rms/mean)^2/Hz), 
phase lag (rad), or time lag (s) for Y-axis, depending on which quantity is calculated.

After the installation of the **PREFLOW** model, 
users can plot a power spectrum with following commands:

		>>> xspec.AllData.dummyrsp(0.1, 50., 1000, scaleType='log')
		>>> xspec.Model('preflow')
		>>> xspec.Plot.device='/xw' # change a plot device if necessary.
		>>> xspec.Plot('emo')
	
or equivalently,

		>>> python -i quick_plot.py`

The power spectrum multiplied by the frequency in units of (rms/mean)^2 is plotted with the 'emo' keyword.

The timing quantity calculated is switched with the parameter **quant** (see 'Parameters' for further details ).
For example, the following commands plot a time lag spectrum.

		>>> xspec.Model('preflow', setPars={40:6})
		>>> xspec.Plot('mo')

Since the time lag spectrum can be negative, it may be more useful to plot it on a linear scale.

		>>> xspec.Plot.iplot('mo')
		>>> log y off
		>>> r y -0.1 0.1

Again, **preflows** and **preflowsCp** are used in the same manner.

## Parameters

There are many parameters because the model is still in a phase of active development.
A large fraction of these parameters may be frozen or tied, which should allow users to fit data with an acceptable number of free parameters.

### **preflow**

| Number | Parameter     | Units      | Definition                                                                              | Comments |
| ----   | ----          | ----       | ----                                                                                    | ----     |
| 1      | **Mass**      | Solar mass | BH mass.                                                                                |          |
| 2      | **rin**       | Rg         | Inner radius of the hard Comptonization.                                                |          |
| 3      | **drh**       | Rg         | Radial width of the hard Comptonization.                                                |          |
| 4      | **drs**       | Rg         | Radial width of the soft Comptonization.                                                |          |
| 5      | **drd**       | Rg         | Radial width of the variable disk.                                                      |          |
| 6      | **Nring**     | -          | Number of rings separating the variable flow.                                           | Frozen.  |
| 7      | **Fvd**       | -          | Fractional variability per radial decade for the variable disk.                         |          |
| 8      | **Fvs**       | -          | Fractional variability per radial decade for the soft Comptonization.                   |          |
| 9      | **Fvh**       | -          | Fractional variability per radial decade for the hard Comptonization.                   |          |
| 10     | **Bd**        | -          | Coefficient of the generator frequency in the variable disk.                            |          |
| 11     | **md**        | -          | Radial index of the generator frequency in the variable disk.                           |          |
| 12     | **Bpd**       | -          | Coefficient of the propagation frequency in the variable disk.                          |          |
| 13     | **mpd**       | -          | Radial index of the propagation frequency in the variable disk.                         |          |
| 14     | **Bf**        | -          | Coefficient of the generator frequency in the hot flow.                                 |          |
| 15     | **mf**        | -          | Radial index of the generator frequency in the hot flow.                                |          |
| 16     | **Bpf**       | -          | Coefficient of the propagation frequency in the hot flow.                               |          |
| 17     | **mpf**       | -          | Radial index of the propagation frequency in the hot flow.                              |          |
| 18     | **D**         | -          | Damping parameter.                                                                      |          |
| 19     | **indexd**    | -          | Radial index of the emissivity for the variable disk.                                   |          |
| 20     | **indexf**    | -          | Radial index of the emissivity for the hot flow.                                        |          |
| 21     | **stress**    | -          | Selection of the boundary condition of the emissivity. 1: stressed, 2: stress-free.     | Frozen.  |
| 22     | **Emin**      | keV        | Lower bound of the subject band.                                                        | Frozen (display purpose). |
| 23     | **Emax**      | keV        | Upper bound of the subject band.                                                        | Frozen (display purpose). |
| 24     | **Sd**        | -          | Fraction of the variable disk in the subject band.                                      |          |
| 25     | **Ss**        | -          | Fraction of the soft Comptonization in the subject band.                                |          |
| 26     | **Sh**        | -          | Fraction of the hard Comptonization in the subject band.                                |          |
| 27     | **Ssr**       | -          | Fraction of the soft reflection in the subject band.                                    |          |
| 28     | **Shr**       | -          | Fraction of the hard reflection in the subject band.                                    |          |
| 29     | **Eminr**     | keV        | Lower bound of the reference band.                                                      | Frozen (display purpose). |
| 30     | **Emaxr**     | keV        | Upper bound of the reference band.                                                      | Frozen (display purpose). |
| 31     | **Srd**       | -          | Fraction of the variable disk in the reference band.                                    |          |
| 32     | **Srs**       | -          | Fraction of the soft Comptonization in the reference band.                              |          |
| 33     | **Srh**       | -          | Fraction of the hard Comptonization in the reference band.                              |          |
| 34     | **Srsr**      | -          | Fraction of the soft reflection in the reference band.                                  |          |
| 35     | **Srhr**      | -          | Fraction of the hard reflection in the reference band.                                  |          |
| 36     | **t0s**       | sec        | Start time of the soft reflection impulse response.                                     |          |
| 37     | **dt0s**      | sec        | Time width of the soft reflection impulse response.                                     |          |
| 38     | **t0h**       | sec        | Start time of the hard reflection impulse response.                                     |          |
| 39     | **dt0h**      | sec        | Time width of the hard reflection impulse response.                                     |          |
| 40     | **quant**     | -          | Selection of timing quantity calculated.                                                | Frozen. See 'Additional comments'.|
| 41     | **invert**    | -          | Definition of cross spectrum. 1: positive lag if the reference band is delayed. 2: opposite.       | Frozen.  |
| 42     | **print**     | -          | Selection of display. 1: details of calculations are displayed, 2: not displayed.       | Frozen (display purpose). |
| 43     | **norm**      | -          | Normalization.                                                                          | Frozen (=1.0). See 'Additional comments'.|

### **preflows** / **preflowsCp**

Many parameters are shared with **preflow**.
Spectral parameters are the same as those in spectral models adopted (**diskbb**, **cutoffpl**, **nthcomp**, **relxill**, or **relxillCp**).
Parameters peculiar to **preflows** and **preflowsCp** are those regulating spectral pivoting:

| Number | Parameter     | Units      | Definition                                                                              | Comments |
| ----   | ----          | ----       | ----                                                                                    | ----     |
| 44     | **eta0d**     | -          | Constant term of the sensitivity in the variable disk.                                  |          |                           
| 45     | **eta1d**     | -          | Coefficient of the sensitivity in the variable disk.                                    |          |
| 46     | **eta0s**     | -          | Constant term of the sensitivity in the soft Comptonization.                            |          |                           
| 47     | **eta1s**     | -          | Coefficient of the sensitivity in the soft Comptonization.                              |          |
| 48     | **eta0h**     | -          | Constant term of the sensitivity in the hard Comptonization.                            |          |                           
| 49     | **eta1h**     | -          | Coefficient of the sensitivity in the hard Comptonization.                              |          |

### Additional comments
- For both the subject band and reference band, the sum of the fractions of each spectral component typically corresponds to unity, i.e.,  **Sd**+**Ss**+**Sh**+**Ssr**+**Shr**=1.
	If other stable spectral components exist, the sum is smaller than the unity.
	If spectral pivoting is taken into account, the sum can take an arbitrary value.
- Impulse response of the reflection is simply described with a top-hat function.
	**t0** is the time at the rising edge, **dt0** is the duration of the hat.
- Photon energy parameters (**Emin**, **Emax**, **Eminr**, **Emaxr**, **Eminrr**, **Emaxrr**) 
	are just used to help users to know the information.
- Normalization parameter **norm** provided by **XSPEC** should be fixed to unity since 
	- **Fvd**, **Fvs**, and **Fvh** are responsible for the amplitude for 
		power spectrum, cross real spectrum, cross imaginary spectrum and cross amplitude spectrum, 
	- Many parameters, such as viscous frequency prescription and accretion flow geometry, 
		determines the amplitude for
		phase lag spectrum and time lag spectrum.
- Parameter **quant** is used for the selection of output quantity:
	| Value  | Quantity calculated       | Units                    |
	| ----   | ----                      | ----                     |
	| 0      | Energy spectrum           | photons cm^-2 s^-1 keV^-1|
	| 1      | Power spectrum            | (rms/mean)^2/Hz          |
	| 2      | Cross real spectrum       | (rms/mean)^2/Hz          |
	| 3      | Cross imaginary spectrum  | (rms/mean)^2/Hz          |
	| 4      | Cross amplitude spectrum  | (rms/mean)^2/Hz          |
	| 5      | Phase lag spectrum        | radian                   |
	| 6      | Time lag spectrum         | sec                      |

	**quant**=0 is available for **preflows** and **preflowsCp**, not for **preflow**.
- If the timing quantity to be calculated is a power spectrum (**quant**=1), 
	the parameters on the reference band are not used.

## Timing data
For fits to timing data with the **PREFLOW** model on **PyXSPEC**, the timing data need to be read by **PyXSPEC**.
**XSPEC** regards any data as normal spectra, as far as it is concerned. 
To load fast variability data such as power spectra and cross spectra in **XSPEC**, users need response files (RMF, ARF) as well as a formatted data file (PHA). 
In PHA format, the Fourier frequency is replaced with the 'channel', which is defined between the minimum and maximum Fourier frequencies of interest. 
The diagonal dummy response (RMF) is created to relate the channel and Fourier frequency in a trivial way. 
A constant effective area is simply given in the ARF file. 
The tool flx2xsp can help to prepare these PHA, RMF, and ARF files.

## Additional notes
### Computation efficiency
The model calculation efficiency is determined by 
the number of rings separating the accretion flow, **Nring**, and 
the frequency range of timing quantities.

**Nring** determines the resolution of the model calculation.
If **Nring** is too small, oscillation is seen at high frequencies in timing quantities.
Users are recommended to set **Nring** to ~20 to make a rough estimate and to ~50 to perform the fitting.

When **PyXSPEC** has timing data, 
the frequency range, over which timing quantities are calculated, 
is determined by the frequency range of the timing data.
Unfortunately, ignoring some frequencies by using the `ignore` method 
does *not* change the frequency range of the model calculation.
Thus, it may make fitting faster to remove unnecessary data from the timing data file.

