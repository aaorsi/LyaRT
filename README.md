Monte Carlo radiative transfer of Ly-alpha photons 

LyaRT
=====


## Table of contents
* [Description](#description)
* [Getting starterd](#getting-started)
* [Running example](#running-example)
* [Model description](#model-descrption)
* [Varying parameters](#varying-parameters)

## Description

This code performs a Monte Carlo calculation of the radiative transfer of Lyman-alpha photons as they diffuse through a cloud of neutral gas. The code is written in `C`, but there is code for running and validating the results written in both `Python`and `ÌDL`. This code was developed as part of my PhD thesis, and it has been used in numerous publications: [Orsi et al. 2012](http://adsabs.harvard.edu/abs/2012MNRAS.425...87O), [Orsi et al., 2016](http://adsabs.harvard.edu/abs/2016MNRAS.456.3827O), [Bielby et al., 2016](http://adsabs.harvard.edu/abs/2016MNRAS.456.4061B), [Gurung-Lopez et al. 2018](http://adsabs.harvard.edu/abs/2018arXiv180700006G) and [Gurung-Lopez et al.2019](http://adsabs.harvard.edu/abs/2018arXiv181109630G). 

## Getting started

This code is pretty straightforward to install and use. After cloning the repository, add it to your python path (e.g. `export PYTHONPATH="${PYTHONPATH}:/home/your-path-to-get_emlines/"`) and use it as in the example below.

### Dependencies

Make sure you have the following packages installed:
- numpy
- matplotlib
- scipy

The grid dataset is included in this repository.

## Running example
```python
import get_emlines as lines
sfr  = 0.1   # Msun/yr
Z    = 0.01  # M_metals/M_gas 
lums = lines.get_emlines(sfr,Z)

# OIII 5007 luminosity:
print lums['OIII_5007']

# All luminosities available should retrieve something like this:
lums.dtype
Out[3]: dtype([('Lyalpha', '<f4'), ('Hbeta', '<f4'), ('Halpha', '<f4'), ('OII_3727', '<f4'), ('OII_3729', '<f4'), 
('OIII_5007', '<f4'), ('OIII_4959', '<f4'), ('OI_6300', '<f4'), ('NII_6548', '<f4'), ('NII_6584', '<f4'), 
('SII_6717', '<f4'), ('SII_6731', '<f4'), ('NeIII_3870', '<f4'), ('CII_158um', '<f4'), ('NII_205um', '<f4')])

# Another example with an array of dummy galaxies:
import numpy as np
sfr = np.logspace(-2,2,1e3) ; z = np.linspace(1e-3,1e-1,1e3)
# Testing the verbose keyword as well
lums = lines.get_emlines(sfr,z,verbose=True)

```
The messages printed in the above example (verbose activated) should look something like this:
```

 DEBUG - get_emlines(): verbose output activated
 DEBUG - get_lumlines(): Rootdir:
 INFO - get_lumlines(): Computing emission lines for 1000 galaxies
 DEBUG - get_2dfunc(): Constructing grid function for all lines
 DEBUG - get_2dfunc(): Done constructing 2D grid function
 DEBUG - get_2dfunc(): Constructing grid function for Halpha
 DEBUG - get_2dfunc(): Done constructing 2D grid function
 DEBUG - get_lumlines(): Running lines with g0=-1.300000
 INFO - get_emlines(): Luminosities computed OK

```

```python
# some output luminosities:
print lums['Halpha'].min(), lums['Hbeta'].max()
39.2671 42.8018
```

You can access the list of lines available in, e.g.: 
`HIImodels/Lines/LineInfo_Levesque10`
By default, the list of lines available (and their names) are:

```
# Emission Lines contained in file LineData_Levesque10
# id	Emission Line			Central Wavelength
0  Lyalpha			1216.0
1  Hbeta			4861.0
2  Halpha			6563.0
3  OII_3727			3727.0
4  OII_3729			3729.0
5  OIII_5007			5007.0
6  OIII_4959			4959.0
7  OI_6300			6300.0
8  NII_6548			6548.0
9  NII_6584			6584.0
10  SII_6717			6717.0
11  SII_6731			6731.0
12  NeIII_3870			3870.0
13  CII_158um			1.5800
14 NII_205um 2.0500
```
The central wavelgnth of the last two FIR lines is in um, the rest in Angstroms.


## Model description

This code computes galaxy line luminosities *log(L [erg s-1])* based on an input *SFR [M_sun/yr]* and metallicity *Z*. To assign emission lines to galaxies, first we assign an ionization parameter *q*, which is assumed to be related to the gas-phase metallicity by a power-law:

![Alt Text](https://github.com/aaorsi/get_emlines/blob/master/eq_gif.gif)

Then the code makes use of the grid of photo-ionization models of [Levesque et al. 2010](https://www.emlevesque.com/model-grids/) to perform a bilinear interpolation over *q* and *Z*, leaving the electron density fixed. Finally, to scale line fluxes to galaxy-wide luminosities, the code uses the input *SFR* to infer the *H-alpha* luminosity, and then use its predicted flux to infer all other line luminosities.


## Varying parameters


. You can change the slope `g0` or normalization `q0` from their standard values (from Orsi+14) by doing:
```python
lum_disk = lines.get_emlines(sfrdisk, zdisk, g0 = -1.5,q0 = 3.5e7)
```
In the above the slope and normalization changed to `-1.5` and `3.5e7`, respectively.




The code is being adapted to be released together with python routines to handle it. 


19/10/14:
The code might suffer some major changes, major debugging taking place

