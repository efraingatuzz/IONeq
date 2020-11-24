# IONeq

IONeq is X-ray high-resolution photoabsorption model which compute the absorption coefficient assuming ionization equilibrium. The model takes into account turbulent broadening.  A complete description of the science behind the model is described in [Gatuzz & Churazov (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..696G/abstract).

OBTAINING IONeq

The model can be downloaded from the Github repository at https://github.com/efraingatuzz/IONeq. The contents of the folder include:

- atomic_data/AtomicData.fits -- atomic database binary fits file. This must reside in the directory atomic_data inside the folder where the model is located.
- ioneq.f90 -- source code for IONeq
- rates.com, rates.h, com-auger.com, com-file.com -- global variables definition
- tab-auger-Wijk.dat -- Auger probabilities
- lmodel.dat, lmodel_ioneq.dat -- local model definition file needed by xspec.
- compile.sh -- installation script written on bash.
- README.md -- this file

INSTALLATION

You can use the compile.sh file to install the model by doing

sh compile.sh

In the  model folder or you can setting up and using this model is as described in the xspec manual:

0) You need to have the heasoft package installed on your machine, but it must be built from source. Local models cannot be installed from a binary installation.

1) untar this directory somewhere in your user area

2) setup your headas environment (eg. 'setenv HEADAS /path/to/architecture',and 'source \$HEADAS/headas-init.csh')

3) start up xspec, and in response to the prompt type 

'initpackage ioneq lmodel_ioneq.dat <path-to-current-directory>',

where <path-to-current-directory> is the full path of the current directory. After the build is complete type 

'lmod ioneq <path-to-current-directory>'

In subsequent  sessions you don't neet to do the initpackage step again, just the lmod.
  
PARAMETERS

The input parameters include the hydrogen column density (in units of 1022 cm-2), the temperature (log T), the ionization parameter (log ξ), the O, Fe and Ne atomic abundances (relative to Grevesse & Sauval, 1998), the Fe metallic iron abundance and the turbulence broadening (in units of km/s). The redshift is also a model parameter.

Note that the IONeq does not solve the thermal equilibrium equations, that is the balance heating rate=cooling rate, from which one obtains the gas temperature. Instead, we treat both, the gas temperature T and the ionization parameter ξ, as independent (free) parameters of the model.

ATOMIC DATA

With the default set up - that is, if you have run compile.sh - the model will look for the cross-section data file in atomic_data/AtomicData.fits, relative to the directory in which the module is located. The XSPEC xset command can be used to set the IONEQROOT variable; if this is set then it is used instead of the path to the module. So after

XSPEC12> xset IONEQROOT /path-to-folder/ioneq/

then the model will use the file atomic_data/AtomicData.fits (the IONEQROOT refers to the directory containing the atomic_data/ directory). Note that IONEQROOT over-rides any changes made by running compile.sh when building the model.

The location of the file can be found by setting the XSPEC chatter level to 20 or higher - e.g.

XSPEC12> chatter 20

before evaluating the model.

We include in our model the following atomic data, in order to compute the ion fractions assuming ionization equilibrium:

- Collisional ionization: ionization rates from [Voronov (1997)](https://ui.adsabs.harvard.edu/abs/1997ADNDT..65....1V/abstract).
- Radiative recombination: rates from [Verner & Ferland (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract).
- Dielectronic recombination: rates from [Arnaud & Rothenflug (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26AS...60..425A/abstract).
- Photoionization: cross sections from [Verner et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract) and Auger probabilities from [Kaastra & Mewe (1993)](https://ui.adsabs.harvard.edu/abs/1993A%26AS...97..443K/abstract). 

NOTE: In order to compute the optical depth, IONeq uses the same high-resolution photoabsorption cross-sections included in the XSTAR photonionization code.
 
CONTACT

This package is still being tested. Please contact me with any reports or questions.

egatuzz@mpe.mpg.de
    
New in version 1.3 (October 2020):
- Abundances are read from XSPEC 
- The model relies on the fftw library provided by HEASOFT

New in version 1.2 (January 2019): 
 - Now the energy grid covers until 1e6 eV
 - Now when doing the interpolation of the photoabsorption-cross sections 
 - We are using Verner cross-sections for high energy

 New in version 1.1 (May 2019): 
- This version includes H,He,C,N,O,Ne,Mg,Ne,S,Si,Ar,Ca,Fe,Ni and Zn  from opacity project
