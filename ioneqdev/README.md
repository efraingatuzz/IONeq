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

First, you need to decompress the "AtomicData.tar.gz" file inside the "atomic_data" folder (e.g. using untar). Then, you can use the compile.sh file to install the model by doing (depending on your OS)

sh compile_linux.sh/sh compile_mac.sh

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


CONTACT

This package is still being tested. Please contact me with any reports or questions.

egatuzz@mpe.mpg.de
    
New in version 1.3 (October 2020):
- Abundances are read from XSPEC 
- The model relies on the fftw library provided by HEASOFT

New in version 1.2 (January 2019): 
 - Now the energy grid covers until 1e6 eV
 
 New in version 1.1 (May 2019): 
- Only cross-sections for ions with relative abundance >1e-3 are interpolated (to increase speed)
