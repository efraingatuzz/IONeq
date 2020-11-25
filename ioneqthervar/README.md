# IONeqthervar

IONeq is an X-ray high-resolution photoabsorption model which compute the absorption coefficient assuming ionization equilibrium. The model takes into account turbulent broadening.  A complete description of the science behind the model is described in [Gatuzz & Churazov (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..696G/abstract). Note that the IONeq does not solve the thermal equilibrium equations, that is the balance heating rate=cooling rate, from which one obtains the gas temperature. Instead, the gas temperature T is as independent (free) parameters of the model.

The IONeqthervar includes the elemental abundances as free parameters in the model

OBTAINING IONeqthervar

The model can be downloaded from the Github repository at https://github.com/efraingatuzz/IONeq. The contents of the folder include:

- atomic_data/AtomicData.fits -- atomic database binary .fits file. This must reside in the directory atomic_data inside the folder where the model is located
- ioneqthervar.f90 -- source code for IONeqthervar
- rates.com, rates.h, com-auger.com, com-file.com -- global variables definition
- tab-auger-Wijk.dat -- Auger probabilities
- lmodel_ioneq.dat -- local model definition file needed by xspec
- compile_linux.sh -- installation script written in bash for LINUX
- compile_mac.sh -- installation script written in bash for MAC

INSTALLATION

You can use the compile_linux(mac).sh file to install the model by doing

sh compile_linux(mac).sh

In the  model folder or you can setting up and using this model is as described in the xspec manual:

0) You need to have the heasoft package installed on your machine, but it must be built from source. Local models cannot be installed from a binary installation.

1) untar this directory somewhere in your user area

2) setup your headas environment (eg. 'setenv HEADAS /path/to/architecture',and 'source \$HEADAS/headas-init.csh')

3) start up xspec, and in response to the prompt type 

'initpackage ioneqthervar lmodel_ioneqthervar.dat <path-to-current-directory>',

where <path-to-current-directory> is the full path of the current directory. After the build is complete type 

'lmod ioneqthervar <path-to-current-directory>'

In subsequent  sessions you don't neet to do the initpackage step again, just the lmod.

PARAMETERS

The input parameters include the hydrogen column density (in units of 10^22 cm-2), the temperature (log T), the Fe metallic iron abundance, the elemental abundances (relative to the solar abundances set in XSPEC), the turbulence broadening (in units of km/s) and the redshift (z)

CONTACT

This package is still being tested. Please contact me with any reports or questions.

egatuzz@mpe.mpg.de









