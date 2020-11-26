# IONeq

IONeq is an X-ray high-resolution photoabsorption model which compute the absorption coefficient assuming ionization equilibrium. The model takes into account turbulent broadening.  A complete description of the science behind the model is described in [Gatuzz & Churazov (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..696G/abstract).

OBTAINING IONeq

Different flavours of the model can be downloaded from the Github repository at https://github.com/efraingatuzz/IONeq. Please note thaat, because of the large files for the atomic data, you must clone the git repository in your local machine (i.e., if you download a zip file the atomic data will not be included). In order to clone the repository first you must install git and lfs

To install git:

Follow the instructions in https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

To install lfs

1- Navigate through https://git-lfs.github.com/ and click in the download option

2- On your computer, locate and unzip the downloaded file.

3- Open Terminal and change the current working directory into the folder you downloaded and unzipped.

4- To install the file, run the command "./install.sh" (without quotes)

5- Verify that the installation was successful by running the command "git lfs install" (without quotes)

To clone the ISMabs respository:

1- Create a folder in your local machine where the model will be located

2- Run the following command:
git clone https://github.com/efraingatuzz/ISMabs

Currently versions include:

- IONeqther: A version of the model that computes optical depths assuming collisional ionization equilibrium conditions and solar abundances.
- IONeqthervar: Same as "IONthervar" but the elemental abundances are parameters of the model
- ioneqdev: latest develop version of the model. 

ATOMIC DATA 


We include in our model the following atomic data, in order to compute the ion fractions assuming ionization equilibrium:

- Collisional ionization: ionization rates from [Voronov (1997)](https://ui.adsabs.harvard.edu/abs/1997ADNDT..65....1V/abstract).
- Radiative recombination: rates from [Verner & Ferland (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract).
- Dielectronic recombination: rates from [Arnaud & Rothenflug (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26AS...60..425A/abstract).
- Photoionization: cross sections from [Verner et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract) and Auger probabilities from [Kaastra & Mewe (1993)](https://ui.adsabs.harvard.edu/abs/1993A%26AS...97..443K/abstract). 

NOTE: In order to compute the optical depth, IONeq uses the same high-resolution photoabsorption cross-sections included in the XSTAR photonionization code.

CONTACT

Please contact me with any reports or questions.

egatuzz@mpe.mpg.de


    
    
    
    
    
