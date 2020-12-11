# IONeq

IONeq is an X-ray high-resolution photoabsorption model which compute the absorption coefficient assuming ionization equilibrium. The model takes into account turbulent broadening.  A complete description of the science behind the model is described in [Gatuzz & Churazov (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.474..696G/abstract).

OBTAINING IONeq

Different flavours of the model can be downloaded from the Github repository at https://github.com/efraingatuzz/IONeq. Currently versions include:

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


    
    
    
    
    
