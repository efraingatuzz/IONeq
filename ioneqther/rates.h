!     rates.h

!     master parameter file; sets physical sizes of arrays 
!     many of these parameters are employed in the common blocks

!     SPECIES/ION GRID
!............................................................................
!     Imax  = maximum number of ions
!     NSHmax = maximum number of shells (Verners)
!     STGmax = maximum number of ionization stages for recombination 
!     Smax   = maximum number of shells for direct collisional ionization 

      integer    Imax,Nmax,NSHmax,STGmax,Smax
      parameter  (Imax   =     30,&
                 Nmax   =    496,&
                 NShmax =      7,&
                 STGmax = Imax-1,&
                 Smax   =      3)



!     UVB GRID
!............................................................................
!     UVBmaxpix = number of lines in Haardt & Madau 1996 SEDs
!     UVBzpix  = number of H&M96 files (each one at a given z)

      integer    UVBmaxpix,UVBzpix
      parameter  (UVBmaxpix = 3750,&
                 UVBzpix   =   26)



!     SB99 GRID
!............................................................................
!     SB99maxpix = number of lines in Starburst99 spectra SEDs
!     SB99Nage   = number of ages of SB99 SEDS
!     SB99NZ     = number of metallicities of SB99 SEDS
!     SB99maxAZ  = maximum of SB99Nage and SB99NZ (tmp storage array)

      integer    SB99maxpix,SB99Nage,SB99NZ,SB99maxAZ
      parameter  (SB99maxpix =     1221,&
                 SB99Nage   =        5,&
                 SB99NZ     =        2,&
                 SB99maxAZ  = SB99Nage)



!     SED IONIZING SPECTRUM GRID
!............................................................................
!     SEDmaxpix = number of lines in final SEDs

      integer    SEDmaxpix
      parameter  (SEDmaxpix = UVBmaxpix)



!     ENERGY AND TEMPERATURE GRIDS
!............................................................................
!     NEmax - maximum number of energy points on E grid
!     NTmax - maximum number of temperature points of T grid

      integer    NEmax,NTmax
      parameter  (NEmax = 9000,&
                 NTmax = 1200)



!     PHYSICAL CONSTANTS AND CONVERSIONS
!............................................................................
!     h      = Planck constant  (erg sec)
!     clight = speed of light (cm/s)
!     pi     = pie, yummm
!     erg2eV = conversion of 1 eV to ergs
!     pc2cm  = conversion of 1 pc to cm

      double precision  h,clight,pi,erg2eV,pc2cm
      parameter         (h      =       6.626075540e-27,&       
                        clight =         2.99792458d10,& 
                        pi     = 3.1415926535897932384,& 
                        erg2eV =        1.60217646d-12,&
                        pc2cm  =         3.08568025d18) 
