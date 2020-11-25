!     com-auger.com

!     Auger ionization data structures

!     global variable declarations for Auger ionization
!     yield probabilities by Kaastra+Mewe (1993,A&ASS,97,443)

!     the dimensions of the arrays in the common blocks are 
!     defined in the file 'rates.h'

!     jAugtot   = number of tabulated ionization stages for speices k
!     nAugshtot = number of tabulated shells for species k,j
!     Xshell    = Kaastra+Mew shell index
!     W_Auger   = the Auger ejection probabilities (k,j,shell,Neject)
!                 (Neject goes up to 10)
 
      integer           jAugktot,nAugshtot,Xshell
      COMMON /iauger/   jAugktot(Imax),nAugshtot(Imax,Imax), &
                       Xshell(Imax,Imax,Nshmax)

      double precision  W_Auger
      COMMON /rauger/   W_Auger(Imax,Imax,Nshmax,10) 
                       

! eof
