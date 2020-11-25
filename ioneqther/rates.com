!     general global variables for rates.f
!     this module is included in virtually every subroutine

!     SPECIES, ABUNDANCES, IONIZATION FRACTIONS, ETC.
!     specflag = logical flag for incluse/exclusions of species k
!     Nspecies = the number of species being modeled
!     kidx     = k index of included species
!     fabund   = abundance fractions X_k/H for species k
!     xfrac    = abundance mass fractions X_k for species k
!     nk       = number density of species k (all ionization stages)
!     nkj      = number density of ionization staje j of species k 
!     eden     = electron density
!     nions    = number density of all ions
!     fion     = ionization fraction of ion k,j
!     fkjT     = stored ionization fractions at each T (for printing)
!     Temp     - temperature grid point

      logical           specflag
      COMMON/logiblk/   specflag(Imax)

      integer           Nspecies,kidx,lastcall
      COMMON/ispecblk/  Nspecies,kidx(Imax),lastcall

      double precision  nH,fabund,xfrac,nk,nkj,eden,nions,fion,fkjT,&
                        Temp,nelectron,natoms,a_ph,a_rec,a_coll
      COMMON/ions/      nH,fabund(Imax),xfrac(Imax),&
                        nk(Imax),nkj(Imax,Imax+1),eden,nions,&
                        fion(Imax,Imax+1),fkjT(Imax,Imax+1,NTmax),&
                        Temp(NTmax),nelectron(NTmax),&
                        natoms(NTmax),&
                        a_ph(Imax,Imax+1,NTmax),&
                        a_rec(Imax,Imax+1,NTmax),&
                        a_coll(Imax,Imax+1,NTmax)

!     INDIVIDUAL OR PARTIAL RATE COEFFICIENTS 
!     R_phs         = photoionization rate for shell s
!     R_ph          = photoionization rate (ejecting 1 electron)
!     alpha_rec     = recombination rate coefficient (rad,loTde,hiTde)
!     alpha_Cdi     = direct collisional ionization rate coefficient
!     alpha_Cea     = ex-autoionization rate coefficient
!     alpha_CTrecH  = charge transfer H recombination rate coefficient
!     alpha_CTionH  = charge transfer H ionization rate coefficient
!     alpha_CTrecHe = charge transfer He recombination rate coefficent
!     R_Agrout      = total Auger ionization rate out of k,j
!     Q             = Auger ionization k,j ejecting >1 electron 

      double precision  R_phs,R_ph,alpha_rec,alpha_rec_ph,alpha_rec_die       
      double precision  alpha_Cdi,alpha_Cea,alpha_CTrecH,alpha_CTionH    
      double precision  alpha_CTrecHe,R_Agrout,Q      
      COMMON/ratblock/  R_phs(NShmax),&
                        R_ph(Imax,Imax+1),alpha_rec(Imax,Imax+1),&
                        alpha_rec_ph(Imax,Imax+1),&
                        alpha_rec_die(Imax,Imax+1),&
                        alpha_Cdi(Imax,Imax+1),&
                        alpha_Cea(Imax,Imax+1),&
                        alpha_CTrecH(Imax,Imax+1),&
                        alpha_CTionH(Imax,Imax+1),&
                        alpha_CTrecHe(Imax,Imax+1),&
                        R_Agrout(Imax,Imax+1),&
                        Q(Imax,Imax+1,Imax+1)

!     RATES FOR SOLVER 
!     R_hat(k,j) = total destruction rate of k,j
!     R_ion(k,j) = total creation rate of k,j-1 by ionization of k,j
!     R_rec(k,j) = total creation rate of k,j+1 by recombination to k,j

      double precision  R_hat,R_ion,R_rec
      COMMON/totrats/   R_hat(Imax,Imax+1),&
                        R_ion(Imax,Imax+1),&
                        R_rec(Imax,Imax+1)

!     CHARACTER DATA FOR COMMUNICATIONS (file name generation, etc)
!     specID   = name of elemental species k,(i.e., k=12 is Mg)
!     ionID    = name of ion k,j (i.e., k=12,j=2 is MgII)
!     shellID  = shell spectroscopic notation (2s, etc)
!     group    = Group type on Periodic Table (i.e., IA, etc)
!     sequence = iso-electronic sequence (i.e., lithium, etc)
!     config   = iso-electronic sequence spectroscpic notation

!      character*80      specID,ionID,shellID,sequence,group,config
      character (len=80)      specID,ionID,shellID,sequence,group,config
      COMMON/chblk/     specID(Imax),ionID(Imax,Imax+1),&
                        shellID(NSHmax),sequence(Imax),&
                        group(Imax),config(Imax)

! eof
