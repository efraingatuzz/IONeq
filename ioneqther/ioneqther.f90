! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ioneqther
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
! XSPEC local model for absorption considering ionization equilibrium
! Version 1.3 January 2025
!
! Notes on version 1.3
! -- Including atom_header(30,30) in a module to be used persistently
!
! Notes on version 1.2
! -- Abundances are read from XSPEC 
! -- The model relies on the fftw library provided by XSPEC
! 
! Notes on version 1.1
! -- Now the energy grid covers until 1e6 eV
!
!
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE atom_header_mod
    INTEGER, SAVE :: atom_header(30,30)
END MODULE atom_header_mod


subroutine ioneqther(ear, ne, param, ifl, photar)
!
! The main routine to call all subroutines
!
use atom_header_mod
implicit none
integer,parameter :: num_param = 9, out_unit=20, nion=168
integer,parameter :: nemod=650000 !Number of elements for each ion cross section.
integer,parameter :: nnemod=350000
integer :: ne, ifl,size_old_cross,kk,k
double precision :: nH, rshift, emod(1:nemod), coemod(nemod),  cion(0:nion), bener(1:nemod), atomabund(1:30)
double precision :: optical_depth(nemod)
double precision :: optical_depth_convolved(0:nemod),bxs_restored(0:nion,nemod)  
double precision :: bxs_crude(0:nion,nnemod),ener_crude(0:nion,nnemod)
double precision :: zfac, Xi,ion_par, AFe_dust,vturb,temp1,Fegr_dust,tol
real :: ear(0:ne), param(num_param),photar(ne)
logical :: startup=.true.
!For the startup of the cross-section
integer :: startup_ion(0:nion)
!Variables for ionization equilibrium
include 'rates.h'
integer :: zn,ii,i
!integer :: atom_header(30,30)
integer,parameter :: znm=100
double precision :: n(znm),t,nee,lt, atom_frac(30,30)
double precision :: W_Auger_tmp(Imax,Imax,Nshmax,10)



!To read abundances from XSPEC
real :: fgabnz
external :: fgabnz 

character (len=40) version
version='1.3'
if(startup)then
 print *, ' '
 print *, 'ioneqther model Version DEV' 
 print *, 'Gatuzz & Churazov (2016)'
 print *, 'Abundaces are provided by XSPEC'  
 print *, ' '
 !To read only once the atomic and auger data
 call readAugerKM93(W_Auger_tmp)
 call read_atomic_data_header_ioneqther()
! call read_atomic_data_header_ioneqther(atom_header)
 call create_energy_grid_ioneqther(1.d1,1.d6,bener,nemod) !Absorption coefficient calculation grid  = cross section grid 
 !To read solid_iron
 call read_one_cross_sections_ioneqther(169,nnemod,bxs_crude,ener_crude,size_old_cross)
 call interpolate_cross_section_ioneqther(bxs_crude,168,ener_crude,size_old_cross,bener,bxs_restored,26,26) 
 startup=.false.   
 do i=0,nion,1
  startup_ion(i)=1
 enddo
endif 

!Tolerance to accept ionic fractions
tol=1e-3

ifl=ifl
! Model parameters
nH = param(1)
lt = param(2)
t=10**lt
ion_par=0.
!Ionization parameter is in units of keV*cm/s, in order to compare with XSTAR a conversion is required --> ion_par(ioneqther) = 8.79 + ion_par(xstar)
!For this version photoionization is not included
Xi = 10**(ion_par) 
AFe_dust = param(3)  
vturb=param(4)
rshift = param(5)
zfac = 1/(1.d0+dble(rshift))  

atomabund(2)=fgabnz(2) 
atomabund(3)=fgabnz(3) 
atomabund(4)=fgabnz(4) 
atomabund(5)=fgabnz(5) 
atomabund(6)=fgabnz(6) 
atomabund(7)=fgabnz(7)
atomabund(8)=fgabnz(8) 
atomabund(9)=fgabnz(9)
atomabund(10)=fgabnz(10) 
atomabund(11)=fgabnz(11)
atomabund(12)=fgabnz(12)
atomabund(13)=fgabnz(13)
atomabund(14)=fgabnz(14)
atomabund(15)=fgabnz(15)
atomabund(16)=fgabnz(16)
atomabund(17)=fgabnz(17)
atomabund(18)=fgabnz(18)
atomabund(19)=fgabnz(19)
atomabund(20)=fgabnz(20)
atomabund(21)=fgabnz(21)
atomabund(22)=fgabnz(22)
atomabund(23)=fgabnz(23)
atomabund(24)=fgabnz(24)
atomabund(25)=fgabnz(25)
atomabund(26)=fgabnz(26) 
atomabund(27)=fgabnz(27)
atomabund(28)=fgabnz(28)
atomabund(29)=fgabnz(29)
atomabund(30)=fgabnz(30)
Fegr_dust=AFe_dust*fgabnz(26)   

nee=1.e4  !Low electron density  
 

k=0
 
do zn=1,30,1 !Nuclear charge
! Ions with special cross-sections (He,O,Ne,Fe and also C,N,Mg,Si,S,Ar,Ca,Ni)
 if(zn.eq.1.or.zn.eq.2.or.zn.eq.6.or.zn.eq.7.or.zn.eq.8.or.zn.eq.10&
 .or.zn.eq.12.or.zn.eq.14.or.zn.eq.16.or.zn.eq.18&
 .or.zn.eq.20.or.zn.eq.26.or.zn.eq.28)then 
  call pceqsub(zn,t,nee,n,Xi,W_Auger_tmp)  
  do ii=1,zn,1 !Ion fraction except "nude core"
   if(n(ii)>tol)then !tolerance to accept ionic fraction
    atom_frac(zn,ii)=n(ii) 
    if((startup_ion(atom_header(zn,ii))).gt.0)then
     startup_ion(atom_header(zn,ii))=0
     call read_one_cross_sections_ioneqther(atom_header(zn,ii),nnemod,bxs_crude,ener_crude,size_old_cross)
     call interpolate_cross_section_ioneqther(bxs_crude,(atom_header(zn,ii)-1),ener_crude,size_old_cross,bener,bxs_restored,zn,ii)
    endif 
   else
    atom_frac(zn,ii)=0.d0
   endif

!Ion fractions for species with special cross-sections
   cion(k)=atomabund(zn)*atom_frac(zn,ii)
   k=k+1
  enddo 
 else 

!For all other elements we use Verner cross-sections
  call pceqsub(zn,t,nee,n,Xi,W_Auger_tmp)
  do ii=1,zn,1 !Ion fraction except "nude core"
   if(n(ii)>tol)then !tolerance to accept ionic fraction		
    atom_frac(zn,ii)=n(ii)
   else
    atom_frac(zn,ii)=0.d0
   endif
  enddo 

!
 endif

enddo  

call absorption_ioneqther(nH,Fegr_dust,atom_frac,atomabund,zfac, emod, nemod, optical_depth,bxs_restored,cion,bener,tol) 

if (vturb.eq.0)then 
 temp1=0.0
 optical_depth(0)=dexp(-temp1)
 do i=1,nemod
  coemod(i)=dexp(-optical_depth(i))
 enddo
else if(vturb.gt.0)then
 call optical_depth_convolved_ioneqther( nemod,optical_depth ,bener,vturb,optical_depth_convolved)
 temp1=0.0
 optical_depth_convolved(0)=dexp(-temp1)
 do i=1,nemod
  coemod(i)=dexp(-optical_depth_convolved(i))
 enddo
endif 
 
!!Move to instrument energy grid
call map_to_grid_ioneqther(dble(ear),ne,emod,nemod,photar,coemod )

return
end subroutine ioneqther
 
! ======================================= !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!TO READ ATOMIC DATA HEADER!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_atomic_data_header_ioneqther()
!subroutine read_atomic_data_header_ioneqther(atom_header)
!
! This routine reads the atomic data IDs from first header
! on atomic data file
!
!
use atom_header_mod
implicit none
integer,parameter :: nion=168, out_unit=20
integer ::   i,  status
double precision :: z(nion), charge(nion), column_id(nion)
!integer :: atom_header(30,30)
character (*), parameter :: fileloc = '/atomic_data/AtomicData.fits'
character (*), parameter :: ismreadchat = 'ioneqther: reading from '
character (len=255 + 29) :: filename2 ! len(fileloc)

character (len=240) :: dir_input,tmp_txt1 
character (len=255) :: ioneqther_root = ''
character (len=len(ismreadchat)+len(filename2)) :: chatmsg = ''
integer inunit,readwrite,blocksize
integer :: hdutype
integer :: felem=1, nulld=0
logical :: anynull
character (len=255) :: fgmstr
external :: fgmstr
character (len=240) :: local_dir = '/media/efrain/DATA/softwares/modelosXSPEC/IonEq/thermal/solar/v1.1'  
 
! Where do we look for the data?
ioneqther_root = trim(fgmstr('ioneqtherROOT'))
if (ioneqther_root .EQ. '') then
ioneqther_root = local_dir
endif

! parameters to specify the opening process
status=0
readwrite=0
blocksize=1
filename2=trim(ioneqther_root) // fileloc
chatmsg=ismreadchat // filename2
call xwrite(chatmsg, 20)
! Get an unused Logical Unit Number to use to open the FITS file.
call ftgiou(inunit,status)
! Open the FITS file
call ftopen(inunit,filename2,readwrite,blocksize,status)
! Move to the extension 3 (the binary table)
call ftmahd(inunit,3,hdutype,status)

do i=1, nion
call ftgcvd(inunit,2,i,felem,1,nulld,z(i),anynull,status)
call ftgcvd(inunit,3,i,felem,1,nulld,charge(i),anynull,status)
call ftgcvd(inunit,4,i,felem,1,nulld,column_id(i),anynull,status)
z(i)=int(z(i))
charge(i)=int(charge(i)) 
column_id(i)=int(column_id(i))
atom_header(int(z(i)),int(charge(i)))=int(column_id(i))
enddo
 
! Report on errors (done before closing the file in case the error
! comes from closing the file). Unfortunately the X-Spec API does not
! provide a way to signal an error to the calling code, so a screen
! message is used, using the same method used to report the model
! the first time it is used. An alternative would be to use xwrite()
! with a low chatter level.
!
! This message could be displayed only once, but it is probaly worth
! repeating each time it is used.
if (status .ne. 0) then
write (*,*) 'ERROR: unable to read cross sections from ', filename2
endif
! Close the file and free the unit number
call ftclos(inunit, status)
call ftfiou(-1, status)
end subroutine read_atomic_data_header_ioneqther
! ======================================= !


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!TO READ ONLY ONE CROSS_SECTION!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine read_one_cross_sections_ioneqther(column_number,bnene,xs,ener,nelemm)
!
! This routine reads all cross sections and puts them on a given grid
!
! It uses the X-Spec local variable/dictionary value ISMABSROOT
! to locate the data file. If this is not set then it uses
! the setting of the local_dir parameter below (which should be
! changed to match the location of the data file). The path -
! however it is given - should point to the location that contains
! the atomic_data/ directory - i.e. it loads
! <path>/atomic_data/AtomicData.fits
!
implicit none
integer,parameter :: nion=168, out_unit=20
integer :: bnene,  i, j, status,column_number
integer :: nemax(0:nion)  
double precision :: ener(0:nion,bnene), xs(0:nion,bnene) 
 
character (*), parameter :: fileloc = '/atomic_data/AtomicData.fits'
character (*), parameter :: ismreadchat = 'ioneqther: reading from '
character (len=255 + 29) :: filename2 ! len(fileloc)
character (len=240) :: dir_input,tmp_txt1 
 
character (len=255) :: ioneqther_root = ''
character (len=len(ismreadchat)+len(filename2)) :: chatmsg = ''
integer inunit,readwrite,blocksize,nelemm,offset
integer :: hdutype
integer :: nulld=0, logical_start(0:nion)
logical :: anynull
character (len=255) :: fgmstr
external :: fgmstr
character (len=240) :: local_dir = '/media/efrain/DATA/softwares/modelosXSPEC/IonEq/thermal/solar/v1.1'
 

!Number of elements for each ion cross section.
do i=0,nion
nemax(i)=650000
enddo
! Where do we look for the data?
ioneqther_root = trim(fgmstr('ioneqtherROOT'))
if (ioneqther_root .EQ. '') then
ioneqther_root = local_dir
endif
! parameters to specify the opening process
status=0
readwrite=0
blocksize=1
filename2=trim(ioneqther_root) // fileloc
chatmsg=ismreadchat // filename2
call xwrite(chatmsg, 20)
! Get an unused Logical Unit Number to use to open the FITS file.
call ftgiou(inunit,status)
! Open the FITS file
call ftopen(inunit,filename2,readwrite,blocksize,status)
! Move to the extension 2 (the binary table)
call ftmahd(inunit,2,hdutype,status)
!Read length of energy array
call ftgdes(inunit,column_number,1,nelemm,offset,status)
i=column_number-1

do j=1,nelemm

!Read cross-section 
call ftgcvd(inunit,column_number,1,j,1,nulld,ener(i,j),anynull,status)
 
!Read energy
call ftgcvd(inunit,column_number,2,j,1,nulld,xs(i,j),anynull,status)
  
enddo
logical_start(i)=1


! Report on errors (done before closing the file in case the error
! comes from closing the file). Unfortunately the X-Spec API does not
! provide a way to signal an error to the calling code, so a screen
! message is used, using the same method used to report the model
! the first time it is used. An alternative would be to use xwrite()
! with a low chatter level.
!
! This message could be displayed only once, but it is probaly worth
! repeating each time it is used.
if (status .ne. 0) then
write (*,*) 'ERROR: unable to read cross sections from ', filename2
endif
! Close the file and free the unit number
call ftclos(inunit, status)
call ftfiou(-1, status)
end subroutine read_one_cross_sections_ioneqther

 

subroutine absorption_ioneqther(col22, Fegr_dust,atom_frac,atomabund, &
zfac, e1, bnene, coeff, bxs2,cion, bener,tol)
!
! This is routine that calculates the optical depth given the column densities
! Finally returns the absorption coefficient exp(-tau)
!
implicit none
include           'rates.h'
integer,parameter :: nion=168, out_unit=20
integer :: bnene 
integer :: i,j,iz,in,is,Z,k
double precision :: col22, col, tmp, cion(0:nion),vernercross, cionverner,tol 
double precision :: bener(0:bnene), bxs2(0:nion,bnene), e1(0:bnene), atomabund(1:30)
double precision, parameter :: c=299792.458 
double precision :: tau, coeff(bnene),Fegr_dust  
double precision :: zfac,atom_frac(30,30)  
!To compute cross-sections from Verner
real :: thresh,S,s_tmp 

!To read abundances from XSPEC
real :: fgabnz
external :: fgabnz  
 
! Calculates the optical depth and the absorption coefficient exp(-tau)
col=col22*1.d22
e1(0)=(bener(0)*zfac)/1.d3
coeff(0)=0.0  
 
!To compute optical depth per energies
do i=1,bnene
 e1(i)=(bener(i)*zfac)/1.d3 
!Hydrogen 
 tmp=bxs2(0,i)*1.0*atom_frac(1,1)*(col)  
! Calculates the optical depth for metals with special cross-sections  (including solid_iron)
 do j=1,nion 
  tmp=tmp+(col)*(cion(j)*bxs2(j,i))
 enddo
!To add solid_iron 
 tmp=tmp+(col)*(Fegr_dust*bxs2(168,i)) 

!!!!Calculates the optical depth for elements for which I am using Verner cross-section!!!
 do Z=1,30
  if(Z.eq.3.or.Z.eq.4.or.Z.eq.5.or.Z.eq.9.or.Z.eq.11.or.Z.eq.13&
  .or.Z.eq.15.or.Z.eq.17.or.Z.eq.19.or.Z.eq.21.or.Z.eq.22.or.Z.eq.23&
  .or.Z.eq.24.or.Z.eq.25.or.Z.eq.27.or.Z.eq.29.or.Z.eq.30)then
   vernercross=0.
   do j=1,Z
    if(atom_frac(Z,j).gt.tol)then
     cionverner=atomabund(Z)*atom_frac(Z,j) 
     do is=1,7
      call phfit2_b(Z,j,is,real(e1(i)),s_tmp,thresh)   
      vernercross=vernercross+s_tmp*1e-18
     enddo
     tmp=tmp+(col)*(cionverner*vernercross)
    endif
   enddo
  endif 
 enddo 
!!!!!

! Calculates the optical depth
 tau=tmp
 coeff(i)=tau
 
enddo 


end subroutine absorption_ioneqther
! ======================================= !



subroutine map_to_grid_ioneqther(new_en,nne,old_en, one, nflux, old_flu )
! This routine map to a given grid
implicit none
integer :: i, j, k, one, nne, bmin, bmax 
double precision :: new_en(0:nne)
double precision :: old_en(0:one), old_flu(one)
double precision :: stemp,etemp, s, etemp2
real :: nflux(nne)
integer,parameter :: out_unit=20
do i=1,nne
nflux(i)=real(0.d0)
call dbinsrch_ioneqther(new_en(i-1),bmin,old_en,one+1)
call dbinsrch_ioneqther(new_en(i),bmax,old_en,one+1)
bmin = bmin-1
bmax = bmax-1
! Linear interpolation
if (bmin.eq.bmax) then
if(new_en(i).le.old_en(1))then
s=real(old_flu(1))
else if(new_en(i).gt.old_en(one))then
s=real(old_flu(one))
else
do j=2,one
if(new_en(i).gt.old_en(j-1).and.new_en(i).le.old_en(j))then
etemp2=(new_en(i)+new_en(i-1))/2
s=old_flu(j-1)+(old_flu(j)-old_flu(j-1))*(etemp2-old_en(j-1))/(old_en(j)-old_en(j-1))
endif
enddo
endif
! Average (integral)
else
stemp=0.d0
etemp=0.d0
do k=bmin,bmax
stemp=stemp+(old_flu(k))*(old_en(k)-old_en(k-1))
etemp=etemp+(old_en(k)-old_en(k-1))
enddo
s=real(stemp/etemp)
endif
nflux(i)=real(s)
enddo
end subroutine map_to_grid_ioneqther
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dbinsrch_ioneqther(e,k,ener,n)
!
! search for energy e in array ener(1:n) and return bin k for
! which ener(k).le.e.lt.ener(k+1)
! adapted from J. Wilms's tbvabs_new.f routine
!
implicit none
double precision :: e
integer :: n,k,klo,khi
double precision :: ener(n)
klo=1
khi=n
1 if (khi-klo.gt.1) then
k=(khi+klo)/2
if(ener(k).gt.e)then
khi=k
else
klo=k
endif
goto 1
endif
if (klo.eq.0) then
print *,'Energy out of bounds. Should not happen'
stop
endif
k=klo
end subroutine dbinsrch_ioneqther
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine create_energy_grid_ioneqther(emin,emax,en,nen)
implicit none
integer :: i, nen
double precision :: en(nen)
double precision :: emin,emax
!
do i=1,nen
en(i)=10**(((log10(emax)-log10(emin))*real(i)/nen)+log10(emin))
enddo
!
end subroutine create_energy_grid_ioneqther



!====================================================================
      function background_ioneqther(e_bkg,e_min)
!--------------------------------------------------------------------
!     Intensity of the ionizing background
!     in photons/cm2/s/keV, e - energy (keV)
!--------------------------------------------------------------------
      implicit none  
      double precision, parameter :: pi=3.14159265358979312d0
      real (kind = 8) :: background_ioneqther
      real (kind = 8) :: norm,alpha,pgamma
      real (kind = 4) :: e_bkg,e_min
      real (kind = 8) :: i0u,bu,xrb
      background_ioneqther=0

!!! Powerlaw spectral shape with gamma = 2.1 (i.e. as warmabs per defect)
      pgamma=2.1
      alpha=pgamma-1
      norm=(alpha-1)*e_min**(alpha-1)
      background_ioneqther=norm*(e_bkg**(-alpha))


      return
      end


!====================================================================

!=================================================================
      subroutine pceqsub(zn,t0,nee,n,Xi,W_Auger_tmp)
!-----------------------------------------------------------------
!     Calculate ionization equilibrium
!     Proc: photo,col_ion,rad_rec,diel_rec
!     Atomic data from Verner's database
!
!     zn - nuclear charge
!     t0 - temperature (K)
!     nee - electron concentration (cm-3)
!     n(1,...,zn+1) - fraction of given (SS) ion 
!     photoionization background is define by a function 
!     background(e) in photons/cm2/s/keV, e - energy (keV)
!
!
!     Notes:
!      1. It is likely that below T~10^4 K the dielectronic recombination 
!         data are not reliable.
!
!     Chur, 12/11/99
!-----------------------------------------------------------------
      implicit none 
      include           'rates.h'
!      double precision, parameter :: pi=3.14159265358979312d0
      integer :: iz,in,i,is,z,zn,iii,m,mmm,j
      integer, parameter :: zm=100
      double precision :: n(*),t0,nee
      double precision :: Xi,tmp2,Sk,tmp1
      real :: t,r,d,c,p,e_pceq,s_tmp,thresh,e1,e2,de,p1,p2,r2,temp,xi_factor
      double precision :: wa,R_auger
      double precision :: background_ioneqther,integration_background_ioneqther
      double precision :: phi(zm,zm),pp(zm),rr(zm),dd(zm),cc(zm),rr2(zm,zm,zm)
      double precision :: R_ion(zm,zm),R_rec(zm,zm),R_hat(zm,zm),fion(zm,zm),Pk(zm,zm)
      double precision :: W_Auger_tmp(Imax,Imax,Nshmax,10)
      integer, dimension(1:30) :: max_shell = (/1,1,2,2,3,3,3,3,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7/)
 
      nee=nee

      thresh=0.0
      iz=zn

!--------------------------------------------------------------------
!     This functions return the value:
!     sqrt(xi/(integral(J(E)dE/ne)))
!     which multiplies the ionization rate (substitutes the 4*pi factor
!     in equation 2 of Gatuzz & Churazov (2017)
!--------------------------------------------------------------------
            e1=0.01360000
            e2=13.6056980
            de=e1*0.01
            e_pceq=e1
            temp=0.
            do while (e_pceq < (e2+de))
              temp=temp+background_ioneqther(e_pceq,e1)*(de)/nee
              e_pceq = e_pceq+de
            end do 
            xi_factor=sqrt(Xi/temp)

!------------------------------
      do i=1,zn+1
         n(i)=0
      end do 

      t=t0
      do in=1,iz



         i=iz-in+1
         r=0
         d=0
         c=0
         p=0
         p1=0
         p2=0
         
!------- radiative recombination
         call rrfit_b(iz,in,t,r)
         rr(i)=r

!------ dielectronic recombination
        call drfit2_a(iz,in,t,d)
        dd(i)=d
 


!------- collisional ionization
         call  cfit_b(iz,in,t,c) 
         cc(i)=c

!------- photoionization
         do is=1,7


            e_pceq=1000.
            call phfit2_b(iz,in,is,e_pceq,s_tmp,thresh)

            e1=thresh/1.d3
            e2=8.*e1
            de=e1*0.01
            wa=0
            

            e_pceq=e1
            do while (e_pceq < (e2+de))

               call phfit2_b(iz,in,is,e_pceq*1000.,s_tmp,thresh)

!Without ionization parameter
               p=p+4.d0*pi*s_tmp*1d-18*background_ioneqther(e_pceq,e1)*(de/(e_pceq)) 
               wa=wa+4.d0*pi*s_tmp*1d-18*background_ioneqther(e_pceq,e1)*(de/(e_pceq)) 

 
               e_pceq = e_pceq+de
            end do 

!Multiply for Auger probability
      mmm = max_shell(i)
        do m=in+2,iz-1
             iii    = m - in
               if(is.le.mmm.or.is.eq.mmm)then
                       if(iii.le.10)then  ! tables stop at 10 ejected electrons
 
                              p2=p2+p*W_Auger_tmp(iz,in,is,iii)
                       endif
               endif
        enddo
 
!Recombination probabilities from Auger effect (Upper levels that recombine)
       do m=1,in-2
             iii    = m 
               if(is.le.mmm.or.is.eq.mmm)then
                   if(iii.le.10)then  ! tables stop at 10 ejected electrons
 
                       r2=r2+p*W_Auger_tmp(iz,in,is,iii)
                       rr2(iz,iii,is)=r2
!                       rr2(iz,iii,is)=0
                    endif
                endif
           enddo
 
         end do 
         p=p+p2 !Auger effect
         p=(Xi*p)/(4*pi) !Ionization parameter 
         pp(i)=p
 
      end do 



!!!ION FRACTIONS !!!
   
!!! FOR HYDROGEN !!!!!
    if(iz.eq.1)then
    phi(1,1)=(pp(1)+cc(1))/(rr(1))
    n(1)=1.d0/(1.d0+phi(1,1))
    n(2)=phi(1,1)*n(1)
    endif

!!! FOR HELIUM !!!!!
    if(iz.eq.2)then
    phi(2,1)=(pp(1)+cc(1))/rr(1)
    phi(2,2)=(pp(2)+cc(2))/rr(2) 
    n(1)=1.d0/(1.d0+phi(2,1)+(phi(2,1)*phi(2,2)))
    n(2)=phi(2,1)*n(1)
    n(3)=phi(2,2)*n(2)
    endif

!!! For Metals!!!
    if(iz.eq.3.or.iz.gt.3)then

!     rate matrix notation: 
!     R_ion(iz,j) = destruction of iz,j to iz,j+1 via ionization
!     R_rec(iz,j) = creation of iz,j via recombination to iz,j-1
!     R_hat(iz,j) = total destruction of iz,j by all processes

!     .................................
!     Population of each matrix
!     .................................

      do z=3,iz
!H
    phi(1,1)=(pp(1)+cc(1))/(rr(1))
    fion(1,1)=1.d0/(1.d0+phi(1,1))
    fion(1,2)=phi(1,1)*n(1)

!He
    phi(2,1)=(pp(1)+cc(1))/rr(1)
    phi(2,2)=(pp(2)+cc(2))/rr(2) 
    fion(2,1)=1.d0/(1.d0+phi(2,1)+(phi(2,1)*phi(2,2)))
    fion(2,2)=phi(2,1)*n(1)
    fion(2,3)=phi(2,2)*n(2)

!Metals
!     j=1 (neutral ion) easy case

       j = 1
       R_ion(z,j)=(pp(j)+cc(j))
       R_rec(z,j) = 0.0d0  ! never accessed   
       R_hat(z,j) = R_ion(z,j)         

       

!     all cases from j=2 to j=z

       do j=2,z
        R_ion(z,j)=pp(j)+(cc(j))
        R_rec(z,j) = rr(j-1) +dd(j-1)
        R_hat(z,j) = R_ion(z,j) + R_rec(z,j)   
       enddo

!     j=k+1 (fully ionized ion) is an special case

       j = z + 1
       R_ion(z,j) = 0.0d0  ! never accessed
       R_rec(z,j) = rr(j-1)+dd(j-1) 
       R_hat(z,j) = 0.0d0  ! never accessed


      enddo  ! next z


!     ......................................................
!     STEP 2: solve the equations
!     ......................................................


      do z=3,iz

!     j=1 and j=2: trivial cases, solve directly

       phi(z,1) = R_hat(z,1)/R_rec(z,2)
       phi(z,2) = (R_hat(z,2)-(R_ion(z,1)/phi(z,1)))/R_rec(z,3)

!     j=3 to j=k  

       do j=3,z       
        R_auger = 0.0d0
        tmp1 = phi(z,j-1)
        do i=j-2,1,-1    
         tmp1   = phi(z,i)*tmp1
         R_auger=R_auger + rr2(z,j,i)/abs(tmp1)
        enddo
        tmp2       = (R_ion(z,j-1)/phi(z,j-1))/R_hat(z,j)+R_auger/R_hat(z,j)
        phi(z,j) = (R_hat(z,j)/R_rec(z,j+1))*tmp2*(1.0d0/tmp2-1.0d0)
       enddo

!     denominator of the ionization fractions

       Sk    = 0.0d0
       Pk(z,1) = 1.0d0
       do j=2,z+1
        Pk(z,j) = Pk(z,j-1)*abs(phi(z,j-1))
       enddo      
       do j=1,z+1
        Sk = Sk + Pk(z,j)
       enddo

!     compute the ionization fractions

       n(1) = 1.0d0/Sk
       do j=2,z+1
        n(j) = abs(phi(z,j-1))*n(j-1)
       enddo


      enddo  ! next k


    endif



      return
      end subroutine pceqsub





!=====================================================================
      subroutine drfit2_a(iz,in,t,d)
!---------------------------------------------------------------------
!     dielectronic recombination rate for T>6.10^4 K
!---------------------------------------------------------------------
      implicit none
      integer :: iz,in
      real :: t,d
      integer :: n,i
      integer, parameter :: nm=10000
      integer :: iza(nm),ina(nm)
      real :: aa(nm),ba(nm),t0a(nm),t1a(nm)
      integer,parameter :: ndr=167
      real :: drec_a(6,ndr)
      common/drec_a/drec_a

!      logical start/.true./
      logical start
      data start/.true./
      save start
      
      if(start) then
         n=ndr
         do i=1,n
            iza(i)=nint(drec_a(1,i))
            ina(i)=nint(drec_a(2,i))
            aa(i)=drec_a(3,i)
            ba(i)=drec_a(4,i)
            t0a(i)=drec_a(5,i)
            t1a(i)=drec_a(6,i)
!>            print 10,iza(i),ina(i),aa(i),ba(i),t0a(i),t1a(i)
!> 10         format(2i4,1p4e9.2)
         end do 
         start=.false.
      endif
      
      d=0
      do i=1,n
         if(iza(i).eq.iz.and.ina(i).eq.in) goto 200
      end do 
      return

 200  continue
      d=aa(i)*t**(-3./2.)*exp(-t0a(i)/t)*(1.+ba(i)*exp(-t1a(i)/t))

      return
      end subroutine drfit2_a
 !==========================================================================
      block data derec_a
!--------------------------------------------------------------------------
!     dielectronic recombination data (T>6*10^4 K)
!--------------------------------------------------------------------------
      integer, parameter :: ndr=167
      real :: drec_a(6,ndr)
      common/drec_a/drec_a
!      DATA(CF(I, 4, 3),I=1,5)/   18.2,1.,2.08E-08,0.4390,0.21/
!     data ((drec(i,j),i=1,6),j=1,100)
      data (drec_a(i,1),i=1,6) /2,  2, 1.90E-03, 3.00E-01, 4.70E+05, 9.40E+04/
      data (drec_a(i,2),i=1,6) /6,  6, 2.54E-03, 4.42E-02, 1.57E+05, 3.74E+05/
      data (drec_a(i,3),i=1,6) /6,  5, 6.15E-03, 5.88E-02, 1.41E+05, 1.41E+05/
      data (drec_a(i,4),i=1,6) /6,  4, 1.62E-03, 3.43E-01, 8.19E+04, 1.59E+05/
      data (drec_a(i,5),i=1,6) /6,  3, 4.78E-02, 3.62E-01, 3.44E+06, 5.87E+05/
      data (drec_a(i,6),i=1,6) /6,  2, 3.22E-02, 3.15E-01, 4.06E+06, 8.31E+05/
      data (drec_a(i,7),i=1,6) /7,  7, 2.98E-03, 0.00E+00, 2.20E+05, 1.00E+05/
      data (drec_a(i,8),i=1,6) /7,  6, 7.41E-03, 7.64E-02, 2.01E+05, 7.37E+04/
      data (drec_a(i,9),i=1,6) /7,  5, 1.13E-02, 1.64E-01, 1.72E+05, 2.25E+05/
      data (drec_a(i,10),i=1,6) /7,  4, 2.62E-03, 2.43E-01, 1.02E+05, 1.25E+05/
      data (drec_a(i,11),i=1,6) /7,  3, 7.50E-02, 3.50E-01, 4.75E+06, 8.35E+05/
      data (drec_a(i,12),i=1,6) /7,  2, 4.61E-02, 3.09E-01, 5.44E+06, 1.14E+06/
      data (drec_a(i,13),i=1,6) /8,  8, 1.11E-03, 9.25E-02, 1.75E+05, 1.45E+05/
      data (drec_a(i,14),i=1,6) /8,  7, 5.07E-03, 1.81E-01, 1.98E+05, 3.35E+05/
      data (drec_a(i,15),i=1,6) /8,  6, 1.48E-02, 3.05E-01, 2.41E+05, 2.83E+05/
      data (drec_a(i,16),i=1,6) /8,  5, 1.84E-02, 1.00E-01, 2.12E+05, 2.83E+05/
      data (drec_a(i,17),i=1,6) /8,  4, 4.13E-03, 1.62E-01, 1.25E+05, 2.27E+05/
      data (drec_a(i,18),i=1,6) /8,  3, 1.06E-01, 3.40E-01, 6.25E+06, 1.12E+06/
      data (drec_a(i,19),i=1,6) /8,  2, 4.72E-02, 0.00E+00, 6.00E+06, 1.00E+05/
      data (drec_a(i,20),i=1,6) /10, 10, 9.77E-04, 7.30E-02, 3.11E+05, 2.06E+05/
      data (drec_a(i,21),i=1,6) /10,  9, 2.65E-03, 2.42E-01, 2.84E+05, 3.07E+05/
      data (drec_a(i,22),i=1,6) /10,  8, 3.69E-03, 1.01E+00, 2.24E+05, 2.94E+05/
      data (drec_a(i,23),i=1,6) /10,  7, 1.18E-02, 3.91E-01, 2.70E+05, 5.50E+05/
      data (drec_a(i,24),i=1,6) /10,  6, 2.44E-02, 2.52E+00, 3.09E+05, 9.91E+05/
      data (drec_a(i,25),i=1,6) /10,  5, 3.02E-02, 4.45E-01, 2.83E+05, 1.73E+06/
      data (drec_a(i,26),i=1,6) /10,  4, 6.10E-03, 2.54E-01, 1.68E+05, 6.13E+05/
      data (drec_a(i,27),i=1,6) /10,  3, 2.52E-01, 3.04E-01, 1.40E+07, 1.80E+06/
      data (drec_a(i,28),i=1,6) /10,  2, 7.14E-02, 2.96E-01, 1.10E+07, 2.24E+06/
      data (drec_a(i,29),i=1,6) /12, 12, 4.49E-04, 2.10E-02, 5.01E+04, 2.81E+04/
      data (drec_a(i,30),i=1,6) /12, 11, 1.95E-03, 7.40E-02, 6.06E+05, 1.44E+06/
      data (drec_a(i,31),i=1,6) /12, 10, 5.12E-03, 3.23E-01, 4.69E+05, 7.55E+05/
      data (drec_a(i,32),i=1,6) /12,  9, 7.74E-03, 6.36E-01, 3.74E+05, 7.88E+05/
      data (drec_a(i,33),i=1,6) /12,  8, 1.17E-02, 8.07E-01, 3.28E+05, 1.02E+06/
      data (drec_a(i,34),i=1,6) /12,  7, 3.69E-02, 3.51E-01, 4.80E+05, 9.73E+05/
      data (drec_a(i,35),i=1,6) /12,  6, 3.63E-02, 5.48E-01, 3.88E+05, 7.36E+05/
      data (drec_a(i,36),i=1,6) /12,  5, 4.15E-02, 2.33E-01, 3.39E+05, 3.82E+05/
      data (drec_a(i,37),i=1,6) /12,  4, 8.86E-03, 3.18E-01, 2.11E+05, 1.54E+06/
      data (drec_a(i,38),i=1,6) /12,  3, 2.52E-01, 3.15E-01, 1.40E+07, 2.64E+06/
      data (drec_a(i,39),i=1,6) /12,  2, 9.28E-02, 0.00E+00, 1.45E+07, 1.00E+05/
      data (drec_a(i,40),i=1,6) /13, 13, 2.00E-03, 1.34E+00, 4.93E+04, 1.01E+05/
      data (drec_a(i,41),i=1,6) /13, 12, 2.35E-03, 1.25E-01, 6.72E+04, 2.45E+04/
      data (drec_a(i,42),i=1,6) /13, 11, 2.97E-03, 3.49E-01, 7.90E+05, 9.20E+05/
      data (drec_a(i,43),i=1,6) /13, 10, 6.97E-03, 1.12E-01, 8.42E+05, 3.34E+05/
      data (drec_a(i,44),i=1,6) /13,  9, 1.22E-02, 1.29E+00, 4.29E+05, 8.42E+05/
      data (drec_a(i,45),i=1,6) /13,  8, 1.84E-02, 1.36E+00, 3.72E+05, 1.45E+06/
      data (drec_a(i,46),i=1,6) /13,  7, 3.32E-02, 1.24E+00, 4.17E+05, 1.55E+06/
      data (drec_a(i,47),i=1,6) /13,  6, 3.78E-02, 8.02E-01, 3.66E+05, 1.41E+06/
      data (drec_a(i,48),i=1,6) /13,  5, 5.12E-02, 1.72E-01, 3.59E+05, 3.39E+05/
      data (drec_a(i,49),i=1,6) /13,  4, 1.15E-02, 9.43E-01, 2.31E+05, 2.97E+06/
      data (drec_a(i,50),i=1,6) /13,  3, 2.85E-01, 3.12E-01, 1.61E+07, 3.11E+06/
      data (drec_a(i,51),i=1,6) /13,  2, 1.02E-01, 1.25E-01, 1.71E+07, 1.91E+06/
      data (drec_a(i,52),i=1,6) /14, 14, 1.10E-03, 0.00E+00, 7.70E+04, 1.00E+05/
      data (drec_a(i,53),i=1,6) /14, 13, 5.87E-03, 7.53E-01, 9.63E+04, 6.46E+04/
      data (drec_a(i,54),i=1,6) /14, 12, 5.03E-03, 1.88E-01, 8.75E+04, 4.71E+04/
      data (drec_a(i,55),i=1,6) /14, 11, 5.43E-03, 4.50E-01, 1.05E+06, 7.98E+05/
      data (drec_a(i,56),i=1,6) /14, 10, 8.86E-03, 0.00E+00, 1.14E+06, 1.00E+05/
      data (drec_a(i,57),i=1,6) /14,  9, 1.68E-02, 1.80E+00, 4.85E+05, 1.03E+06/
      data (drec_a(i,58),i=1,6) /14,  8, 2.49E-02, 1.88E+00, 4.15E+05, 1.91E+06/
      data (drec_a(i,59),i=1,6) /14,  7, 3.13E-02, 2.01E+00, 3.66E+05, 2.11E+06/
      data (drec_a(i,60),i=1,6) /14,  6, 4.25E-02, 1.22E+00, 3.63E+05, 2.14E+06/
      data (drec_a(i,61),i=1,6) /14,  5, 6.18E-02, 3.03E-01, 3.88E+05, 1.12E+06/
      data (drec_a(i,62),i=1,6) /14,  4, 1.38E-02, 1.42E+00, 2.51E+05, 3.93E+06/
      data (drec_a(i,63),i=1,6) /14,  3, 3.27E-01, 3.06E-01, 1.88E+07, 3.60E+06/
      data (drec_a(i,64),i=1,6) /14,  2, 1.13E-01, 2.86E-01, 1.99E+07, 4.14E+06/
      data (drec_a(i,65),i=1,6) /16, 16, 1.62E-03, 0.00E+00, 1.25E+05, 1.00E+05/
      data (drec_a(i,66),i=1,6) /16, 15, 1.09E-02, 1.20E-02, 1.92E+05, 1.80E+04/
      data (drec_a(i,67),i=1,6) /16, 14, 3.35E-02, 6.59E-02, 1.89E+05, 1.59E+05/
      data (drec_a(i,68),i=1,6) /16, 13, 3.14E-02, 6.89E-02, 1.68E+05, 8.04E+04/
      data (drec_a(i,69),i=1,6) /16, 12, 1.27E-02, 1.87E-01, 1.38E+05, 1.71E+05/
      data (drec_a(i,70),i=1,6) /16, 11, 1.47E-02, 1.29E-01, 1.80E+06, 1.75E+06/
      data (drec_a(i,71),i=1,6) /16, 10, 1.34E-02, 1.04E+00, 6.90E+05, 2.15E+06/
      data (drec_a(i,72),i=1,6) /16,  9, 2.38E-02, 1.12E+00, 5.84E+05, 2.59E+06/
      data (drec_a(i,73),i=1,6) /16,  8, 3.19E-02, 1.40E+00, 5.17E+05, 2.91E+06/
      data (drec_a(i,74),i=1,6) /16,  7, 7.13E-02, 1.00E+00, 6.66E+05, 2.32E+06/
      data (drec_a(i,75),i=1,6) /16,  6, 8.00E-02, 5.55E-01, 6.00E+05, 2.41E+06/
      data (drec_a(i,76),i=1,6) /16,  5, 7.96E-02, 1.63E+00, 5.09E+05, 6.37E+06/
      data (drec_a(i,77),i=1,6) /16,  4, 1.34E-02, 3.04E-01, 2.91E+05, 1.04E+06/
      data (drec_a(i,78),i=1,6) /16,  3, 4.02E-01, 2.98E-01, 2.41E+07, 4.67E+06/
      data (drec_a(i,79),i=1,6) /16,  2, 1.45E-01, 2.81E-01, 2.54E+07, 5.30E+06/
      data (drec_a(i,80),i=1,6) /18, 18, 1.00E-03, 5.00E-03, 3.20E+05, 3.10E+05/
      data (drec_a(i,81),i=1,6) /18, 17, 1.10E-02, 4.50E-02, 2.90E+05, 5.50E+05/
      data (drec_a(i,82),i=1,6) /18, 16, 3.40E-02, 5.70E-02, 2.39E+05, 6.00E+05/
      data (drec_a(i,83),i=1,6) /18, 15, 6.85E-02, 8.70E-02, 2.56E+05, 3.81E+05/
      data (drec_a(i,84),i=1,6) /18, 14, 9.00E-02, 7.69E-02, 2.50E+05, 3.30E+05/
      data (drec_a(i,85),i=1,6) /18, 13, 6.35E-02, 1.40E-01, 2.10E+05, 2.15E+05/
      data (drec_a(i,86),i=1,6) /18, 12, 2.60E-02, 1.20E-01, 1.80E+05, 2.15E+05/
      data (drec_a(i,87),i=1,6) /18, 11, 1.70E-02, 1.00E-01, 2.70E+06, 3.30E+06/
      data (drec_a(i,88),i=1,6) /18, 10, 2.10E-02, 1.92E+00, 8.30E+05, 3.50E+06/
      data (drec_a(i,89),i=1,6) /18,  9, 3.50E-02, 1.66E+00, 6.95E+05, 3.60E+06/
      data (drec_a(i,90),i=1,6) /18,  8, 4.30E-02, 1.67E+00, 6.05E+05, 3.80E+06/
      data (drec_a(i,91),i=1,6) /18,  7, 7.13E-02, 1.40E+00, 6.68E+05, 2.90E+06/
      data (drec_a(i,92),i=1,6) /18,  6, 9.60E-02, 1.31E+00, 6.50E+05, 3.60E+06/
      data (drec_a(i,93),i=1,6) /18,  5, 8.50E-02, 1.02E+00, 5.30E+05, 2.80E+06/
      data (drec_a(i,94),i=1,6) /18,  4, 1.70E-02, 2.45E-01, 3.55E+05, 1.10E+06/
      data (drec_a(i,95),i=1,6) /18,  3, 4.76E-01, 2.94E-01, 3.01E+07, 6.05E+06/
      data (drec_a(i,96),i=1,6) /18,  2, 2.97E-01, 2.77E-01, 3.13E+07, 6.54E+06/
      data (drec_a(i,97),i=1,6) /20, 20, 3.28E-04, 9.07E-02, 3.46E+04, 1.64E+04/
      data (drec_a(i,98),i=1,6) /20, 19, 5.84E-02, 1.10E-01, 3.85E+05, 2.45E+05/
      data (drec_a(i,99),i=1,6) /20, 18, 1.12E-01, 1.74E-02, 4.08E+05, 4.27E+05/
      data (drec_a(i,100),i=1,6) /20, 17, 1.32E-01, 1.32E-01, 3.82E+05, 6.92E+05/
      data (drec_a(i,101),i=1,6) /20, 16, 1.33E-01, 1.14E-01, 3.53E+05, 8.78E+05/
      data (drec_a(i,102),i=1,6) /20, 15, 1.26E-01, 1.62E-01, 3.19E+05, 7.43E+05/
      data (drec_a(i,103),i=1,6) /20, 14, 1.39E-01, 8.78E-02, 3.22E+05, 6.99E+05/
      data (drec_a(i,104),i=1,6) /20, 13, 9.55E-02, 2.63E-01, 2.47E+05, 4.43E+05/
      data (drec_a(i,105),i=1,6) /20, 12, 4.02E-02, 6.27E-02, 2.29E+05, 2.81E+05/
      data (drec_a(i,106),i=1,6) /20, 11, 4.19E-02, 6.16E-02, 3.73E+06, 5.84E+06/
      data (drec_a(i,107),i=1,6) /20, 10, 2.57E-02, 2.77E+00, 9.26E+05, 4.89E+06/
      data (drec_a(i,108),i=1,6) /20,  9, 4.45E-02, 2.23E+00, 7.96E+05, 4.62E+06/
      data (drec_a(i,109),i=1,6) /20,  8, 5.48E-02, 2.00E+00, 6.90E+05, 4.52E+06/
      data (drec_a(i,110),i=1,6) /20,  7, 7.13E-02, 1.82E+00, 6.70E+05, 3.32E+06/
      data (drec_a(i,111),i=1,6) /20,  6, 1.09E-01, 1.74E+00, 7.00E+05, 4.93E+06/
      data (drec_a(i,112),i=1,6) /20,  5, 1.10E-01, 2.43E-01, 5.67E+05, 4.41E+06/
      data (drec_a(i,113),i=1,6) /20,  4, 2.05E-02, 1.85E-01, 4.21E+05, 2.27E+06/
      data (drec_a(i,114),i=1,6) /20,  3, 5.49E-01, 2.92E-01, 3.65E+07, 7.25E+06/
      data (drec_a(i,115),i=1,6) /20,  2, 2.68E-01, 0.00E+00, 3.74E+07, 1.00E+05/
      data (drec_a(i,116),i=1,6) /26, 26, 1.58E-03, 4.56E-01, 6.00E+04, 8.97E+04/
      data (drec_a(i,117),i=1,6) /26, 25, 8.38E-03, 3.23E-01, 1.94E+05, 1.71E+05/
      data (drec_a(i,118),i=1,6) /26, 24, 1.54E-02, 3.10E-01, 3.31E+05, 2.73E+05/
      data (drec_a(i,119),i=1,6) /26, 23, 3.75E-02, 4.11E-01, 4.32E+05, 3.49E+05/
      data (drec_a(i,120),i=1,6) /26, 22, 1.17E-01, 3.59E-01, 6.28E+05, 5.29E+05/
      data (drec_a(i,121),i=1,6) /26, 21, 2.54E-01, 9.75E-02, 7.50E+05, 4.69E+05/
      data (drec_a(i,122),i=1,6) /26, 20, 2.91E-01, 2.29E-01, 7.73E+05, 6.54E+05/
      data (drec_a(i,123),i=1,6) /26, 19, 1.50E-01, 4.20E+00, 2.62E+05, 1.32E+06/
      data (drec_a(i,124),i=1,6) /26, 18, 1.40E-01, 3.30E+00, 2.50E+05, 1.33E+06/
      data (drec_a(i,125),i=1,6) /26, 17, 1.00E-01, 5.30E+00, 2.57E+05, 1.41E+06/
      data (drec_a(i,126),i=1,6) /26, 16, 2.00E-01, 1.50E+00, 2.84E+05, 1.52E+06/
      data (drec_a(i,127),i=1,6) /26, 15, 2.40E-01, 7.00E-01, 8.69E+05, 1.51E+06/
      data (drec_a(i,128),i=1,6) /26, 14, 2.60E-01, 6.00E-01, 4.21E+05, 1.82E+06/
      data (drec_a(i,129),i=1,6) /26, 13, 1.90E-01, 5.00E-01, 4.57E+05, 1.84E+06/
      data (drec_a(i,130),i=1,6) /26, 12, 1.20E-01, 1.00E+00, 2.85E+05, 2.31E+06/
      data (drec_a(i,131),i=1,6) /26, 11, 3.50E-01, 0.00E+00, 8.18E+06, 1.00E+05/
      data (drec_a(i,132),i=1,6) /26, 10, 6.60E-02, 7.80E+00, 1.51E+06, 9.98E+06/
      data (drec_a(i,133),i=1,6) /26,  9, 1.00E-01, 6.30E+00, 1.30E+06, 9.98E+06/
      data (drec_a(i,134),i=1,6) /26,  8, 1.30E-01, 5.50E+00, 1.19E+06, 1.00E+07/
      data (drec_a(i,135),i=1,6) /26,  7, 2.30E-01, 3.60E+00, 1.09E+06, 1.10E+07/
      data (drec_a(i,136),i=1,6)/26,  6, 1.40E-01, 4.90E+00, 9.62E+05, 8.34E+06/
      data (drec_a(i,137),i=1,6) /26,  5, 1.10E-01, 1.60E+00, 7.23E+05, 1.01E+07/
      data (drec_a(i,138),i=1,6) /26,  4, 4.10E-02, 4.20E+00, 4.23E+05, 1.07E+07/
      data (drec_a(i,139),i=1,6) /26,  3, 7.47E-01, 2.84E-01, 5.84E+07, 1.17E+07/
      data (drec_a(i,140),i=1,6) /26,  2, 3.69E-01, 0.00E+00, 6.00E+07, 1.00E+05/
      data (drec_a(i,141),i=1,6) /28, 28, 1.41E-03, 4.69E-01, 9.82E+04, 1.01E+05/
      data (drec_a(i,142),i=1,6) /28, 27, 5.20E-03, 3.57E-01, 2.01E+05, 1.91E+05/
      data (drec_a(i,143),i=1,6) /28, 26, 1.38E-02, 2.81E-01, 3.05E+05, 2.32E+05/
      data (drec_a(i,144),i=1,6) /28, 25, 2.30E-02, 1.28E-01, 4.20E+05, 3.18E+05/
      data (drec_a(i,145),i=1,6) /28, 24, 4.19E-02, 4.17E-02, 5.56E+05, 4.55E+05/
      data (drec_a(i,146),i=1,6) /28, 23, 6.83E-02, 5.58E-02, 6.72E+05, 5.51E+05/
      data (drec_a(i,147),i=1,6) /28, 22, 1.22E-01, 3.46E-02, 7.93E+05, 5.28E+05/
      data (drec_a(i,148),i=1,6) /28, 21, 3.00E-01, 0.00E+00, 9.00E+05, 1.00E+05/
      data (drec_a(i,149),i=1,6) /28, 20, 1.50E-01, 1.90E+00, 1.00E+06, 5.50E+05/
      data (drec_a(i,150),i=1,6) /28, 19, 6.97E-01, 2.77E-01, 7.81E+05, 8.87E+05/
      data (drec_a(i,151),i=1,6) /28, 18, 7.09E-01, 1.35E-01, 7.64E+05, 1.80E+06/
      data (drec_a(i,152),i=1,6) /28, 17, 6.44E-01, 1.34E-01, 7.44E+05, 1.25E+06/
      data (drec_a(i,153),i=1,6) /28, 16, 5.25E-01, 1.92E-01, 6.65E+05, 1.89E+06/
      data (drec_a(i,154),i=1,6) /28, 15, 4.46E-01, 3.32E-01, 5.97E+05, 8.84E+05/
      data (drec_a(i,155),i=1,6) /28, 14, 3.63E-01, 3.37E-01, 5.24E+05, 1.29E+06/
      data (drec_a(i,156),i=1,6) /28, 13, 3.02E-01, 1.21E-01, 4.96E+05, 6.24E+05/
      data (drec_a(i,157),i=1,6) /28, 12, 1.02E-01, 5.14E-02, 4.46E+05, 1.59E+06/
      data (drec_a(i,158),i=1,6) /28, 11, 2.70E-01, 1.83E-01, 8.49E+06, 8.01E+06/
      data (drec_a(i,159),i=1,6) /28, 10, 4.67E-02, 7.56E+00, 1.36E+06, 9.32E+06/
      data (drec_a(i,160),i=1,6) /28,  9, 8.35E-02, 4.55E+00, 1.23E+06, 9.45E+06/
      data (drec_a(i,161),i=1,6) /28,  8, 9.96E-02, 4.87E+00, 1.06E+06, 9.45E+06/
      data (drec_a(i,162),i=1,6) /28,  7, 1.99E-01, 2.19E+00, 1.25E+06, 8.01E+06/
      data (drec_a(i,163),i=1,6) /28,  6, 2.40E-01, 1.15E+00, 1.23E+06, 7.57E+06/
      data (drec_a(i,164),i=1,6) /28,  5, 1.15E-01, 1.23E+00, 3.32E+05, 2.64E+06/
      data (drec_a(i,165),i=1,6) /28,  4, 3.16E-02, 1.32E-01, 6.45E+05, 1.93E+06/
      data (drec_a(i,166),i=1,6) /28,  3, 8.03E-01, 2.89E-01, 6.65E+07, 1.19E+07/
      data (drec_a(i,167),i=1,6) /28,  2, 5.75E-01, 2.86E-01, 6.81E+07, 9.08E+06/
       end 



!=====================================================================
      subroutine drfit2_b(iz,in,t,d)
!---------------------------------------------------------------------
!     dielectronic recombination rate for 2*10^4 K <  T < 6*10^4 K
!---------------------------------------------------------------------
      implicit none
      integer :: iz,in
      real :: t,d,tt
      integer :: n,i
      integer, parameter :: nm=10000
      integer :: iza(nm),ina(nm)
      real :: aa(nm),ba(nm),t0a(nm),t1a(nm),fa(nm)
      integer,parameter :: ndr=23
      real :: drec_b(7,ndr)
      common/drec_b/drec_b

!      logical start/.true./
      logical start
      data start/.true./
      save start
         
      !To use units of t = T[K]/10**4
      tt = t/10**4

      if(start) then
         n=ndr
         do i=1,n
            iza(i)=nint(drec_b(1,i))
            ina(i)=nint(drec_b(2,i))
            aa(i)=drec_b(3,i)
            ba(i)=drec_b(4,i)
            t0a(i)=drec_b(5,i)
            t1a(i)=drec_b(6,i)
            fa(i)=drec_b(7,i)
         end do 
         start=.false.
      endif
      
      d=0
      do i=1,n
         if(iza(i).eq.iz.and.ina(i).eq.in) goto 200
      end do 
      return


 200  continue
      d=10**(-12)*(aa(i)/tt+ba(i)+t0a(i)*tt+t1a(i)*tt**2)*tt**(-3./2.)*exp(-fa(i)/tt)

      return
      end subroutine drfit2_b


 !==========================================================================
      block data derec_b
!--------------------------------------------------------------------------
!     dielectronic recombination data (2*10^4 K <  T < 6*10^4 K)
!--------------------------------------------------------------------------
      integer, parameter :: ndr=23
      real :: drec_b(7,ndr)
      common/drec_b/drec_b
data (drec_b(i,1),i=1,7) /6,  6,   0.0108,  -0.1075,   0.2810,  -0.0193,  -0.1127/
data (drec_b(i,2),i=1,7) /6,  5,   1.8267,   4.1012,   4.8443,   0.2261,   0.5960/
data (drec_b(i,3),i=1,7) /6,  4,   2.3196,  10.7328,   6.8830,  -0.1824,   0.4101/
data (drec_b(i,4),i=1,7) /7,  6,   0.0320,  -0.6624,   4.3191,   0.0003,   0.5946/
data (drec_b(i,5),i=1,7) /7,  5,  -0.8806,  11.2406,  30.7066,  -1.1721,   0.6127/
data (drec_b(i,6),i=1,7) /7,  4,   0.4134,  -4.6319,  25.9172,  -2.2290,   0.2360/
data (drec_b(i,7),i=1,7) /8,  8,   0.3715,  -0.0239,  -0.0597,   0.0678,   0.7993/
data (drec_b(i,8),i=1,7) /8,  7,  -0.0036,   0.7519,   1.5252,  -0.0838,   0.2769/
data (drec_b(i,9),i=1,7) /8,  6,   0.0000,  21.8790,  16.2730,  -0.7020,   1.1899/
data (drec_b(i,10),i=1,7) /8,  5,  -2.5053,   3.4903,  67.4128,  -3.4450,   0.8501/
data (drec_b(i,11),i=1,7) /8,  4,  -2.8425,   0.2283,  40.4072,  -3.4956,   1.7558/
data (drec_b(i,12),i=1,7) /10,  9,   0.0129,  -0.1779,   0.9353,  -0.0682,   0.4516/
data (drec_b(i,13),i=1,7) /10,  8,   3.6781,  14.1481,  17.1175,  -0.5017,   0.2313/
data (drec_b(i,14),i=1,7) /10,  7,  -0.0254,   5.5365,  17.0727,  -0.7225,   0.1702/
data (drec_b(i,15),i=1,7) /10,  6,  -0.0141,  33.8479,  43.1608,  -1.6072,   0.1942/
data (drec_b(i,16),i=1,7) /10,  5,  19.9280, 235.0536, 152.5096,   9.1413,   0.1282/
data (drec_b(i,17),i=1,7) /10,  4,   5.4751, 203.9751,  86.9016,  -7.4568,   2.5145/
data (drec_b(i,18),i=1,7) /12, 12,   1.2044,  -4.6836,   7.6620,  -0.5930,   1.6260/
data (drec_b(i,19),i=1,7) /13, 13,   0.0219,  -0.4528,   2.5427,  -0.1678,   0.2276/
data (drec_b(i,20),i=1,7) /13, 12,   0.7086,  -3.1083,   7.0422,   0.5998,   0.4194/
data (drec_b(i,21),i=1,7) /14, 14,  -0.0219,   0.4364,   0.0684,  -0.0032,   0.1342/
data (drec_b(i,22),i=1,7) /14, 13,   3.2163, -12.0571,  16.2118,  -0.5886,   0.5613/
data (drec_b(i,23),i=1,7) /14, 12,   0.1203,  -2.6900,  19.1943,  -0.1479,   0.1118/
       end 


!=====================================================================
      subroutine drfit2_c(iz,in,t,d)
!---------------------------------------------------------------------
!     dielectronic recombination rate for 1*10^3 K <  T < 2*10^4 K
!---------------------------------------------------------------------
      implicit none
      integer :: iz,in
      real :: t,d,tt
      integer :: n,i
      integer, parameter :: nm=10000
      integer :: iza(nm),ina(nm)
      real :: aa(nm),ba(nm),t0a(nm),t1a(nm),fa(nm)
      integer,parameter :: ndr=23
      real :: drec_c(7,ndr)
      common/drec_c/drec_c

!      logical start/.true./
      logical start
      data start/.true./
      save start
         
      !To use units of t = T[K]/10**4
      tt = t/10**4

      if(start) then
         n=ndr
         do i=1,n
            iza(i)=nint(drec_c(1,i))
            ina(i)=nint(drec_c(2,i))
            aa(i)=drec_c(3,i)
            ba(i)=drec_c(4,i)
            t0a(i)=drec_c(5,i)
            t1a(i)=drec_c(6,i)
            fa(i)=drec_c(7,i)
         end do 
         start=.false.
      endif
      
      d=0
      do i=1,n
         if(iza(i).eq.iz.and.ina(i).eq.in) goto 200
      end do 
      return


 200  continue
      d=10**(-12)*(aa(i)/tt+ba(i)+t0a(i)*tt+t1a(i)*tt**2)*tt**(-3./2.)*exp(-fa(i)/tt)

      return
      end subroutine drfit2_c


 !==========================================================================
      block data derec_c
!--------------------------------------------------------------------------
!     dielectronic recombination data (1*10^3 K <  T < 2*10^4 K)
!--------------------------------------------------------------------------
      integer, parameter :: ndr=23
      real :: drec_c(7,ndr)
      common/drec_c/drec_c
data (drec_c(i,1),i=1,7) /6,  6,   0.0108,  -0.1075,   0.2810,  -0.0193,  -0.1127/
data (drec_c(i,2),i=1,7) /6,  5,   1.8267,   4.1012,   4.8443,   0.2261,   0.5960/
data (drec_c(i,3),i=1,7) /6,  4,   2.3196,  10.7328,   6.8830,  -0.1824,   0.4101/
data (drec_c(i,4),i=1,7) /7,  6,   0.0320,  -0.6624,   4.3191,   0.0003,   0.5946/
data (drec_c(i,5),i=1,7) /7,  5,  -0.8806,  11.2406,  30.7066,  -1.1721,   0.6127/
data (drec_c(i,6),i=1,7) /7,  4,   0.4134,  -4.6319,  25.9172,  -2.2290,   0.2360/
data (drec_c(i,7),i=1,7) /8,  8,  -0.0001,   0.0001,   0.0956,   0.0193,   0.4106/
data (drec_c(i,8),i=1,7) /8,  7,  -0.0036,   0.7519,   1.5252,  -0.0838,   0.2769/
data (drec_c(i,9),i=1,7) /8,  6,   0.0000,  21.8790,  16.2730,  -0.7020,   1.1899/
data (drec_c(i,10),i=1,7) /8,  5,  -0.3648,   7.2698,  17.2187,   9.8335,  -0.0166/
data (drec_c(i,11),i=1,7) /8,  4,  -2.8425,   0.2283,  40.4072,  -3.4956,   1.7558/
data (drec_c(i,12),i=1,7) /10,  9,   0.0129,  -0.1779,   0.9353,  -0.0682,   0.4516/
data (drec_c(i,13),i=1,7) /10,  8,   3.6781,  14.1481,  17.1175,  -0.5017,   0.2313/
data (drec_c(i,14),i=1,7) /10,  7,  -0.0254,   5.5365,  17.0727,  -0.7225,   0.1702/
data (drec_c(i,15),i=1,7) /10,  6,  -0.0141,  33.8479,  43.1608,  -1.6072,   0.1942/
data (drec_c(i,16),i=1,7) /10,  5,  19.9280, 235.0536, 152.5096,   9.1413,   0.1282/
data (drec_c(i,17),i=1,7) /10,  4,   5.4751, 203.9751,  86.9016,  -7.4568,   2.5145/
data (drec_c(i,18),i=1,7) /12, 12,   1.2044,  -4.6836,   7.6620,  -0.5930,   1.6260/
data (drec_c(i,19),i=1,7) /13, 13,   0.0219,  -0.4528,   2.5427,  -0.1678,   0.2276/
data (drec_c(i,20),i=1,7) /13, 12,   0.7086,  -3.1083,   7.0422,   0.5998,   0.4194/
data (drec_c(i,21),i=1,7) /14, 14,  -0.0219,   0.4364,   0.0684,  -0.0032,   0.1342/
data (drec_c(i,22),i=1,7) /14, 13,   3.2163, -12.0571,  16.2118,  -0.5886,   0.5613/
data (drec_c(i,23),i=1,7) /14, 12,   0.1203,  -2.6900,  19.1943,  -0.1479,   0.1118/
       end    

                        
!============================================================================
      subroutine cfit_b(iz,in,t,c)
!* Version 2, March 24, 1997
!******************************************************************************
!*** This subroutine calculates rates of direct collisional ionization 
!*** for all ionization stages of all elements from H to Ni (Z=28)
!*** by use of the fits from G. S. Voronov, 1997, ADNDT, 65, 1
!*** Input parameters:  iz - atomic number 
!***                    in - number of electrons from 1 to iz 
!***                    t  - temperature, K
!*** Output parameter:  c  - rate coefficient, cm^3 s^(-1)
!******************************************************************************
      common/CF/CF(5,28,28)
      c=0.0
      if(iz.lt.1.or.iz.gt.28)return
      if(in.lt.1.or.in.gt.iz)return
      te=t*8.617385e-05
      u=cf(1,iz,in)/te
      if(u.gt.80.0)return
      c=cf(3,iz,in)*(1.0+cf(2,iz,in)*sqrt(u))/(cf(4,iz,in)+u)* &
      u**cf(5,iz,in)*exp(-u)
      return
      end subroutine cfit_b
!*********************************************
      block data cifit
      common/CF/CF(5,28,28)
      DATA(CF(I, 1, 1),I=1,5)/   13.6,0.,2.91E-08,0.2320,0.39/
      DATA(CF(I, 2, 2),I=1,5)/   24.6,0.,1.75E-08,0.1800,0.35/
      DATA(CF(I, 2, 1),I=1,5)/   54.4,1.,2.05E-09,0.2650,0.25/
      DATA(CF(I, 3, 3),I=1,5)/    5.4,0.,1.39E-07,0.4380,0.41/
      DATA(CF(I, 3, 2),I=1,5)/   75.6,1.,2.01E-09,0.2090,0.23/
      DATA(CF(I, 3, 1),I=1,5)/  122.4,1.,9.60E-10,0.5820,0.17/
      DATA(CF(I, 4, 4),I=1,5)/    9.3,0.,1.02E-07,0.3750,0.27/
      DATA(CF(I, 4, 3),I=1,5)/   18.2,1.,2.08E-08,0.4390,0.21/
      DATA(CF(I, 4, 2),I=1,5)/  153.9,0.,2.67E-09,0.6120,0.27/
      DATA(CF(I, 4, 1),I=1,5)/  217.7,1.,4.27E-10,0.6580,0.15/
      DATA(CF(I, 5, 5),I=1,5)/    8.3,0.,6.49E-08,0.2000,0.26/
      DATA(CF(I, 5, 4),I=1,5)/   25.2,1.,1.24E-08,0.2670,0.22/
      DATA(CF(I, 5, 3),I=1,5)/   37.9,1.,3.27E-09,0.2950,0.23/
      DATA(CF(I, 5, 2),I=1,5)/  259.4,1.,4.95E-10,0.4890,0.09/
      DATA(CF(I, 5, 1),I=1,5)/  340.2,1.,2.19E-10,0.6570,0.15/
      DATA(CF(I, 6, 6),I=1,5)/   11.3,0.,6.85E-08,0.1930,0.25/
      DATA(CF(I, 6, 5),I=1,5)/   24.4,1.,1.86E-08,0.2860,0.24/
      DATA(CF(I, 6, 4),I=1,5)/   47.9,1.,6.35E-09,0.4270,0.21/
      DATA(CF(I, 6, 3),I=1,5)/   64.5,1.,1.50E-09,0.4160,0.13/
      DATA(CF(I, 6, 2),I=1,5)/  392.1,1.,2.99E-10,0.6660,0.02/
      DATA(CF(I, 6, 1),I=1,5)/  490.0,1.,1.23E-10,0.6200,0.16/
      DATA(CF(I, 7, 7),I=1,5)/   14.5,0.,4.82E-08,0.0652,0.42/
      DATA(CF(I, 7, 6),I=1,5)/   29.6,0.,2.98E-08,0.3100,0.30/
      DATA(CF(I, 7, 5),I=1,5)/   47.5,1.,8.10E-09,0.3500,0.24/
      DATA(CF(I, 7, 4),I=1,5)/   77.5,1.,3.71E-09,0.5490,0.18/
      DATA(CF(I, 7, 3),I=1,5)/   97.9,0.,1.51E-09,0.0167,0.74/
      DATA(CF(I, 7, 2),I=1,5)/  552.1,0.,3.71E-10,0.5460,0.29/
      DATA(CF(I, 7, 1),I=1,5)/  667.0,1.,7.77E-11,0.6240,0.16/
      DATA(CF(I, 8, 8),I=1,5)/   13.6,0.,3.59E-08,0.0730,0.34/
      DATA(CF(I, 8, 7),I=1,5)/   35.1,1.,1.39E-08,0.2120,0.22/
      DATA(CF(I, 8, 6),I=1,5)/   54.9,1.,9.31E-09,0.2700,0.27/
      DATA(CF(I, 8, 5),I=1,5)/   77.4,0.,1.02E-08,0.6140,0.27/
      DATA(CF(I, 8, 4),I=1,5)/  113.9,1.,2.19E-09,0.6300,0.17/
      DATA(CF(I, 8, 3),I=1,5)/  138.1,0.,1.95E-09,0.3600,0.54/
      DATA(CF(I, 8, 2),I=1,5)/  739.3,0.,2.12E-10,0.3960,0.35/
      DATA(CF(I, 8, 1),I=1,5)/  871.4,1.,5.21E-11,0.6290,0.16/
      DATA(CF(I, 9, 9),I=1,5)/   17.4,1.,7.00E-08,0.1780,0.29/
      DATA(CF(I, 9, 8),I=1,5)/   35.0,0.,5.41E-08,0.5710,0.27/
      DATA(CF(I, 9, 7),I=1,5)/   62.7,1.,9.37E-09,0.3190,0.20/
      DATA(CF(I, 9, 6),I=1,5)/   87.1,1.,4.92E-09,0.3230,0.24/
      DATA(CF(I, 9, 5),I=1,5)/  114.2,0.,7.06E-09,0.6840,0.27/
      DATA(CF(I, 9, 4),I=1,5)/  157.2,1.,1.28E-09,0.6480,0.16/
      DATA(CF(I, 9, 3),I=1,5)/  185.2,1.,5.61E-10,0.7380,0.16/
      DATA(CF(I, 9, 2),I=1,5)/  953.9,0.,1.66E-10,0.5420,0.29/
      DATA(CF(I, 9, 1),I=1,5)/ 1103.1,1.,3.74E-11,0.6590,0.15/
      DATA(CF(I,10,10),I=1,5)/   21.6,1.,1.50E-08,0.0329,0.43/
      DATA(CF(I,10, 9),I=1,5)/   41.0,0.,1.98E-08,0.2950,0.20/
      DATA(CF(I,10, 8),I=1,5)/   63.5,1.,7.03E-09,0.0677,0.39/
      DATA(CF(I,10, 7),I=1,5)/   97.1,1.,4.24E-09,0.0482,0.58/
      DATA(CF(I,10, 6),I=1,5)/  126.2,1.,2.79E-09,0.3050,0.25/
      DATA(CF(I,10, 5),I=1,5)/  157.9,0.,3.45E-09,0.5810,0.28/
      DATA(CF(I,10, 4),I=1,5)/  207.3,1.,9.56E-10,0.7490,0.14/
      DATA(CF(I,10, 3),I=1,5)/  239.1,1.,4.73E-10,0.9920,0.04/
      DATA(CF(I,10, 2),I=1,5)/ 1196.0,1.,3.92E-11,0.2620,0.20/
      DATA(CF(I,10, 1),I=1,5)/ 1360.6,1.,2.77E-11,0.6610,0.13/
      DATA(CF(I,11,11),I=1,5)/    5.1,1.,1.01E-07,0.2750,0.23/
      DATA(CF(I,11,10),I=1,5)/   47.3,1.,7.35E-09,0.0560,0.35/
      DATA(CF(I,11, 9),I=1,5)/   71.6,1.,8.10E-09,0.1480,0.32/
      DATA(CF(I,11, 8),I=1,5)/   98.9,0.,1.14E-08,0.5530,0.28/
      DATA(CF(I,11, 7),I=1,5)/  138.4,1.,2.63E-09,0.2300,0.29/
      DATA(CF(I,11, 6),I=1,5)/  172.2,1.,1.85E-09,0.3630,0.22/
      DATA(CF(I,11, 5),I=1,5)/  208.5,0.,2.82E-09,0.6740,0.27/
      DATA(CF(I,11, 4),I=1,5)/  264.2,1.,6.72E-10,0.7520,0.14/
      DATA(CF(I,11, 3),I=1,5)/  299.9,1.,2.80E-10,0.7810,0.15/
      DATA(CF(I,11, 2),I=1,5)/ 1465.1,1.,4.63E-11,0.5580,0.16/
      DATA(CF(I,11, 1),I=1,5)/ 1648.7,1.,2.16E-11,0.7430,0.13/
      DATA(CF(I,12,12),I=1,5)/    7.6,0.,6.21E-07,0.5920,0.39/
      DATA(CF(I,12,11),I=1,5)/   15.2,0.,1.92E-08,0.0027,0.85/
      DATA(CF(I,12,10),I=1,5)/   80.1,1.,5.56E-09,0.1070,0.30/
      DATA(CF(I,12, 9),I=1,5)/  109.3,1.,4.35E-09,0.1590,0.31/
      DATA(CF(I,12, 8),I=1,5)/  141.3,0.,7.10E-09,0.6580,0.25/
      DATA(CF(I,12, 7),I=1,5)/  186.5,1.,1.70E-09,0.2420,0.28/
      DATA(CF(I,12, 6),I=1,5)/  224.9,1.,1.22E-09,0.3430,0.23/
      DATA(CF(I,12, 5),I=1,5)/  266.0,0.,2.20E-09,0.8970,0.22/
      DATA(CF(I,12, 4),I=1,5)/  328.2,1.,4.86E-10,0.7510,0.14/
      DATA(CF(I,12, 3),I=1,5)/  367.5,1.,2.35E-10,1.0300,0.10/
      DATA(CF(I,12, 2),I=1,5)/ 1761.8,1.,2.06E-11,0.1960,0.25/
      DATA(CF(I,12, 1),I=1,5)/ 1962.7,1.,1.75E-11,0.8350,0.11/
      DATA(CF(I,13,13),I=1,5)/    6.0,1.,2.28E-07,0.3870,0.25/
      DATA(CF(I,13,12),I=1,5)/   18.8,0.,1.18E-07,2.2100,0.25/
      DATA(CF(I,13,11),I=1,5)/   28.5,1.,4.40E-09,0.1060,0.24/
      DATA(CF(I,13,10),I=1,5)/  120.0,0.,1.75E-08,0.8720,0.22/
      DATA(CF(I,13, 9),I=1,5)/  153.8,1.,2.61E-09,0.1590,0.31/
      DATA(CF(I,13, 8),I=1,5)/  198.5,1.,1.85E-09,0.1520,0.36/
      DATA(CF(I,13, 7),I=1,5)/  241.4,1.,1.14E-09,0.2280,0.29/
      DATA(CF(I,13, 6),I=1,5)/  284.6,1.,8.00E-10,0.4170,0.16/
      DATA(CF(I,13, 5),I=1,5)/  390.2,1.,5.83E-10,0.4970,0.23/
      DATA(CF(I,13, 4),I=1,5)/  399.4,0.,4.93E-10,0.7060,0.16/
      DATA(CF(I,13, 3),I=1,5)/  442.0,1.,9.77E-11,0.2780,0.17/
      DATA(CF(I,13, 2),I=1,5)/ 2086.6,0.,3.94E-11,0.2860,0.36/
      DATA(CF(I,13, 1),I=1,5)/ 2304.1,1.,1.38E-11,0.8350,0.11/
      DATA(CF(I,14,14),I=1,5)/    8.2,1.,1.88E-07,0.3760,0.25/
      DATA(CF(I,14,13),I=1,5)/   16.4,1.,6.43E-08,0.6320,0.20/
      DATA(CF(I,14,12),I=1,5)/   33.5,1.,2.01E-08,0.4730,0.22/
      DATA(CF(I,14,11),I=1,5)/   54.0,1.,4.94E-09,0.1720,0.23/
      DATA(CF(I,14,10),I=1,5)/  166.8,1.,1.76E-09,0.1020,0.31/
      DATA(CF(I,14, 9),I=1,5)/  205.3,1.,1.74E-09,0.1800,0.29/
      DATA(CF(I,14, 8),I=1,5)/  246.5,1.,1.23E-09,0.5180,0.07/
      DATA(CF(I,14, 7),I=1,5)/  303.5,1.,8.27E-10,0.2390,0.28/
      DATA(CF(I,14, 6),I=1,5)/  351.1,1.,6.01E-10,0.3050,0.25/
      DATA(CF(I,14, 5),I=1,5)/  401.4,1.,4.65E-10,0.6660,0.04/
      DATA(CF(I,14, 4),I=1,5)/  476.4,1.,2.63E-10,0.6660,0.16/
      DATA(CF(I,14, 3),I=1,5)/  523.5,1.,1.18E-10,0.7340,0.16/
      DATA(CF(I,14, 2),I=1,5)/ 2437.7,0.,3.36E-11,0.3360,0.37/
      DATA(CF(I,14, 1),I=1,5)/ 2673.2,1.,1.19E-11,0.9890,0.08/
      DATA(CF(I,15,15),I=1,5)/   10.5,1.,1.99E-07,0.5350,0.24/
      DATA(CF(I,15,14),I=1,5)/   19.8,1.,5.88E-08,0.5370,0.21/
      DATA(CF(I,15,13),I=1,5)/   30.2,1.,2.96E-08,0.8650,0.16/
      DATA(CF(I,15,12),I=1,5)/   51.4,1.,1.01E-08,0.5460,0.20/
      DATA(CF(I,15,11),I=1,5)/   65.0,1.,2.36E-09,0.1920,0.17/
      DATA(CF(I,15,10),I=1,5)/  220.4,0.,6.66E-09,1.0000,0.18/
      DATA(CF(I,15, 9),I=1,5)/  263.2,1.,1.24E-09,0.2150,0.26/
      DATA(CF(I,15, 8),I=1,5)/  309.4,0.,2.27E-09,0.7340,0.23/
      DATA(CF(I,15, 7),I=1,5)/  371.7,1.,6.14E-10,0.2560,0.27/
      DATA(CF(I,15, 6),I=1,5)/  424.5,1.,4.69E-10,0.3420,0.23/
      DATA(CF(I,15, 5),I=1,5)/  479.6,0.,6.14E-10,0.3340,0.39/
      DATA(CF(I,15, 4),I=1,5)/  560.4,0.,3.22E-10,0.8500,0.12/
      DATA(CF(I,15, 3),I=1,5)/  611.9,1.,9.32E-11,0.7340,0.16/
      DATA(CF(I,15, 2),I=1,5)/ 2816.9,0.,3.79E-11,0.8050,0.22/
      DATA(CF(I,15, 1),I=1,5)/ 3069.9,1.,9.73E-12,0.9910,0.08/
      DATA(CF(I,16,16),I=1,5)/   10.4,1.,5.49E-08,0.1000,0.25/
      DATA(CF(I,16,15),I=1,5)/   23.3,1.,6.81E-08,0.6930,0.21/
      DATA(CF(I,16,14),I=1,5)/   34.8,1.,2.14E-08,0.3530,0.24/
      DATA(CF(I,16,13),I=1,5)/   47.3,1.,1.66E-08,1.0300,0.14/
      DATA(CF(I,16,12),I=1,5)/   72.6,1.,6.12E-09,0.5800,0.19/
      DATA(CF(I,16,11),I=1,5)/   88.1,1.,1.33E-09,0.0688,0.35/
      DATA(CF(I,16,10),I=1,5)/  280.9,0.,4.93E-09,1.1300,0.16/
      DATA(CF(I,16, 9),I=1,5)/  328.2,1.,8.73E-10,0.1930,0.28/
      DATA(CF(I,16, 8),I=1,5)/  379.1,0.,1.35E-09,0.4310,0.32/
      DATA(CF(I,16, 7),I=1,5)/  447.1,1.,4.59E-10,0.2420,0.28/
      DATA(CF(I,16, 6),I=1,5)/  504.8,1.,3.49E-10,0.3050,0.25/
      DATA(CF(I,16, 5),I=1,5)/  564.7,0.,5.23E-10,0.4280,0.35/
      DATA(CF(I,16, 4),I=1,5)/  651.6,0.,2.59E-10,0.8540,0.12/
      DATA(CF(I,16, 3),I=1,5)/  707.2,1.,7.50E-11,0.7340,0.16/
      DATA(CF(I,16, 2),I=1,5)/ 3223.9,0.,2.67E-11,0.5720,0.28/
      DATA(CF(I,16, 1),I=1,5)/ 3494.2,1.,6.32E-12,0.5850,0.17/
      DATA(CF(I,17,17),I=1,5)/   13.0,1.,1.69E-07,0.4300,0.24/
      DATA(CF(I,17,16),I=1,5)/   23.8,1.,6.96E-08,0.6700,0.20/
      DATA(CF(I,17,15),I=1,5)/   39.6,1.,3.40E-08,0.8650,0.18/
      DATA(CF(I,17,14),I=1,5)/   53.5,1.,1.10E-08,0.3280,0.25/
      DATA(CF(I,17,13),I=1,5)/   67.8,1.,1.11E-08,1.3700,0.10/
      DATA(CF(I,17,12),I=1,5)/   97.0,1.,3.17E-09,0.3300,0.24/
      DATA(CF(I,17,11),I=1,5)/  114.2,1.,1.01E-09,0.1960,0.16/
      DATA(CF(I,17,10),I=1,5)/  348.3,0.,2.11E-09,0.3130,0.37/
      DATA(CF(I,17, 9),I=1,5)/  400.1,1.,6.32E-10,0.1730,0.30/
      DATA(CF(I,17, 8),I=1,5)/  455.6,0.,9.48E-10,0.3440,0.36/
      DATA(CF(I,17, 7),I=1,5)/  529.3,1.,3.69E-10,0.2730,0.26/
      DATA(CF(I,17, 6),I=1,5)/  592.0,1.,2.85E-10,0.3430,0.23/
      DATA(CF(I,17, 5),I=1,5)/  656.7,0.,4.81E-10,0.6580,0.27/
      DATA(CF(I,17, 4),I=1,5)/  749.8,1.,1.31E-10,0.6230,0.16/
      DATA(CF(I,17, 3),I=1,5)/  809.4,1.,6.13E-11,0.7360,0.16/
      DATA(CF(I,17, 2),I=1,5)/ 3658.4,0.,1.90E-11,0.3790,0.36/
      DATA(CF(I,17, 1),I=1,5)/ 3946.3,1.,5.14E-12,0.5530,0.18/
      DATA(CF(I,18,18),I=1,5)/   15.8,1.,5.99E-08,0.1360,0.26/
      DATA(CF(I,18,17),I=1,5)/   27.6,1.,6.07E-08,0.5440,0.21/
      DATA(CF(I,18,16),I=1,5)/   40.9,1.,3.43E-08,0.8340,0.17/
      DATA(CF(I,18,15),I=1,5)/   52.3,0.,3.00E-08,1.0300,0.25/
      DATA(CF(I,18,14),I=1,5)/   75.0,1.,8.73E-09,0.3660,0.31/
      DATA(CF(I,18,13),I=1,5)/   91.0,1.,5.78E-09,0.3140,0.34/
      DATA(CF(I,18,12),I=1,5)/  124.3,1.,2.98E-09,0.7030,0.16/
      DATA(CF(I,18,11),I=1,5)/  143.5,1.,7.25E-10,0.2070,0.15/
      DATA(CF(I,18,10),I=1,5)/  422.4,1.,1.40E-09,0.6960,0.13/
      DATA(CF(I,18, 9),I=1,5)/  478.7,1.,4.78E-10,0.1640,0.31/
      DATA(CF(I,18, 8),I=1,5)/  539.0,0.,8.02E-10,0.4390,0.32/
      DATA(CF(I,18, 7),I=1,5)/  618.3,1.,2.88E-10,0.2590,0.27/
      DATA(CF(I,18, 6),I=1,5)/  686.1,1.,2.32E-10,0.3620,0.22/
      DATA(CF(I,18, 5),I=1,5)/  755.7,0.,3.33E-10,0.4120,0.36/
      DATA(CF(I,18, 4),I=1,5)/  854.8,1.,1.27E-10,0.9100,0.13/
      DATA(CF(I,18, 3),I=1,5)/  918.0,1.,5.21E-11,0.7810,0.15/
      DATA(CF(I,18, 2),I=1,5)/ 4120.7,0.,1.66E-11,0.4350,0.33/
      DATA(CF(I,18, 1),I=1,5)/ 4426.2,1.,4.32E-12,0.5540,0.18/
      DATA(CF(I,19,19),I=1,5)/    4.3,1.,2.02E-07,0.2720,0.31/
      DATA(CF(I,19,18),I=1,5)/   31.6,1.,4.01E-08,0.3710,0.22/
      DATA(CF(I,19,17),I=1,5)/   45.8,1.,1.50E-08,0.4330,0.21/
      DATA(CF(I,19,16),I=1,5)/   60.9,1.,1.94E-08,0.8890,0.16/
      DATA(CF(I,19,15),I=1,5)/   82.7,1.,6.95E-09,0.4940,0.18/
      DATA(CF(I,19,14),I=1,5)/   99.4,1.,4.11E-09,0.5400,0.17/
      DATA(CF(I,19,13),I=1,5)/  117.6,1.,2.23E-09,0.5190,0.16/
      DATA(CF(I,19,12),I=1,5)/  154.7,1.,2.15E-09,0.8280,0.14/
      DATA(CF(I,19,11),I=1,5)/  175.8,0.,1.61E-09,0.6420,0.13/
      DATA(CF(I,19,10),I=1,5)/  504.0,1.,1.07E-09,0.6950,0.13/
      DATA(CF(I,19, 9),I=1,5)/  564.7,1.,3.78E-10,0.1730,0.30/
      DATA(CF(I,19, 8),I=1,5)/  629.4,0.,6.24E-10,0.4180,0.33/
      DATA(CF(I,19, 7),I=1,5)/  714.6,1.,2.29E-10,0.2450,0.28/
      DATA(CF(I,19, 6),I=1,5)/  786.6,1.,1.86E-10,0.3440,0.23/
      DATA(CF(I,19, 5),I=1,5)/  861.1,0.,2.69E-10,0.3960,0.37/
      DATA(CF(I,19, 4),I=1,5)/  968.0,1.,1.06E-10,0.9120,0.13/
      DATA(CF(I,19, 3),I=1,5)/ 1053.4,1.,4.24E-11,0.7370,0.16/
      DATA(CF(I,19, 2),I=1,5)/ 4610.9,0.,1.38E-11,0.4160,0.34/
      DATA(CF(I,19, 1),I=1,5)/ 4934.1,1.,3.67E-12,0.5550,0.18/
      DATA(CF(I,20,20),I=1,5)/    6.1,0.,4.40E-07,0.8480,0.33/
      DATA(CF(I,20,19),I=1,5)/   11.9,0.,5.22E-08,0.1510,0.34/
      DATA(CF(I,20,18),I=1,5)/   50.9,1.,2.06E-08,0.4180,0.20/
      DATA(CF(I,20,17),I=1,5)/   67.3,1.,1.72E-08,0.6380,0.19/
      DATA(CF(I,20,16),I=1,5)/   84.5,1.,1.26E-08,1.0100,0.14/
      DATA(CF(I,20,15),I=1,5)/  108.8,1.,4.72E-09,0.5260,0.17/
      DATA(CF(I,20,14),I=1,5)/  127.2,1.,2.89E-09,0.5480,0.17/
      DATA(CF(I,20,13),I=1,5)/  147.2,1.,1.64E-09,0.5520,0.15/
      DATA(CF(I,20,12),I=1,5)/  188.3,1.,1.57E-09,0.7990,0.14/
      DATA(CF(I,20,11),I=1,5)/  211.3,1.,4.32E-10,0.2320,0.14/
      DATA(CF(I,20,10),I=1,5)/  591.9,0.,9.47E-10,0.3110,0.38/
      DATA(CF(I,20, 9),I=1,5)/  657.2,1.,2.98E-10,0.1630,0.31/
      DATA(CF(I,20, 8),I=1,5)/  726.6,0.,4.78E-10,0.3590,0.36/
      DATA(CF(I,20, 7),I=1,5)/  817.6,1.,1.86E-10,0.2440,0.28/
      DATA(CF(I,20, 6),I=1,5)/  894.5,1.,1.56E-10,0.3640,0.22/
      DATA(CF(I,20, 5),I=1,5)/  974.0,0.,2.16E-10,0.3570,0.39/
      DATA(CF(I,20, 4),I=1,5)/ 1087.0,1.,7.70E-11,0.6550,0.15/
      DATA(CF(I,20, 3),I=1,5)/ 1157.0,1.,3.58E-11,0.7360,0.16/
      DATA(CF(I,20, 2),I=1,5)/ 5128.9,0.,1.28E-11,0.5200,0.30/
      DATA(CF(I,20, 1),I=1,5)/ 5469.9,1.,3.08E-12,0.5280,0.19/
      DATA(CF(I,21,21),I=1,5)/    6.6,1.,3.16E-07,0.2040,0.28/
      DATA(CF(I,21,20),I=1,5)/   12.8,1.,8.61E-08,0.1810,0.25/
      DATA(CF(I,21,19),I=1,5)/   24.8,1.,5.08E-08,0.3570,0.24/
      DATA(CF(I,21,18),I=1,5)/   73.5,1.,1.00E-08,0.4530,0.15/
      DATA(CF(I,21,17),I=1,5)/   91.9,1.,6.76E-09,0.4600,0.15/
      DATA(CF(I,21,16),I=1,5)/  110.7,1.,5.27E-09,0.5610,0.17/
      DATA(CF(I,21,15),I=1,5)/  138.0,1.,3.40E-09,0.5600,0.16/
      DATA(CF(I,21,14),I=1,5)/  158.1,1.,2.18E-09,0.6120,0.15/
      DATA(CF(I,21,13),I=1,5)/  180.0,1.,1.26E-09,0.6100,0.14/
      DATA(CF(I,21,12),I=1,5)/  225.1,1.,1.24E-09,0.8520,0.13/
      DATA(CF(I,21,11),I=1,5)/  249.8,1.,3.62E-10,0.3490,0.05/
      DATA(CF(I,21,10),I=1,5)/  687.4,1.,5.52E-10,0.3750,0.28/
      DATA(CF(I,21, 9),I=1,5)/  756.7,1.,5.64E-10,0.8730,0.15/
      DATA(CF(I,21, 8),I=1,5)/  830.8,1.,4.50E-10,1.0500,0.13/
      DATA(CF(I,21, 7),I=1,5)/  927.5,1.,2.73E-10,0.8660,0.15/
      DATA(CF(I,21, 6),I=1,5)/ 1009.0,1.,1.56E-10,0.7150,0.17/
      DATA(CF(I,21, 5),I=1,5)/ 1094.0,0.,1.81E-10,1.1400,0.36/
      DATA(CF(I,21, 4),I=1,5)/ 1213.0,1.,4.29E-11,0.7840,0.15/
      DATA(CF(I,21, 3),I=1,5)/ 1288.0,0.,2.21E-11,0.0270,0.82/
      DATA(CF(I,21, 2),I=1,5)/ 5674.9,1.,4.51E-12,0.9180,0.04/
      DATA(CF(I,21, 1),I=1,5)/ 6033.8,0.,2.03E-12,0.0170,0.70/
      DATA(CF(I,22,22),I=1,5)/    6.8,1.,1.60E-07,0.3600,0.28/
      DATA(CF(I,22,21),I=1,5)/   13.6,0.,2.14E-07,0.8800,0.28/
      DATA(CF(I,22,20),I=1,5)/   27.5,1.,2.85E-08,0.2270,0.21/
      DATA(CF(I,22,19),I=1,5)/   43.3,1.,3.48E-08,0.3900,0.23/
      DATA(CF(I,22,18),I=1,5)/   99.3,1.,1.00E-08,0.5790,0.18/
      DATA(CF(I,22,17),I=1,5)/  119.5,1.,7.01E-09,0.6380,0.17/
      DATA(CF(I,22,16),I=1,5)/  140.8,1.,4.95E-09,0.7170,0.16/
      DATA(CF(I,22,15),I=1,5)/  170.4,1.,2.99E-09,0.6930,0.17/
      DATA(CF(I,22,14),I=1,5)/  192.1,1.,2.10E-09,0.7220,0.16/
      DATA(CF(I,22,13),I=1,5)/  215.9,1.,1.62E-09,0.7650,0.14/
      DATA(CF(I,22,12),I=1,5)/  265.0,1.,1.11E-09,0.8850,0.12/
      DATA(CF(I,22,11),I=1,5)/  291.5,0.,9.09E-10,0.9720,0.06/
      DATA(CF(I,22,10),I=1,5)/  787.8,1.,4.41E-10,0.3590,0.29/
      DATA(CF(I,22, 9),I=1,5)/  863.1,1.,4.39E-10,0.7810,0.17/
      DATA(CF(I,22, 8),I=1,5)/  941.9,1.,3.73E-10,1.0500,0.13/
      DATA(CF(I,22, 7),I=1,5)/ 1044.0,1.,2.28E-10,0.8580,0.15/
      DATA(CF(I,22, 6),I=1,5)/ 1131.0,1.,1.34E-10,0.7570,0.16/
      DATA(CF(I,22, 5),I=1,5)/ 1221.0,0.,1.55E-10,1.1500,0.36/
      DATA(CF(I,22, 4),I=1,5)/ 1346.0,1.,3.80E-11,0.8350,0.14/
      DATA(CF(I,22, 3),I=1,5)/ 1426.0,0.,1.89E-11,0.0280,0.82/
      DATA(CF(I,22, 2),I=1,5)/ 6249.1,1.,4.01E-12,0.9680,0.03/
      DATA(CF(I,22, 1),I=1,5)/ 6625.0,1.,1.62E-12,0.6570,0.14/
      DATA(CF(I,23,23),I=1,5)/    6.7,0.,8.82E-07,0.3590,0.32/
      DATA(CF(I,23,22),I=1,5)/   14.7,0.,3.11E-07,0.4320,0.29/
      DATA(CF(I,23,21),I=1,5)/   29.3,1.,3.50E-08,0.2470,0.25/
      DATA(CF(I,23,20),I=1,5)/   46.7,0.,5.32E-08,1.1100,0.16/
      DATA(CF(I,23,19),I=1,5)/   65.3,1.,8.98E-09,0.1400,0.37/
      DATA(CF(I,23,18),I=1,5)/  128.1,1.,5.87E-09,0.5170,0.17/
      DATA(CF(I,23,17),I=1,5)/  150.6,1.,5.11E-09,0.6790,0.16/
      DATA(CF(I,23,16),I=1,5)/  173.4,1.,3.71E-09,0.7610,0.15/
      DATA(CF(I,23,15),I=1,5)/  205.8,1.,2.24E-09,0.7110,0.17/
      DATA(CF(I,23,14),I=1,5)/  230.5,1.,1.65E-09,0.7640,0.15/
      DATA(CF(I,23,13),I=1,5)/  256.0,1.,1.26E-09,0.7620,0.14/
      DATA(CF(I,23,12),I=1,5)/  308.0,1.,8.86E-10,0.8860,0.12/
      DATA(CF(I,23,11),I=1,5)/  336.3,0.,3.89E-10,0.1420,0.39/
      DATA(CF(I,23,10),I=1,5)/  896.0,1.,3.80E-10,0.4090,0.27/
      DATA(CF(I,23, 9),I=1,5)/  976.0,0.,4.84E-10,0.1730,0.57/
      DATA(CF(I,23, 8),I=1,5)/ 1060.0,1.,2.49E-10,0.6500,0.14/
      DATA(CF(I,23, 7),I=1,5)/ 1168.0,0.,5.91E-10,1.6100,0.18/
      DATA(CF(I,23, 6),I=1,5)/ 1260.0,0.,5.02E-10,2.1200,0.15/
      DATA(CF(I,23, 5),I=1,5)/ 1355.0,1.,5.38E-11,0.1370,0.40/
      DATA(CF(I,23, 4),I=1,5)/ 1486.0,1.,5.56E-11,0.7080,0.10/
      DATA(CF(I,23, 3),I=1,5)/ 1571.0,0.,2.84E-11,0.0240,0.79/
      DATA(CF(I,23, 2),I=1,5)/ 6851.3,0.,2.54E-11,2.9200,0.09/
      DATA(CF(I,23, 1),I=1,5)/ 7246.1,0.,1.32E-11,3.5100,0.07/
      DATA(CF(I,24,24),I=1,5)/    6.8,1.,1.03E-07,0.2170,0.27/
      DATA(CF(I,24,23),I=1,5)/   16.5,0.,2.45E-07,0.3810,0.32/
      DATA(CF(I,24,22),I=1,5)/   31.0,0.,1.09E-07,0.5180,0.27/
      DATA(CF(I,24,21),I=1,5)/   49.1,1.,1.52E-08,0.1820,0.30/
      DATA(CF(I,24,20),I=1,5)/   69.5,0.,3.25E-08,1.3600,0.13/
      DATA(CF(I,24,19),I=1,5)/   90.6,1.,5.50E-09,0.1430,0.37/
      DATA(CF(I,24,18),I=1,5)/  160.2,1.,5.13E-09,0.6570,0.16/
      DATA(CF(I,24,17),I=1,5)/  184.7,1.,3.85E-09,0.7220,0.15/
      DATA(CF(I,24,16),I=1,5)/  209.3,1.,2.81E-09,0.7590,0.15/
      DATA(CF(I,24,15),I=1,5)/  244.4,1.,1.76E-09,0.7320,0.16/
      DATA(CF(I,24,14),I=1,5)/  271.0,1.,1.30E-09,0.7640,0.15/
      DATA(CF(I,24,13),I=1,5)/  298.0,1.,1.02E-09,0.8100,0.13/
      DATA(CF(I,24,12),I=1,5)/  354.8,1.,7.19E-10,0.8870,0.12/
      DATA(CF(I,24,11),I=1,5)/  384.2,1.,1.61E-10,0.1500,0.22/
      DATA(CF(I,24,10),I=1,5)/ 1011.0,1.,4.64E-10,0.9710,0.12/
      DATA(CF(I,24, 9),I=1,5)/ 1097.0,1.,3.31E-10,0.9240,0.14/
      DATA(CF(I,24, 8),I=1,5)/ 1185.0,1.,2.49E-10,0.9310,0.15/
      DATA(CF(I,24, 7),I=1,5)/ 1299.0,1.,1.68E-10,0.9100,0.14/
      DATA(CF(I,24, 6),I=1,5)/ 1396.0,1.,1.01E-10,0.8050,0.15/
      DATA(CF(I,24, 5),I=1,5)/ 1496.0,0.,1.17E-10,1.2100,0.35/
      DATA(CF(I,24, 4),I=1,5)/ 1634.0,1.,2.91E-11,0.8840,0.13/
      DATA(CF(I,24, 3),I=1,5)/ 1721.0,0.,1.45E-11,0.0350,0.80/
      DATA(CF(I,24, 2),I=1,5)/ 7482.0,1.,3.07E-12,0.9670,0.03/
      DATA(CF(I,24, 1),I=1,5)/ 7894.8,1.,1.46E-12,0.1830,0.39/
      DATA(CF(I,25,25),I=1,5)/    7.4,1.,8.56E-08,0.1320,0.26/
      DATA(CF(I,25,24),I=1,5)/   15.6,0.,1.18E-07,0.3590,0.19/
      DATA(CF(I,25,23),I=1,5)/   33.7,0.,8.54E-08,0.3970,0.32/
      DATA(CF(I,25,22),I=1,5)/   51.2,1.,1.80E-08,0.2720,0.18/
      DATA(CF(I,25,21),I=1,5)/   72.4,1.,8.22E-09,0.1610,0.32/
      DATA(CF(I,25,20),I=1,5)/   95.0,0.,2.15E-08,1.5400,0.11/
      DATA(CF(I,25,19),I=1,5)/  119.3,1.,3.65E-09,0.1470,0.37/
      DATA(CF(I,25,18),I=1,5)/  194.5,1.,3.91E-09,0.6990,0.15/
      DATA(CF(I,25,17),I=1,5)/  221.8,1.,2.92E-09,0.7190,0.15/
      DATA(CF(I,25,16),I=1,5)/  248.3,1.,2.23E-09,0.8060,0.14/
      DATA(CF(I,25,15),I=1,5)/  286.0,1.,1.39E-09,0.7350,0.16/
      DATA(CF(I,25,14),I=1,5)/  314.4,1.,1.04E-09,0.7610,0.15/
      DATA(CF(I,25,13),I=1,5)/  343.6,1.,8.28E-10,0.8090,0.13/
      DATA(CF(I,25,12),I=1,5)/  403.0,1.,5.60E-10,0.7870,0.14/
      DATA(CF(I,25,11),I=1,5)/  435.2,1.,1.52E-10,0.2990,0.08/
      DATA(CF(I,25,10),I=1,5)/ 1133.0,1.,4.03E-10,1.0400,0.11/
      DATA(CF(I,25, 9),I=1,5)/ 1244.0,1.,2.74E-10,0.9230,0.14/
      DATA(CF(I,25, 8),I=1,5)/ 1317.0,1.,2.18E-10,0.9900,0.14/
      DATA(CF(I,25, 7),I=1,5)/ 1437.0,1.,1.49E-10,0.9680,0.13/
      DATA(CF(I,25, 6),I=1,5)/ 1539.0,1.,8.70E-11,0.8020,0.15/
      DATA(CF(I,25, 5),I=1,5)/ 1644.0,0.,1.02E-10,1.2200,0.35/
      DATA(CF(I,25, 4),I=1,5)/ 1788.0,1.,2.54E-11,0.8830,0.13/
      DATA(CF(I,25, 3),I=1,5)/ 1880.0,0.,1.28E-11,0.0330,0.81/
      DATA(CF(I,25, 2),I=1,5)/ 8141.0,1.,2.77E-12,1.0100,0.02/
      DATA(CF(I,25, 1),I=1,5)/ 8571.9,1.,1.32E-12,0.2190,0.37/
      DATA(CF(I,26,26),I=1,5)/    7.9,0.,2.52E-07,0.7010,0.25/
      DATA(CF(I,26,25),I=1,5)/   16.2,1.,2.21E-08,0.0330,0.45/
      DATA(CF(I,26,24),I=1,5)/   30.6,0.,4.10E-08,0.3660,0.17/
      DATA(CF(I,26,23),I=1,5)/   54.8,0.,3.53E-08,0.2430,0.39/
      DATA(CF(I,26,22),I=1,5)/   75.0,1.,1.04E-08,0.2850,0.17/
      DATA(CF(I,26,21),I=1,5)/   99.0,1.,1.23E-08,0.4110,0.21/
      DATA(CF(I,26,20),I=1,5)/  125.0,1.,9.47E-09,0.4580,0.21/
      DATA(CF(I,26,19),I=1,5)/  151.1,1.,4.71E-09,0.2800,0.28/
      DATA(CF(I,26,18),I=1,5)/  233.6,1.,3.02E-09,0.6970,0.15/
      DATA(CF(I,26,17),I=1,5)/  262.1,1.,2.34E-09,0.7640,0.14/
      DATA(CF(I,26,16),I=1,5)/  290.0,1.,1.76E-09,0.8050,0.14/
      DATA(CF(I,26,15),I=1,5)/  331.0,1.,1.14E-09,0.7730,0.15/
      DATA(CF(I,26,14),I=1,5)/  361.0,1.,8.66E-10,0.8050,0.14/
      DATA(CF(I,26,13),I=1,5)/  392.0,1.,6.61E-10,0.7620,0.14/
      DATA(CF(I,26,12),I=1,5)/  457.0,1.,4.41E-10,0.6980,0.16/
      DATA(CF(I,26,11),I=1,5)/  489.3,1.,1.18E-10,0.2110,0.15/
      DATA(CF(I,26,10),I=1,5)/ 1262.0,1.,3.61E-10,1.1600,0.09/
      DATA(CF(I,26, 9),I=1,5)/ 1360.0,1.,2.45E-10,0.9780,0.13/
      DATA(CF(I,26, 8),I=1,5)/ 1470.0,1.,1.87E-10,0.9880,0.14/
      DATA(CF(I,26, 7),I=1,5)/ 1582.0,1.,1.33E-10,1.0300,0.12/
      DATA(CF(I,26, 6),I=1,5)/ 1690.0,1.,7.84E-11,0.8480,0.14/
      DATA(CF(I,26, 5),I=1,5)/ 1800.0,0.,8.90E-11,1.2000,0.35/
      DATA(CF(I,26, 4),I=1,5)/ 1960.0,1.,2.29E-11,0.9360,0.12/
      DATA(CF(I,26, 3),I=1,5)/ 2046.0,0.,1.12E-11,0.0340,0.81/
      DATA(CF(I,26, 2),I=1,5)/ 8828.0,1.,2.46E-12,1.0200,0.02/
      DATA(CF(I,26, 1),I=1,5)/ 9277.7,1.,9.79E-13,0.6640,0.14/
      DATA(CF(I,27,27),I=1,5)/    7.9,1.,8.89E-08,0.1270,0.24/
      DATA(CF(I,27,26),I=1,5)/   17.1,1.,5.65E-08,0.1940,0.23/
      DATA(CF(I,27,25),I=1,5)/   33.5,1.,3.06E-08,0.2010,0.22/
      DATA(CF(I,27,24),I=1,5)/   51.3,0.,2.27E-08,0.5740,0.10/
      DATA(CF(I,27,23),I=1,5)/   79.5,0.,1.93E-08,0.1950,0.42/
      DATA(CF(I,27,22),I=1,5)/  102.0,0.,1.27E-08,0.1260,0.47/
      DATA(CF(I,27,21),I=1,5)/  129.0,1.,3.58E-09,0.1940,0.29/
      DATA(CF(I,27,20),I=1,5)/  158.0,0.,1.17E-08,1.9800,0.07/
      DATA(CF(I,27,19),I=1,5)/  186.1,1.,1.78E-09,0.1120,0.42/
      DATA(CF(I,27,18),I=1,5)/  275.0,1.,2.41E-09,0.7390,0.14/
      DATA(CF(I,27,17),I=1,5)/  305.0,1.,1.86E-09,0.7620,0.14/
      DATA(CF(I,27,16),I=1,5)/  336.0,1.,1.41E-09,0.8040,0.14/
      DATA(CF(I,27,15),I=1,5)/  379.0,1.,9.54E-10,0.8130,0.14/
      DATA(CF(I,27,14),I=1,5)/  411.0,1.,7.12E-10,0.8030,0.14/
      DATA(CF(I,27,13),I=1,5)/  444.0,1.,5.34E-10,0.7180,0.15/
      DATA(CF(I,27,12),I=1,5)/  512.0,1.,3.62E-10,0.6580,0.17/
      DATA(CF(I,27,11),I=1,5)/  546.6,1.,1.05E-10,0.2520,0.12/
      DATA(CF(I,27,10),I=1,5)/ 1397.0,1.,3.10E-10,1.1700,0.09/
      DATA(CF(I,27, 9),I=1,5)/ 1486.0,1.,1.56E-10,0.5720,0.15/
      DATA(CF(I,27, 8),I=1,5)/ 1603.0,1.,1.32E-10,0.6820,0.13/
      DATA(CF(I,27, 7),I=1,5)/ 1735.0,1.,9.08E-11,0.5110,0.17/
      DATA(CF(I,27, 6),I=1,5)/ 1846.0,0.,3.45E-10,2.8400,0.11/
      DATA(CF(I,27, 5),I=1,5)/ 1962.0,1.,5.01E-11,0.7140,0.11/
      DATA(CF(I,27, 4),I=1,5)/ 2119.0,1.,1.92E-11,0.1170,0.42/
      DATA(CF(I,27, 3),I=1,5)/ 2219.0,1.,1.95E-11,1.2000,0.00/
      DATA(CF(I,27, 2),I=1,5)/ 9544.0,0.,1.68E-11,3.5200,0.05/
      DATA(CF(I,27, 1),I=1,5)/10012.1,1.,1.45E-12,0.6350,0.15/
      DATA(CF(I,28,28),I=1,5)/    7.6,0.,1.65E-07,0.4520,0.28/
      DATA(CF(I,28,27),I=1,5)/   18.2,0.,8.42E-08,0.6190,0.16/
      DATA(CF(I,28,26),I=1,5)/   35.3,1.,1.89E-08,0.2200,0.21/
      DATA(CF(I,28,25),I=1,5)/   54.9,1.,1.48E-08,0.2160,0.21/
      DATA(CF(I,28,24),I=1,5)/   76.0,0.,1.13E-08,0.5180,0.09/
      DATA(CF(I,28,23),I=1,5)/  108.0,0.,1.16E-08,0.1850,0.44/
      DATA(CF(I,28,22),I=1,5)/  133.0,0.,8.68E-09,0.1380,0.46/
      DATA(CF(I,28,21),I=1,5)/  162.0,1.,2.45E-09,0.1630,0.32/
      DATA(CF(I,28,20),I=1,5)/  193.0,0.,9.24E-09,2.2500,0.05/
      DATA(CF(I,28,19),I=1,5)/  225.0,0.,2.41E-09,0.0270,0.79/
      DATA(CF(I,28,18),I=1,5)/  321.0,1.,1.92E-09,0.7380,0.14/
      DATA(CF(I,28,17),I=1,5)/  352.0,1.,1.50E-09,0.7610,0.14/
      DATA(CF(I,28,16),I=1,5)/  384.0,1.,1.16E-09,0.8030,0.14/
      DATA(CF(I,28,15),I=1,5)/  430.0,1.,8.08E-10,0.8560,0.13/
      DATA(CF(I,28,14),I=1,5)/  464.0,1.,6.09E-10,0.8500,0.13/
      DATA(CF(I,28,13),I=1,5)/  499.0,1.,4.48E-10,0.7180,0.15/
      DATA(CF(I,28,12),I=1,5)/  571.3,1.,3.00E-10,0.6220,0.18/
      DATA(CF(I,28,11),I=1,5)/  607.0,1.,7.90E-11,0.1600,0.19/
      DATA(CF(I,28,10),I=1,5)/ 1541.0,1.,2.78E-10,1.2500,0.08/
      DATA(CF(I,28, 9),I=1,5)/ 1648.0,1.,1.92E-10,1.0400,0.12/
      DATA(CF(I,28, 8),I=1,5)/ 1756.0,1.,1.51E-10,1.1100,0.12/
      DATA(CF(I,28, 7),I=1,5)/ 1894.0,1.,1.05E-10,1.0900,0.11/
      DATA(CF(I,28, 6),I=1,5)/ 2011.0,1.,6.04E-11,0.8490,0.14/
      DATA(CF(I,28, 5),I=1,5)/ 2131.0,0.,6.91E-11,1.2100,0.35/
      DATA(CF(I,28, 4),I=1,5)/ 2295.0,1.,1.84E-11,0.9910,0.11/
      DATA(CF(I,28, 3),I=1,5)/ 2399.0,0.,9.03E-12,0.0420,0.79/
      DATA(CF(I,28, 2),I=1,5)/10290.0,1.,2.61E-12,0.5680,0.16/
      DATA(CF(I,28, 1),I=1,5)/10775.3,1.,1.39E-12,0.7290,0.13/
      end
!====================================================================
      subroutine rrfit_b(iz,in,t,r)
!*** Version 4. June 29, 1999.
!*** Written by D. A. Verner, verner@pa.uky.edu 
!******************************************************************************
!*** This subroutine calculates rates of radiative recombination for all ions
!*** of all elements from H through Zn by use of the following fits:
!*** H-like, He-like, Li-like, Na-like - Verner & Ferland, 1996, ApJS, 103, 467
!*** Other ions of C, N, O, Ne - Pequignot et al. 1991, A&A, 251, 680,
!***    refitted by Verner & Ferland formula to ensure correct asymptotes
!*** Fe XVII-XXIII - Arnaud & Raymond, 1992, ApJ, 398, 394
!*** Fe I-XV - refitted by Verner & Ferland formula to ensure correct asymptotes
!*** Other ions of Mg, Si, S, Ar, Ca, Fe, Ni - 
!***                      - Shull & Van Steenberg, 1982, ApJS, 48, 95
!*** Other ions of Na, Al - Landini & Monsignori Fossi, 1990, A&AS, 82, 229
!*** Other ions of F, P, Cl, K, Ti, Cr, Mn, Co (excluding Ti I-II, Cr I-IV,
!*** Mn I-V, Co I)        - Landini & Monsignori Fossi, 1991, A&AS, 91, 183
!*** All other species    - interpolations of the power-law fits
!*** Input parameters:  iz - atomic number 
!***                    in - number of electrons from 1 to iz 
!***                    t  - temperature, K
!*** Output parameter:  r  - rate coefficient, cm^3 s^(-1)
!******************************************************************************
      implicit none
      real :: rrec,rnew,fe,t,r,tt
      integer :: iz,in
      common/rrec/rrec(2,30,30)
      common/rnew/rnew(4,30,30)
      common/ferr/fe(3,13)
      r=0.0
      if(iz.lt.1.or.iz.gt.30)then
!         write(6,"(rrfit called with insane atomic number, =",iz")")
         stop
      endif
      if(in.lt.1.or.in.gt.iz)then
!         write(6,"(rrfit called with insane number elec ",in")")      
        stop
      endif
      if(in.le.3.or.in.eq.11.or.(iz.gt.5.and.iz.lt.9).or.iz.eq.10.or.(iz.eq.26.and.in.gt.11))then         
         tt=sqrt(t/rnew(3,iz,in))
         r=rnew(1,iz,in)/(tt*(tt+1.0)**(1.0-rnew(2,iz,in))*(1.0+sqrt(t/rnew(4,iz,in)))**(1.0+rnew(2,iz,in)))
      else
         tt=t*1.0e-04
         if(iz.eq.26.and.in.le.13)then
            r=fe(1,in)/tt**(fe(2,in)+fe(3,in)*log10(tt))
         else
            r=rrec(1,iz,in)/tt**rrec(2,iz,in)
         endif
      endif
      return 
      end subroutine rrfit_b
!*********************************************
      block data radrec
      common/rrec/rrec(2,30,30)
      common/rnew/rnew(4,30,30)
      common/ferr/fe(3,13)
      data(rrec(i, 4, 4),i=1,2)/4.500e-13,0.6480/
      data(rrec(i, 5, 4),i=1,2)/2.680e-12,0.7250/
      data(rrec(i, 6, 4),i=1,2)/4.900e-12,0.8030/
      data(rrec(i, 7, 4),i=1,2)/9.400e-12,0.7650/
      data(rrec(i, 8, 4),i=1,2)/1.590e-11,0.7590/
      data(rrec(i, 9, 4),i=1,2)/1.919e-11,0.7849/
      data(rrec(i,10, 4),i=1,2)/2.800e-11,0.7710/
      data(rrec(i,11, 4),i=1,2)/4.030e-11,0.7873/
      data(rrec(i,12, 4),i=1,2)/6.830e-11,0.7650/
      data(rrec(i,13, 4),i=1,2)/8.343e-11,0.8291/
      data(rrec(i,14, 4),i=1,2)/1.430e-10,0.8230/
      data(rrec(i,15, 4),i=1,2)/1.460e-10,0.8391/
      data(rrec(i,16, 4),i=1,2)/2.000e-10,0.8060/
      data(rrec(i,17, 4),i=1,2)/2.091e-10,0.8691/
      data(rrec(i,18, 4),i=1,2)/3.070e-10,0.8190/
      data(rrec(i,19, 4),i=1,2)/3.210e-10,0.8934/
      data(rrec(i,20, 4),i=1,2)/4.610e-10,0.8330/
      data(rrec(i,21, 4),i=1,2)/4.870e-10,0.8060/
      data(rrec(i,22, 4),i=1,2)/5.139e-10,0.7781/
      data(rrec(i,23, 4),i=1,2)/5.850e-10,0.7570/
      data(rrec(i,24, 4),i=1,2)/6.556e-10,0.7359/
      data(rrec(i,25, 4),i=1,2)/7.238e-10,0.7242/
      data(rrec(i,27, 4),i=1,2)/8.404e-10,0.7167/
      data(rrec(i,28, 4),i=1,2)/1.360e-09,0.8420/
      data(rrec(i,29, 4),i=1,2)/1.880e-09,0.8420/
      data(rrec(i,30, 4),i=1,2)/2.400e-09,0.8420/
      data(rrec(i, 5, 5),i=1,2)/4.600e-13,0.6360/
      data(rrec(i, 6, 5),i=1,2)/2.300e-12,0.6450/
      data(rrec(i, 7, 5),i=1,2)/5.000e-12,0.6760/
      data(rrec(i, 8, 5),i=1,2)/9.600e-12,0.6700/
      data(rrec(i, 9, 5),i=1,2)/1.558e-11,0.6816/
      data(rrec(i,10, 5),i=1,2)/2.300e-11,0.7040/
      data(rrec(i,11, 5),i=1,2)/3.253e-11,0.7075/
      data(rrec(i,12, 5),i=1,2)/4.600e-11,0.7110/
      data(rrec(i,13, 5),i=1,2)/5.951e-11,0.7125/
      data(rrec(i,14, 5),i=1,2)/7.700e-11,0.7140/
      data(rrec(i,15, 5),i=1,2)/1.042e-10,0.7330/
      data(rrec(i,16, 5),i=1,2)/1.400e-10,0.7550/
      data(rrec(i,17, 5),i=1,2)/1.760e-10,0.7682/
      data(rrec(i,18, 5),i=1,2)/2.140e-10,0.7740/
      data(rrec(i,19, 5),i=1,2)/2.629e-10,0.7772/
      data(rrec(i,20, 5),i=1,2)/3.240e-10,0.7800/
      data(rrec(i,21, 5),i=1,2)/3.970e-10,0.7840/
      data(rrec(i,22, 5),i=1,2)/4.700e-10,0.7873/
      data(rrec(i,23, 5),i=1,2)/5.500e-10,0.7910/
      data(rrec(i,24, 5),i=1,2)/6.301e-10,0.7952/
      data(rrec(i,25, 5),i=1,2)/7.058e-10,0.7982/
      data(rrec(i,27, 5),i=1,2)/8.252e-10,0.8004/
      data(rrec(i,28, 5),i=1,2)/8.710e-10,0.8000/
      data(rrec(i,29, 5),i=1,2)/9.170e-10,0.8000/
      data(rrec(i,30, 5),i=1,2)/9.630e-10,0.7990/
      data(rrec(i, 6, 6),i=1,2)/4.700e-13,0.6240/
      data(rrec(i, 7, 6),i=1,2)/2.200e-12,0.6390/
      data(rrec(i, 8, 6),i=1,2)/5.100e-12,0.6600/
      data(rrec(i, 9, 6),i=1,2)/9.171e-12,0.6757/
      data(rrec(i,10, 6),i=1,2)/1.500e-11,0.6840/
      data(rrec(i,11, 6),i=1,2)/2.191e-11,0.6875/
      data(rrec(i,12, 6),i=1,2)/3.200e-11,0.6910/
      data(rrec(i,13, 6),i=1,2)/4.308e-11,0.6970/
      data(rrec(i,14, 6),i=1,2)/5.800e-11,0.7030/
      data(rrec(i,15, 6),i=1,2)/7.316e-11,0.7027/
      data(rrec(i,16, 6),i=1,2)/9.200e-11,0.7140/
      data(rrec(i,17, 6),i=1,2)/1.198e-10,0.7508/
      data(rrec(i,18, 6),i=1,2)/1.580e-10,0.7900/
      data(rrec(i,19, 6),i=1,2)/2.048e-10,0.8032/
      data(rrec(i,20, 6),i=1,2)/2.600e-10,0.8000/
      data(rrec(i,21, 6),i=1,2)/3.280e-10,0.7990/
      data(rrec(i,22, 6),i=1,2)/3.966e-10,0.7973/
      data(rrec(i,23, 6),i=1,2)/4.760e-10,0.8000/
      data(rrec(i,24, 6),i=1,2)/5.547e-10,0.8027/
      data(rrec(i,25, 6),i=1,2)/6.313e-10,0.8058/
      data(rrec(i,27, 6),i=1,2)/7.503e-10,0.8085/
      data(rrec(i,28, 6),i=1,2)/7.940e-10,0.8080/
      data(rrec(i,29, 6),i=1,2)/8.380e-10,0.8080/
      data(rrec(i,30, 6),i=1,2)/8.810e-10,0.8070/
      data(rrec(i, 7, 7),i=1,2)/4.100e-13,0.6080/
      data(rrec(i, 8, 7),i=1,2)/2.000e-12,0.6460/
      data(rrec(i, 9, 7),i=1,2)/5.231e-12,0.6615/
      data(rrec(i,10, 7),i=1,2)/9.100e-12,0.6680/
      data(rrec(i,11, 7),i=1,2)/1.447e-11,0.6814/
      data(rrec(i,12, 7),i=1,2)/2.300e-11,0.6950/
      data(rrec(i,13, 7),i=1,2)/3.145e-11,0.6915/
      data(rrec(i,14, 7),i=1,2)/4.300e-11,0.6880/
      data(rrec(i,15, 7),i=1,2)/5.659e-11,0.7023/
      data(rrec(i,16, 7),i=1,2)/7.400e-11,0.7160/
      data(rrec(i,17, 7),i=1,2)/9.561e-11,0.7102/
      data(rrec(i,18, 7),i=1,2)/1.230e-10,0.7020/
      data(rrec(i,19, 7),i=1,2)/1.587e-10,0.7105/
      data(rrec(i,20, 7),i=1,2)/2.040e-10,0.7300/
      data(rrec(i,21, 7),i=1,2)/2.630e-10,0.7490/
      data(rrec(i,22, 7),i=1,2)/3.220e-10,0.7683/
      data(rrec(i,23, 7),i=1,2)/3.950e-10,0.7830/
      data(rrec(i,24, 7),i=1,2)/4.671e-10,0.7967/
      data(rrec(i,25, 7),i=1,2)/5.407e-10,0.8058/
      data(rrec(i,27, 7),i=1,2)/6.611e-10,0.8121/
      data(rrec(i,28, 7),i=1,2)/7.080e-10,0.8110/
      data(rrec(i,29, 7),i=1,2)/7.550e-10,0.8100/
      data(rrec(i,30, 7),i=1,2)/8.020e-10,0.8090/
      data(rrec(i, 8, 8),i=1,2)/3.100e-13,0.6780/
      data(rrec(i, 9, 8),i=1,2)/1.344e-12,0.6708/
      data(rrec(i,10, 8),i=1,2)/4.400e-12,0.6750/
      data(rrec(i,11, 8),i=1,2)/7.849e-12,0.6952/
      data(rrec(i,12, 8),i=1,2)/1.400e-11,0.7160/
      data(rrec(i,13, 8),i=1,2)/2.049e-11,0.7090/
      data(rrec(i,14, 8),i=1,2)/3.000e-11,0.7020/
      data(rrec(i,15, 8),i=1,2)/4.125e-11,0.6965/
      data(rrec(i,16, 8),i=1,2)/5.500e-11,0.7110/
      data(rrec(i,17, 8),i=1,2)/7.280e-11,0.7518/
      data(rrec(i,18, 8),i=1,2)/9.550e-11,0.7930/
      data(rrec(i,19, 8),i=1,2)/1.235e-10,0.8052/
      data(rrec(i,20, 8),i=1,2)/1.580e-10,0.8000/
      data(rrec(i,21, 8),i=1,2)/2.060e-10,0.7990/
      data(rrec(i,22, 8),i=1,2)/2.537e-10,0.7977/
      data(rrec(i,23, 8),i=1,2)/3.190e-10,0.8020/
      data(rrec(i,24, 8),i=1,2)/3.844e-10,0.8071/
      data(rrec(i,25, 8),i=1,2)/4.564e-10,0.8124/
      data(rrec(i,27, 8),i=1,2)/5.842e-10,0.8168/
      data(rrec(i,28, 8),i=1,2)/6.380e-10,0.8160/
      data(rrec(i,29, 8),i=1,2)/6.920e-10,0.8150/
      data(rrec(i,30, 8),i=1,2)/7.460e-10,0.8140/
      data(rrec(i, 9, 9),i=1,2)/6.273e-13,0.6798/
      data(rrec(i,10, 9),i=1,2)/1.500e-12,0.6930/
      data(rrec(i,11, 9),i=1,2)/3.399e-12,0.7054/
      data(rrec(i,12, 9),i=1,2)/7.700e-12,0.7180/
      data(rrec(i,13, 9),i=1,2)/1.275e-11,0.7170/
      data(rrec(i,14, 9),i=1,2)/2.110e-11,0.7160/
      data(rrec(i,15, 9),i=1,2)/2.975e-11,0.6945/
      data(rrec(i,16, 9),i=1,2)/4.000e-11,0.6960/
      data(rrec(i,17, 9),i=1,2)/5.281e-11,0.7491/
      data(rrec(i,18, 9),i=1,2)/6.920e-11,0.8110/
      data(rrec(i,19, 9),i=1,2)/9.044e-11,0.8251/
      data(rrec(i,20, 9),i=1,2)/1.180e-10,0.8100/
      data(rrec(i,21, 9),i=1,2)/1.580e-10,0.8040/
      data(rrec(i,22, 9),i=1,2)/1.983e-10,0.7980/
      data(rrec(i,23, 9),i=1,2)/2.570e-10,0.8040/
      data(rrec(i,24, 9),i=1,2)/3.154e-10,0.8101/
      data(rrec(i,25, 9),i=1,2)/3.837e-10,0.8183/
      data(rrec(i,27, 9),i=1,2)/5.147e-10,0.8253/
      data(rrec(i,28, 9),i=1,2)/5.750e-10,0.8240/
      data(rrec(i,29, 9),i=1,2)/6.350e-10,0.8230/
      data(rrec(i,30, 9),i=1,2)/6.960e-10,0.8210/
      data(rrec(i,10,10),i=1,2)/2.200e-13,0.7590/
      data(rrec(i,11,10),i=1,2)/8.775e-13,0.7467/
      data(rrec(i,12,10),i=1,2)/3.500e-12,0.7340/
      data(rrec(i,13,10),i=1,2)/6.481e-12,0.7345/
      data(rrec(i,14,10),i=1,2)/1.200e-11,0.7350/
      data(rrec(i,15,10),i=1,2)/1.834e-11,0.7285/
      data(rrec(i,16,10),i=1,2)/2.700e-11,0.7330/
      data(rrec(i,17,10),i=1,2)/3.711e-11,0.7641/
      data(rrec(i,18,10),i=1,2)/4.900e-11,0.8010/
      data(rrec(i,19,10),i=1,2)/6.444e-11,0.8175/
      data(rrec(i,20,10),i=1,2)/8.510e-11,0.8200/
      data(rrec(i,21,10),i=1,2)/1.170e-10,0.8220/
      data(rrec(i,22,10),i=1,2)/1.494e-10,0.8242/
      data(rrec(i,23,10),i=1,2)/2.010e-10,0.8280/
      data(rrec(i,24,10),i=1,2)/2.525e-10,0.8311/
      data(rrec(i,25,10),i=1,2)/3.177e-10,0.8341/
      data(rrec(i,27,10),i=1,2)/4.552e-10,0.8364/
      data(rrec(i,28,10),i=1,2)/5.250e-10,0.8360/
      data(rrec(i,29,10),i=1,2)/5.950e-10,0.8360/
      data(rrec(i,30,10),i=1,2)/6.650e-10,0.8350/
      data(rrec(i,12,12),i=1,2)/1.400e-13,0.8550/
      data(rrec(i,13,12),i=1,2)/7.197e-13,0.7697/
      data(rrec(i,14,12),i=1,2)/3.700e-12,0.6930/
      data(rrec(i,15,12),i=1,2)/7.980e-12,0.6829/
      data(rrec(i,16,12),i=1,2)/1.200e-11,0.7010/
      data(rrec(i,17,12),i=1,2)/1.800e-11,0.7232/
      data(rrec(i,18,12),i=1,2)/2.690e-11,0.7440/
      data(rrec(i,19,12),i=1,2)/3.748e-11,0.7628/
      data(rrec(i,20,12),i=1,2)/5.040e-11,0.7800/
      data(rrec(i,21,12),i=1,2)/7.240e-11,0.7950/
      data(rrec(i,22,12),i=1,2)/9.440e-11,0.8107/
      data(rrec(i,23,12),i=1,2)/1.350e-10,0.8220/
      data(rrec(i,24,12),i=1,2)/1.751e-10,0.8340/
      data(rrec(i,25,12),i=1,2)/2.298e-10,0.8417/
      data(rrec(i,27,12),i=1,2)/3.461e-10,0.8469/
      data(rrec(i,28,12),i=1,2)/4.030e-10,0.8460/
      data(rrec(i,29,12),i=1,2)/4.600e-10,0.8450/
      data(rrec(i,30,12),i=1,2)/5.170e-10,0.8440/
      data(rrec(i,13,13),i=1,2)/3.980e-13,0.8019/
      data(rrec(i,14,13),i=1,2)/1.000e-12,0.7860/
      data(rrec(i,15,13),i=1,2)/2.558e-12,0.7629/
      data(rrec(i,16,13),i=1,2)/5.700e-12,0.7550/
      data(rrec(i,17,13),i=1,2)/1.011e-11,0.7703/
      data(rrec(i,18,13),i=1,2)/1.580e-11,0.7930/
      data(rrec(i,19,13),i=1,2)/2.448e-11,0.8052/
      data(rrec(i,20,13),i=1,2)/3.760e-11,0.8100/
      data(rrec(i,21,13),i=1,2)/5.870e-11,0.8150/
      data(rrec(i,22,13),i=1,2)/7.972e-11,0.8206/
      data(rrec(i,23,13),i=1,2)/1.130e-10,0.8270/
      data(rrec(i,24,13),i=1,2)/1.470e-10,0.8325/
      data(rrec(i,25,13),i=1,2)/1.908e-10,0.8372/
      data(rrec(i,27,13),i=1,2)/2.976e-10,0.8406/
      data(rrec(i,28,13),i=1,2)/3.630e-10,0.8400/
      data(rrec(i,29,13),i=1,2)/4.280e-10,0.8390/
      data(rrec(i,30,13),i=1,2)/4.940e-10,0.8390/
      data(rrec(i,14,14),i=1,2)/5.900e-13,0.6010/
      data(rrec(i,15,14),i=1,2)/1.294e-12,0.6766/
      data(rrec(i,16,14),i=1,2)/2.700e-12,0.7450/
      data(rrec(i,17,14),i=1,2)/5.165e-12,0.7893/
      data(rrec(i,18,14),i=1,2)/9.120e-12,0.8110/
      data(rrec(i,19,14),i=1,2)/1.513e-11,0.8186/
      data(rrec(i,20,14),i=1,2)/2.400e-11,0.8200/
      data(rrec(i,21,14),i=1,2)/3.960e-11,0.8220/
      data(rrec(i,22,14),i=1,2)/5.518e-11,0.8245/
      data(rrec(i,23,14),i=1,2)/8.370e-11,0.8280/
      data(rrec(i,24,14),i=1,2)/1.123e-10,0.8313/
      data(rrec(i,25,14),i=1,2)/1.525e-10,0.8342/
      data(rrec(i,26,14),i=1,2)/2.000e-10,0.8360/
      data(rrec(i,27,14),i=1,2)/2.537e-10,0.8364/
      data(rrec(i,28,14),i=1,2)/3.160e-10,0.8360/
      data(rrec(i,29,14),i=1,2)/3.780e-10,0.8360/
      data(rrec(i,30,14),i=1,2)/4.410e-10,0.8350/
      data(rrec(i,15,15),i=1,2)/9.761e-13,0.6209/
      data(rrec(i,16,15),i=1,2)/1.800e-12,0.6860/
      data(rrec(i,17,15),i=1,2)/3.320e-12,0.7579/
      data(rrec(i,18,15),i=1,2)/6.030e-12,0.8120/
      data(rrec(i,19,15),i=1,2)/1.063e-11,0.8269/
      data(rrec(i,20,15),i=1,2)/1.800e-11,0.8200/
      data(rrec(i,21,15),i=1,2)/3.130e-11,0.8180/
      data(rrec(i,22,15),i=1,2)/4.451e-11,0.8153/
      data(rrec(i,23,15),i=1,2)/6.830e-11,0.8200/
      data(rrec(i,24,15),i=1,2)/9.206e-11,0.8246/
      data(rrec(i,25,15),i=1,2)/1.250e-10,0.8302/
      data(rrec(i,26,15),i=1,2)/1.640e-10,0.8340/
      data(rrec(i,27,15),i=1,2)/2.093e-10,0.8349/
      data(rrec(i,28,15),i=1,2)/2.630e-10,0.8340/
      data(rrec(i,29,15),i=1,2)/3.170e-10,0.8330/
      data(rrec(i,30,15),i=1,2)/3.700e-10,0.8320/
      data(rrec(i,16,16),i=1,2)/4.100e-13,0.6300/
      data(rrec(i,17,16),i=1,2)/1.248e-12,0.7663/
      data(rrec(i,18,16),i=1,2)/3.230e-12,0.8690/
      data(rrec(i,19,16),i=1,2)/6.384e-12,0.8790/
      data(rrec(i,20,16),i=1,2)/1.070e-11,0.8400/
      data(rrec(i,21,16),i=1,2)/1.920e-11,0.8210/
      data(rrec(i,22,16),i=1,2)/2.765e-11,0.8012/
      data(rrec(i,23,16),i=1,2)/4.650e-11,0.8060/
      data(rrec(i,24,16),i=1,2)/6.539e-11,0.8099/
      data(rrec(i,25,16),i=1,2)/9.539e-11,0.8202/
      data(rrec(i,26,16),i=1,2)/1.330e-10,0.8280/
      data(rrec(i,27,16),i=1,2)/1.769e-10,0.8299/
      data(rrec(i,28,16),i=1,2)/2.290e-10,0.8280/
      data(rrec(i,29,16),i=1,2)/2.810e-10,0.8260/
      data(rrec(i,30,16),i=1,2)/3.330e-10,0.8240/
      data(rrec(i,17,17),i=1,2)/1.010e-12,0.7380/
      data(rrec(i,18,17),i=1,2)/1.950e-12,0.7520/
      data(rrec(i,19,17),i=1,2)/3.766e-12,0.7662/
      data(rrec(i,20,17),i=1,2)/7.080e-12,0.7800/
      data(rrec(i,21,17),i=1,2)/1.430e-11,0.7920/
      data(rrec(i,22,17),i=1,2)/2.152e-11,0.8038/
      data(rrec(i,23,17),i=1,2)/3.740e-11,0.8120/
      data(rrec(i,24,17),i=1,2)/5.335e-11,0.8207/
      data(rrec(i,25,17),i=1,2)/7.807e-11,0.8260/
      data(rrec(i,26,17),i=1,2)/1.090e-10,0.8290/
      data(rrec(i,27,17),i=1,2)/1.459e-10,0.8296/
      data(rrec(i,28,17),i=1,2)/1.910e-10,0.8290/
      data(rrec(i,29,17),i=1,2)/2.360e-10,0.8280/
      data(rrec(i,30,17),i=1,2)/2.810e-10,0.8280/
      data(rrec(i,18,18),i=1,2)/3.770e-13,0.6510/
      data(rrec(i,19,18),i=1,2)/1.304e-12,0.6753/
      data(rrec(i,20,18),i=1,2)/3.960e-12,0.7000/
      data(rrec(i,21,18),i=1,2)/1.130e-11,0.7240/
      data(rrec(i,22,18),i=1,2)/1.857e-11,0.7484/
      data(rrec(i,23,18),i=1,2)/3.170e-11,0.7680/
      data(rrec(i,24,18),i=1,2)/4.479e-11,0.7883/
      data(rrec(i,25,18),i=1,2)/6.106e-11,0.8020/
      data(rrec(i,26,18),i=1,2)/8.130e-11,0.8100/
      data(rrec(i,27,18),i=1,2)/1.098e-10,0.8118/
      data(rrec(i,28,18),i=1,2)/1.500e-10,0.8100/
      data(rrec(i,29,18),i=1,2)/1.900e-10,0.8080/
      data(rrec(i,30,18),i=1,2)/2.300e-10,0.8060/
      data(rrec(i,19,19),i=1,2)/2.762e-13,0.8023/
      data(rrec(i,20,19),i=1,2)/6.780e-13,0.8000/
      data(rrec(i,21,19),i=1,2)/2.330e-12,0.7980/
      data(rrec(i,22,19),i=1,2)/3.983e-12,0.7955/
      data(rrec(i,23,19),i=1,2)/1.150e-11,0.7940/
      data(rrec(i,24,19),i=1,2)/1.906e-11,0.7919/
      data(rrec(i,25,19),i=1,2)/3.620e-11,0.7907/
      data(rrec(i,26,19),i=1,2)/6.050e-11,0.7900/
      data(rrec(i,27,19),i=1,2)/8.818e-11,0.7898/
      data(rrec(i,28,19),i=1,2)/1.190e-10,0.7900/
      data(rrec(i,29,19),i=1,2)/1.500e-10,0.7900/
      data(rrec(i,30,19),i=1,2)/1.810e-10,0.7900/
      data(rrec(i,20,20),i=1,2)/1.120e-13,0.9000/
      data(rrec(i,21,20),i=1,2)/6.540e-13,0.8670/
      data(rrec(i,22,20),i=1,2)/1.196e-12,0.8344/
      data(rrec(i,23,20),i=1,2)/5.330e-12,0.8090/
      data(rrec(i,24,20),i=1,2)/9.471e-12,0.7846/
      data(rrec(i,25,20),i=1,2)/2.169e-11,0.7683/
      data(rrec(i,26,20),i=1,2)/4.120e-11,0.7590/
      data(rrec(i,27,20),i=1,2)/6.409e-11,0.7570/
      data(rrec(i,28,20),i=1,2)/8.910e-11,0.7590/
      data(rrec(i,29,20),i=1,2)/1.140e-10,0.7610/
      data(rrec(i,30,20),i=1,2)/1.390e-10,0.7630/
      data(rrec(i,21,21),i=1,2)/1.170e-13,0.8980/
      data(rrec(i,22,21),i=1,2)/5.330e-12,0.8640/
      data(rrec(i,23,21),i=1,2)/1.060e-11,0.8300/
      data(rrec(i,24,21),i=1,2)/1.580e-11,0.7960/
      data(rrec(i,25,21),i=1,2)/2.100e-11,0.7620/
      data(rrec(i,26,21),i=1,2)/2.620e-11,0.7280/
      data(rrec(i,27,21),i=1,2)/2.822e-11,0.7280/
      data(rrec(i,28,21),i=1,2)/3.040e-11,0.7280/
      data(rrec(i,29,21),i=1,2)/3.260e-11,0.7280/
      data(rrec(i,30,21),i=1,2)/3.480e-11,0.7280/
      data(rrec(i,22,22),i=1,2)/1.220e-13,0.8970/
      data(rrec(i,23,22),i=1,2)/3.870e-12,0.8480/
      data(rrec(i,24,22),i=1,2)/7.610e-12,0.7980/
      data(rrec(i,25,22),i=1,2)/1.140e-11,0.7480/
      data(rrec(i,26,22),i=1,2)/1.510e-11,0.6990/
      data(rrec(i,27,22),i=1,2)/1.626e-11,0.6990/
      data(rrec(i,28,22),i=1,2)/1.750e-11,0.6990/
      data(rrec(i,29,22),i=1,2)/1.870e-11,0.6990/
      data(rrec(i,30,22),i=1,2)/2.000e-11,0.6990/
      data(rrec(i,23,23),i=1,2)/1.270e-13,0.8950/
      data(rrec(i,24,23),i=1,2)/2.680e-12,0.8240/
      data(rrec(i,25,23),i=1,2)/5.240e-12,0.7530/
      data(rrec(i,26,23),i=1,2)/7.800e-12,0.6820/
      data(rrec(i,27,23),i=1,2)/8.402e-12,0.6820/
      data(rrec(i,28,23),i=1,2)/9.050e-12,0.6820/
      data(rrec(i,29,23),i=1,2)/9.700e-12,0.6820/
      data(rrec(i,30,23),i=1,2)/1.030e-11,0.6820/
      data(rrec(i,24,24),i=1,2)/1.320e-13,0.8940/
      data(rrec(i,25,24),i=1,2)/1.730e-12,0.8200/
      data(rrec(i,26,24),i=1,2)/3.320e-12,0.7460/
      data(rrec(i,27,24),i=1,2)/3.575e-12,0.7460/
      data(rrec(i,28,24),i=1,2)/3.580e-12,0.7460/
      data(rrec(i,29,24),i=1,2)/3.580e-12,0.7460/
      data(rrec(i,30,24),i=1,2)/3.590e-12,0.7460/
      data(rrec(i,25,25),i=1,2)/1.370e-13,0.8920/
      data(rrec(i,26,25),i=1,2)/1.020e-12,0.8430/
      data(rrec(i,27,25),i=1,2)/1.278e-12,0.7682/
      data(rrec(i,28,25),i=1,2)/1.600e-12,0.7000/
      data(rrec(i,29,25),i=1,2)/1.920e-12,0.7000/
      data(rrec(i,30,25),i=1,2)/2.240e-12,0.7000/
      data(rrec(i,26,26),i=1,2)/1.420e-13,0.8910/
      data(rrec(i,27,26),i=1,2)/4.459e-13,0.7897/
      data(rrec(i,28,26),i=1,2)/1.400e-12,0.7000/
      data(rrec(i,29,26),i=1,2)/1.500e-12,0.7000/
      data(rrec(i,30,26),i=1,2)/1.600e-12,0.7000/
      data(rrec(i,27,27),i=1,2)/2.510e-13,0.7950/
      data(rrec(i,28,27),i=1,2)/1.000e-12,0.7000/
      data(rrec(i,29,27),i=1,2)/1.000e-12,0.7000/
      data(rrec(i,30,27),i=1,2)/1.100e-12,0.7000/
      data(rrec(i,28,28),i=1,2)/3.600e-13,0.7000/
      data(rrec(i,29,28),i=1,2)/3.600e-13,0.7000/
      data(rrec(i,30,28),i=1,2)/3.600e-13,0.7000/
      data(rrec(i,29,29),i=1,2)/3.600e-13,0.7000/
      data(rrec(i,30,29),i=1,2)/3.600e-13,0.7000/
      data(rrec(i,30,30),i=1,2)/3.600e-13,0.7000/
      data(rnew(i, 1, 1),i=1,4)/7.982e-11,0.7480,3.148e+00,7.036e+05/
      data(rnew(i, 2, 1),i=1,4)/1.891e-10,0.7524,9.370e+00,2.774e+06/
      data(rnew(i, 3, 1),i=1,4)/3.039e-10,0.7539,1.871e+01,6.209e+06/
      data(rnew(i, 4, 1),i=1,4)/4.290e-10,0.7557,3.000e+01,1.093e+07/
      data(rnew(i, 5, 1),i=1,4)/5.437e-10,0.7560,4.576e+01,1.706e+07/
      data(rnew(i, 6, 1),i=1,4)/6.556e-10,0.7567,6.523e+01,2.446e+07/
      data(rnew(i, 7, 1),i=1,4)/7.586e-10,0.7563,9.015e+01,3.338e+07/
      data(rnew(i, 8, 1),i=1,4)/8.616e-10,0.7563,1.191e+02,4.352e+07/
      data(rnew(i, 9, 1),i=1,4)/9.712e-10,0.7566,1.499e+02,5.498e+07/
      data(rnew(i,10, 1),i=1,4)/1.085e-09,0.7570,1.834e+02,6.776e+07/
      data(rnew(i,11, 1),i=1,4)/1.163e-09,0.7558,2.328e+02,8.262e+07/
      data(rnew(i,12, 1),i=1,4)/1.317e-09,0.7574,2.585e+02,9.769e+07/
      data(rnew(i,13, 1),i=1,4)/1.419e-09,0.7578,3.057e+02,1.143e+08/
      data(rnew(i,14, 1),i=1,4)/1.517e-09,0.7574,3.601e+02,1.329e+08/
      data(rnew(i,15, 1),i=1,4)/1.586e-09,0.7560,4.327e+02,1.534e+08/
      data(rnew(i,16, 1),i=1,4)/1.729e-09,0.7568,4.725e+02,1.746e+08/
      data(rnew(i,17, 1),i=1,4)/1.791e-09,0.7565,5.591e+02,1.972e+08/
      data(rnew(i,18, 1),i=1,4)/1.913e-09,0.7567,6.175e+02,2.212e+08/
      data(rnew(i,19, 1),i=1,4)/2.033e-09,0.7569,6.797e+02,2.463e+08/
      data(rnew(i,20, 1),i=1,4)/2.129e-09,0.7570,7.591e+02,2.739e+08/
      data(rnew(i,21, 1),i=1,4)/2.262e-09,0.7578,8.186e+02,3.000e+08/
      data(rnew(i,22, 1),i=1,4)/2.370e-09,0.7574,9.002e+02,3.307e+08/
      data(rnew(i,23, 1),i=1,4)/2.415e-09,0.7565,1.032e+03,3.635e+08/
      data(rnew(i,24, 1),i=1,4)/2.537e-09,0.7571,1.108e+03,3.954e+08/
      data(rnew(i,25, 1),i=1,4)/2.618e-09,0.7565,1.225e+03,4.307e+08/
      data(rnew(i,26, 1),i=1,4)/2.735e-09,0.7568,1.314e+03,4.659e+08/
      data(rnew(i,27, 1),i=1,4)/2.809e-09,0.7565,1.444e+03,5.042e+08/
      data(rnew(i,28, 1),i=1,4)/3.002e-09,0.7581,1.467e+03,5.409e+08/
      data(rnew(i,29, 1),i=1,4)/3.022e-09,0.7564,1.666e+03,5.855e+08/
      data(rnew(i,30, 1),i=1,4)/3.127e-09,0.7567,1.779e+03,6.246e+08/
      data(rnew(i, 2, 2),i=1,4)/9.356e-10,0.7892,4.266e-02,4.677e+06/
      data(rnew(i, 3, 2),i=1,4)/1.112e-10,0.6926,2.437e+01,8.323e+06/
      data(rnew(i, 4, 2),i=1,4)/1.317e-10,0.6691,8.473e+01,1.412e+07/
      data(rnew(i, 5, 2),i=1,4)/1.922e-10,0.6717,1.272e+02,1.975e+07/
      data(rnew(i, 6, 2),i=1,4)/2.765e-10,0.6858,1.535e+02,2.556e+07/
      data(rnew(i, 7, 2),i=1,4)/3.910e-10,0.6988,1.611e+02,3.271e+07/
      data(rnew(i, 8, 2),i=1,4)/4.897e-10,0.7048,1.906e+02,4.093e+07/
      data(rnew(i, 9, 2),i=1,4)/5.602e-10,0.7052,2.476e+02,5.077e+07/
      data(rnew(i,10, 2),i=1,4)/6.161e-10,0.7029,3.274e+02,6.243e+07/
      data(rnew(i,11, 2),i=1,4)/6.833e-10,0.7018,4.060e+02,7.491e+07/
      data(rnew(i,12, 2),i=1,4)/7.510e-10,0.7020,4.921e+02,8.643e+07/
      data(rnew(i,13, 2),i=1,4)/8.182e-10,0.7008,5.875e+02,1.007e+08/
      data(rnew(i,14, 2),i=1,4)/8.722e-10,0.6996,7.098e+02,1.155e+08/
      data(rnew(i,15, 2),i=1,4)/9.142e-10,0.6961,8.682e+02,1.335e+08/
      data(rnew(i,16, 2),i=1,4)/9.692e-10,0.6945,1.017e+03,1.517e+08/
      data(rnew(i,17, 2),i=1,4)/1.021e-09,0.6932,1.184e+03,1.695e+08/
      data(rnew(i,18, 2),i=1,4)/1.087e-09,0.6936,1.329e+03,1.880e+08/
      data(rnew(i,19, 2),i=1,4)/1.145e-09,0.6921,1.503e+03,2.098e+08/
      data(rnew(i,20, 2),i=1,4)/1.179e-09,0.6893,1.757e+03,2.344e+08/
      data(rnew(i,21, 2),i=1,4)/1.265e-09,0.6902,1.877e+03,2.555e+08/
      data(rnew(i,22, 2),i=1,4)/1.322e-09,0.6885,2.092e+03,2.829e+08/
      data(rnew(i,23, 2),i=1,4)/1.375e-09,0.6885,2.321e+03,3.056e+08/
      data(rnew(i,24, 2),i=1,4)/1.422e-09,0.6874,2.589e+03,3.336e+08/
      data(rnew(i,25, 2),i=1,4)/1.488e-09,0.6867,2.802e+03,3.623e+08/
      data(rnew(i,26, 2),i=1,4)/1.542e-09,0.6859,3.073e+03,3.926e+08/
      data(rnew(i,27, 2),i=1,4)/1.589e-09,0.6846,3.373e+03,4.267e+08/
      data(rnew(i,28, 2),i=1,4)/1.676e-09,0.6861,3.530e+03,4.538e+08/
      data(rnew(i,29, 2),i=1,4)/1.686e-09,0.6824,4.031e+03,4.948e+08/
      data(rnew(i,30, 2),i=1,4)/1.758e-09,0.6834,4.254e+03,5.258e+08/
      data(rnew(i, 3, 3),i=1,4)/1.036e-11,0.3880,1.077e+02,1.177e+07/
      data(rnew(i, 4, 3),i=1,4)/2.338e-11,0.4211,3.647e+02,1.215e+07/
      data(rnew(i, 5, 3),i=1,4)/4.487e-11,0.4644,5.371e+02,1.465e+07/
      data(rnew(i, 6, 3),i=1,4)/8.540e-11,0.5247,5.014e+02,1.479e+07/
      data(rnew(i, 7, 3),i=1,4)/1.169e-10,0.5470,6.793e+02,1.650e+07/
      data(rnew(i, 8, 3),i=1,4)/2.053e-10,0.6019,4.772e+02,1.711e+07/
      data(rnew(i, 9, 3),i=1,4)/2.739e-10,0.6188,5.033e+02,2.064e+07/
      data(rnew(i,10, 3),i=1,4)/3.200e-10,0.6198,6.329e+02,2.616e+07/
      data(rnew(i,11, 3),i=1,4)/3.873e-10,0.6295,7.000e+02,2.989e+07/
      data(rnew(i,12, 3),i=1,4)/4.284e-10,0.6287,8.748e+02,3.586e+07/
      data(rnew(i,13, 3),i=1,4)/4.881e-10,0.6326,9.941e+02,4.085e+07/
      data(rnew(i,14, 3),i=1,4)/5.373e-10,0.6337,1.164e+03,4.677e+07/
      data(rnew(i,15, 3),i=1,4)/5.876e-10,0.6354,1.341e+03,5.292e+07/
      data(rnew(i,16, 3),i=1,4)/6.571e-10,0.6400,1.452e+03,5.796e+07/
      data(rnew(i,17, 3),i=1,4)/7.076e-10,0.6397,1.653e+03,6.555e+07/
      data(rnew(i,18, 3),i=1,4)/7.538e-10,0.6388,1.889e+03,7.306e+07/
      data(rnew(i,19, 3),i=1,4)/8.182e-10,0.6411,2.044e+03,8.057e+07/
      data(rnew(i,20, 3),i=1,4)/8.577e-10,0.6403,2.334e+03,8.850e+07/
      data(rnew(i,21, 3),i=1,4)/9.162e-10,0.6413,2.543e+03,9.690e+07/
      data(rnew(i,22, 3),i=1,4)/9.844e-10,0.6440,2.708e+03,1.044e+08/
      data(rnew(i,23, 3),i=1,4)/1.020e-09,0.6427,3.057e+03,1.140e+08/
      data(rnew(i,24, 3),i=1,4)/1.091e-09,0.6445,3.225e+03,1.229e+08/
      data(rnew(i,25, 3),i=1,4)/1.151e-09,0.6451,3.461e+03,1.334e+08/
      data(rnew(i,26, 3),i=1,4)/1.198e-09,0.6443,3.789e+03,1.437e+08/
      data(rnew(i,27, 3),i=1,4)/1.211e-09,0.6406,4.357e+03,1.572e+08/
      data(rnew(i,28, 3),i=1,4)/1.288e-09,0.6440,4.506e+03,1.651e+08/
      data(rnew(i,29, 3),i=1,4)/1.372e-09,0.6472,4.627e+03,1.740e+08/
      data(rnew(i,30, 3),i=1,4)/1.412e-09,0.6454,5.053e+03,1.891e+08/
      data(rnew(i,11,11),i=1,4)/5.641e-12,0.1749,3.077e+02,2.617e+06/
      data(rnew(i,12,11),i=1,4)/1.920e-11,0.3028,4.849e+02,5.890e+06/
      data(rnew(i,13,11),i=1,4)/3.753e-11,0.3585,6.848e+02,9.035e+06/
      data(rnew(i,14,11),i=1,4)/5.942e-11,0.3930,8.962e+02,1.213e+07/
      data(rnew(i,15,11),i=1,4)/1.721e-10,0.5429,2.848e+02,3.975e+07/
      data(rnew(i,16,11),i=1,4)/3.502e-10,0.6266,1.532e+02,1.755e+07/
      data(rnew(i,17,11),i=1,4)/2.502e-10,0.5580,5.303e+02,4.558e+07/
      data(rnew(i,18,11),i=1,4)/2.862e-10,0.5621,7.002e+02,4.885e+07/
      data(rnew(i,19,11),i=1,4)/2.757e-10,0.5364,1.204e+03,7.013e+07/
      data(rnew(i,20,11),i=1,4)/5.273e-10,0.6281,5.329e+02,3.188e+07/
      data(rnew(i,21,11),i=1,4)/3.890e-10,0.5645,1.391e+03,6.295e+07/
      data(rnew(i,22,11),i=1,4)/4.207e-10,0.5646,1.688e+03,6.872e+07/
      data(rnew(i,23,11),i=1,4)/4.605e-10,0.5659,1.949e+03,7.419e+07/
      data(rnew(i,24,11),i=1,4)/4.975e-10,0.5655,2.257e+03,8.072e+07/
      data(rnew(i,25,11),i=1,4)/5.349e-10,0.5658,2.577e+03,8.710e+07/
      data(rnew(i,26,11),i=1,4)/7.688e-10,0.6173,1.653e+03,6.161e+07/
      data(rnew(i,27,11),i=1,4)/5.850e-10,0.5598,3.538e+03,1.052e+08/
      data(rnew(i,28,11),i=1,4)/6.347e-10,0.5631,3.780e+03,1.116e+08/
      data(rnew(i,29,11),i=1,4)/6.619e-10,0.5602,4.322e+03,1.210e+08/
      data(rnew(i,30,11),i=1,4)/7.002e-10,0.5612,4.726e+03,1.287e+08/
      data(rnew(i, 6, 4),i=1,4)/2.020E-09,0.7798,6.690E-01,2.425E+06/
      data(rnew(i, 7, 4),i=1,4)/1.595E-11,0.3529,9.870E+03,2.584E+07/
      data(rnew(i, 8, 4),i=1,4)/2.008E-11,0.3567,1.520E+04,3.843E+07/
      data(rnew(i,10, 4),i=1,4)/2.793E-11,0.3533,3.017E+04,7.872E+07/
      data(rnew(i, 6, 5),i=1,4)/8.577E-10,0.7837,7.286E-01,1.140E+07/
      data(rnew(i, 7, 5),i=1,4)/7.039E-10,0.8607,2.203E+00,3.029E+06/
      data(rnew(i, 8, 5),i=1,4)/1.542E-10,0.6712,1.775E+02,1.535E+08/
      data(rnew(i,10, 5),i=1,4)/2.515E-10,0.7011,3.028E+02,8.903E+06/
      data(rnew(i, 6, 6),i=1,4)/7.651E-09,0.8027,1.193E-03,9.334E+12/
      data(rnew(i, 7, 6),i=1,4)/5.989E-10,0.7560,1.377E+00,1.517E+09/
      data(rnew(i, 8, 6),i=1,4)/3.672E-09,0.7676,2.900E-01,8.521E+07/
      data(rnew(i,10, 6),i=1,4)/9.353E-11,0.6270,9.031E+02,2.387E+07/
      data(rnew(i, 7, 7),i=1,4)/1.243E-08,0.7022,1.136E-03,1.015E+13/
      data(rnew(i, 8, 7),i=1,4)/1.816E-08,0.7170,6.717E-03,3.286E+12/
      data(rnew(i,10, 7),i=1,4)/4.227E-11,0.5395,1.687E+03,1.491E+17/
      data(rnew(i, 8, 8),i=1,4)/1.341E-10,0.6159,1.673E+00,6.366E+16/
      data(rnew(i,10, 8),i=1,4)/9.563E-12,0.3067,9.768E+03,4.851E+17/
      data(rnew(i,10, 9),i=1,4)/5.417E-08,0.6930,1.179E-03,1.060E+07/
      data(rnew(i,10,10),i=1,4)/5.023E-12,0.2420,3.181E+02,1.450E+18/
      data(rnew(i,26,12),i=1,4)/8.110E-10,0.5442,1.755E+03,6.799E+07/ 
      data(rnew(i,26,13),i=1,4)/6.967E-10,0.5602,1.727E+03,5.618E+07/ 
      data(rnew(i,26,14),i=1,4)/3.236E-08,0.3247,2.338E+01,8.337E+12/
      data(rnew(i,26,15),i=1,4)/2.664E-08,0.3285,2.297E+01,6.672E+12/
      data(rnew(i,26,16),i=1,4)/2.165E-08,0.3403,2.195E+01,6.383E+12/
      data(rnew(i,26,17),i=1,4)/2.460E-08,0.3387,1.487E+01,5.228E+12/
      data(rnew(i,26,18),i=1,4)/1.907E-08,0.3768,1.216E+01,5.431E+12/
      data(rnew(i,26,19),i=1,4)/1.439E-08,0.4170,1.006E+01,4.898E+12/
      data(rnew(i,26,20),i=1,4)/1.184E-08,0.4798,5.883E+00,2.582E+12/
      data(rnew(i,26,21),i=1,4)/1.036E-08,0.5428,2.743E+00,1.014E+12/
      data(rnew(i,26,22),i=1,4)/8.288E-09,0.6012,1.216E+00,1.182E+12/
      data(rnew(i,26,23),i=1,4)/6.330E-09,0.6355,5.457E-01,8.545E+11/
      data(rnew(i,26,24),i=1,4)/5.422E-09,0.5067,4.998E-01,8.079E+11/
      data(rnew(i,26,25),i=1,4)/6.076E-09,0.3112,3.401E-01,1.960E+12/
      data(rnew(i,26,26),i=1,4)/8.945E-09,0.2156,4.184E-02,5.353E+13/
      data(fe(1,i),i=4,13)/4.33e-10,3.91e-10,3.49e-10,3.16e-10,2.96e-10,2.59e-10,2.24e-10,1.91e-10,1.68e-10,1.46e-10/
      data(fe(2,i),i=4,13)/0.531,0.523,0.521,0.534,0.557,0.567,0.579,0.601,0.602,0.597/
      data(fe(3,i),i=4,13)/5.77e-02,6.15e-02,6.22e-02,6.02e-02,5.79e-02,5.65e-02,5.49e-02,5.10e-02,5.07e-02,5.22e-02/
      end
!==========================================================================
      subroutine phfit2_b(nz,ne_int,is,e_photo,s_photo,thresh)
!*** the parameters list modified by MG; Sun Nov 10 17:50:58 MSK 1996
!*** thresh - ionization energy in eV
!***
!*** Version 2. March 25, 1996.
!*** Written by D. A. Verner, verner@pa.uky.edu
!*** Inner-shell ionization energies of some low-ionized species are slightly
!*** improved to fit smoothly the experimental inner-shell ionization energies 
!*** of neutral atoms.
!******************************************************************************
!*** This subroutine calculates partial photoionization cross sections
!*** for all ionization stages of all atoms from H to Zn (Z=30) by use of
!*** the following fit parameters:
!*** Outer shells of the Opacity Project (OP) elements:
!***    Verner, Ferland, Korista, Yakovlev, 1996, ApJ, in press.
!*** Inner shells of all elements, and outer shells of the non-OP elements:
!***    Verner and Yakovlev, 1995, A&AS, 109, 125
!*** Input parameters:  nz - atomic number from 1 to 30 (integer) 
!***                    ne_int - number of electrons from 1 to iz (integer)
!***                    is - shell number (integer)
!***                    e - photon energy, eV 
!*** Output parameter:  s - photoionization cross section, Mb
!*** Shell numbers:
!*** 1 - 1s, 2 - 2s, 3 - 2p, 4 - 3s, 5 - 3p, 6 - 3d, 7 - 4s. 
!*** If a species in the ground state has no electrons on the given shell,
!*** the subroutine returns s=0.
!******************************************************************************
!! ADDED TO THE ORIGINAL CODE !!

      implicit none 
      integer :: l,ninn,ntot,nint,nout
      real :: ph1,ph2,a,b,einn,p1,q,x,y,z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: nz,ne_int,is
      real :: e_photo,s_photo,thresh
      common/l/l(7)
      common/ninn/ninn(30)
      common/ntot/ntot(30)
      common/ph1/ph1(6,30,30,7)
      common/ph2/ph2(7,30,30)

      s_photo=0.0
!      thresh=1e-20
      thresh=0.0

      if(nz.lt.1.or.nz.gt.30)return
      if(ne_int.lt.1.or.ne_int.gt.nz)return
      nout=ntot(ne_int)
      if(nz.eq.ne_int.and.nz.gt.18)nout=7
      if(nz.eq.(ne_int+1).and.(nz.eq.20.or.nz.eq.21.or.nz.eq.22.or.nz&
      .eq.25.or.nz.eq.26))nout=7
      if(is.gt.nout)return
      thresh=ph1(1,nz,ne_int,is)
      if(e_photo.lt.ph1(1,nz,ne_int,is))return
!c>      print*,'passed',thresh
      nint=ninn(ne_int)
      if(nz.eq.15.or.nz.eq.17.or.nz.eq.19.or.&
      (nz.gt.20.and.nz.ne.26))then
         einn=0.0
      else
         if(ne_int.lt.3)then
            einn=1.0e+30
         else
            einn=ph1(1,nz,ne_int,nint)
         endif
      endif
!c>      print*,' is,nout,nint,e_photo,einn:',is,nout,nint,e,einn
      if(is.lt.nout.and.is.gt.nint.and.e_photo.lt.einn)return
!c>      print*,'passed 2'
      if(is.le.nint.or.e_photo.ge.einn)then
         p1=-ph1(5,nz,ne_int,is)
         y=e_photo/ph1(2,nz,ne_int,is)
         q=-0.5*p1-l(is)-5.5
         a=ph1(3,nz,ne_int,is)*((y-1.0)**2+ph1(6,nz,ne_int,is)**2)
         b=sqrt(y/ph1(4,nz,ne_int,is))+1.0
         s_photo=a*y**q*b**p1
!         print *,ph1(2,nz,ne_int,is),nz,ne_int,is
      else
         p1=-ph2(4,nz,ne_int)
         q=-0.5*p1-5.5
         x=e_photo/ph2(1,nz,ne_int)-ph2(6,nz,ne_int)
         z=sqrt(x*x+ph2(7,nz,ne_int)**2)
         a=ph2(2,nz,ne_int)*((x-1.0)**2+ph2(5,nz,ne_int)**2)
         b=1.0+sqrt(z/ph2(3,nz,ne_int))
         s_photo=a*z**q*b**p1
!         print *,ph2(2,nz,ne_int),nz,ne_int
      endif
      return
      end subroutine phfit2_b

!=======================================================
      BLOCK DATA BDATA
      COMMON/L/L(7)
      COMMON/NINN/NINN(30)
      COMMON/NTOT/NTOT(30)
      COMMON/PH1/PH1(6,30,30,7)
      COMMON/PH2/PH2(7,30,30)
      DATA (L(I),I=1,7) /0,0,1,0,1,2,0/
      DATA (NINN(I),I=1,30) /0,0,1,1,1,1,1,1,1,1,3,3, &
      3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5/
      DATA (NTOT(I),I=1,30) /1,1,2,2,3,3,3,3,3,3,4,4, &
      5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7/
      DATA (PH1(I, 1, 1, 1),I=1,6) /1.360E+01, 4.298E-01, 5.475E+04, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 2, 1, 1),I=1,6) /5.442E+01, 1.720E+00, 1.369E+04, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 2, 2, 1),I=1,6) /2.459E+01, 5.996E+00, 4.470E+03, 2.199E+00, 6.098E+00, 0.000E+00/
      DATA (PH1(I, 3, 1, 1),I=1,6) /1.225E+02, 3.871E+00, 6.083E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 3, 2, 1),I=1,6) /7.564E+01, 2.006E+01, 3.201E+02, 7.391E+00, 2.916E+00, 0.000E+00/
      DATA (PH1(I, 3, 3, 1),I=1,6) /6.439E+01, 2.740E+01, 1.564E+02, 3.382E+01, 1.490E+00, 0.000E+00/
      DATA (PH1(I, 3, 3, 2),I=1,6) /5.392E+00, 3.466E+00, 4.774E+01, 2.035E+01, 4.423E+00, 0.000E+00/
      DATA (PH1(I, 4, 1, 1),I=1,6) /2.177E+02, 6.879E+00, 3.422E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 4, 2, 1),I=1,6) /1.539E+02, 1.760E+01, 5.458E+02, 1.719E+01, 3.157E+00, 0.000E+00/
      DATA (PH1(I, 4, 3, 1),I=1,6) /1.299E+02, 4.759E+01, 9.796E+01, 1.166E+03, 1.022E+00, 0.000E+00/
      DATA (PH1(I, 4, 3, 2),I=1,6) /1.821E+01, 6.685E+00, 4.904E+01, 3.419E+01, 3.738E+00, 0.000E+00/
      DATA (PH1(I, 4, 4, 1),I=1,6) /1.193E+02, 4.093E+01, 1.306E+02, 1.212E+02, 1.348E+00, 0.000E+00/
      DATA (PH1(I, 4, 4, 2),I=1,6) /9.323E+00, 3.427E+00, 1.423E+03, 2.708E+00, 1.064E+01, 0.000E+00/
      DATA (PH1(I, 5, 1, 1),I=1,6) /3.402E+02, 1.075E+01, 2.190E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 5, 2, 1),I=1,6) /2.594E+02, 3.336E+01, 2.846E+02, 2.163E+01, 2.624E+00, 0.000E+00/
      DATA (PH1(I, 5, 3, 1),I=1,6) /2.274E+02, 6.592E+01, 8.605E+01, 1.906E+02, 1.210E+00, 0.000E+00/
      DATA (PH1(I, 5, 3, 2),I=1,6) /3.793E+01, 1.054E+01, 4.273E+01, 4.433E+01, 3.501E+00, 0.000E+00/
      DATA (PH1(I, 5, 4, 1),I=1,6) /2.098E+02, 5.984E+01, 1.037E+02, 7.915E+01, 1.436E+00, 0.000E+00/
      DATA (PH1(I, 5, 4, 2),I=1,6) /2.516E+01, 2.403E+00, 1.530E+02, 9.203E+00, 9.374E+00, 0.000E+00/
      DATA (PH1(I, 5, 5, 1),I=1,6) /1.940E+02, 6.155E+01, 9.698E+01, 7.354E+01, 1.438E+00, 0.000E+00/
      DATA (PH1(I, 5, 5, 2),I=1,6) /1.405E+01, 6.495E+00, 1.525E+04, 1.352E+00, 1.218E+01, 0.000E+00/
      DATA (PH1(I, 5, 5, 3),I=1,6) /8.298E+00, 6.658E+00, 3.643E+02, 9.500E+00, 5.469E+00, 5.608E-01/
      DATA (PH1(I, 6, 1, 1),I=1,6) /4.900E+02, 1.548E+01, 1.521E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 6, 2, 1),I=1,6) /3.921E+02, 4.624E+01, 2.344E+02, 2.183E+01, 2.581E+00, 0.000E+00/
      DATA (PH1(I, 6, 3, 1),I=1,6) /3.522E+02, 8.412E+01, 8.111E+01, 7.459E+01, 1.428E+00, 0.000E+00/
      DATA (PH1(I, 6, 3, 2),I=1,6) /6.449E+01, 7.843E+00, 7.109E+01, 2.828E+01, 4.754E+00, 0.000E+00/
      DATA (PH1(I, 6, 4, 1),I=1,6) /3.289E+02, 8.370E+01, 8.067E+01, 7.471E+01, 1.442E+00, 0.000E+00/
      DATA (PH1(I, 6, 4, 2),I=1,6) /4.789E+01, 3.942E+00, 1.219E+02, 1.499E+01, 7.489E+00, 0.000E+00/
      DATA (PH1(I, 6, 5, 1),I=1,6) /3.076E+02, 9.113E+01, 6.649E+01, 9.609E+01, 1.338E+00, 0.000E+00/
      DATA (PH1(I, 6, 5, 2),I=1,6) /3.047E+01, 2.991E+00, 1.184E+03, 3.085E+00, 1.480E+01, 0.000E+00/
      DATA (PH1(I, 6, 5, 3),I=1,6) /2.438E+01, 1.094E+01, 1.792E+02, 3.308E+01, 4.150E+00, 5.276E-01/
      DATA (PH1(I, 6, 6, 1),I=1,6) /2.910E+02, 8.655E+01, 7.421E+01, 5.498E+01, 1.503E+00, 0.000E+00/
      DATA (PH1(I, 6, 6, 2),I=1,6) /1.939E+01, 1.026E+01, 4.564E+03, 1.568E+00, 1.085E+01, 0.000E+00/
      DATA (PH1(I, 6, 6, 3),I=1,6) /1.126E+01, 9.435E+00, 1.152E+03, 5.687E+00, 6.336E+00, 4.474E-01/
      DATA (PH1(I, 7, 1, 1),I=1,6) /6.671E+02, 2.108E+01, 1.117E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 7, 2, 1),I=1,6) /5.521E+02, 6.943E+01, 1.519E+02, 2.627E+01, 2.315E+00, 0.000E+00/
      DATA (PH1(I, 7, 3, 1),I=1,6) /5.043E+02, 1.060E+02, 7.304E+01, 5.547E+01, 1.538E+00, 0.000E+00/
      DATA (PH1(I, 7, 3, 2),I=1,6) /9.789E+01, 1.862E+01, 3.447E+01, 4.231E+01, 3.606E+00, 0.000E+00/
      DATA (PH1(I, 7, 4, 1),I=1,6) /4.753E+02, 1.070E+02, 7.046E+01, 5.342E+01, 1.552E+00, 0.000E+00/
      DATA (PH1(I, 7, 4, 2),I=1,6) /7.747E+01, 6.225E+00, 1.110E+02, 1.733E+01, 6.719E+00, 0.000E+00/
      DATA (PH1(I, 7, 5, 1),I=1,6) /4.473E+02, 1.220E+02, 5.235E+01, 9.428E+01, 1.335E+00, 0.000E+00/
      DATA (PH1(I, 7, 5, 2),I=1,6) /5.545E+01, 5.853E+00, 1.908E+02, 6.264E+00, 9.711E+00, 0.000E+00/
      DATA (PH1(I, 7, 5, 3),I=1,6) /4.745E+01, 1.925E+01, 9.400E+01, 1.152E+02, 3.194E+00, 5.496E-01/
      DATA (PH1(I, 7, 6, 1),I=1,6) /4.236E+02, 1.242E+02, 5.002E+01, 9.100E+01, 1.335E+00, 0.000E+00/
      DATA (PH1(I, 7, 6, 2),I=1,6) /3.796E+01, 1.094E+01, 7.483E+02, 2.793E+00, 9.956E+00, 0.000E+00/
      DATA (PH1(I, 7, 6, 3),I=1,6) /2.960E+01, 1.827E+01, 1.724E+02, 8.893E+01, 3.348E+00, 4.209E-01/
      DATA (PH1(I, 7, 7, 1),I=1,6) /4.048E+02, 1.270E+02, 4.748E+01, 1.380E+02, 1.252E+00, 0.000E+00/
      DATA (PH1(I, 7, 7, 2),I=1,6) /2.541E+01, 1.482E+01, 7.722E+02, 2.306E+00, 9.139E+00, 0.000E+00/
      DATA (PH1(I, 7, 7, 3),I=1,6) /1.453E+01, 1.164E+01, 1.029E+04, 2.361E+00, 8.821E+00, 4.239E-01/
      DATA (PH1(I, 8, 1, 1),I=1,6) /8.714E+02, 2.754E+01, 8.554E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 8, 2, 1),I=1,6) /7.393E+02, 8.709E+01, 1.329E+02, 2.535E+01, 2.336E+00, 0.000E+00/
      DATA (PH1(I, 8, 3, 1),I=1,6) /6.837E+02, 1.354E+02, 6.029E+01, 5.682E+01, 1.533E+00, 0.000E+00/
      DATA (PH1(I, 8, 3, 2),I=1,6) /1.381E+02, 9.141E+00, 6.896E+01, 3.896E+01, 4.943E+00, 0.000E+00/
      DATA (PH1(I, 8, 4, 1),I=1,6) /6.491E+02, 1.377E+02, 5.735E+01, 5.486E+01, 1.540E+00, 0.000E+00/
      DATA (PH1(I, 8, 4, 2),I=1,6) /1.139E+02, 8.467E+00, 9.641E+01, 2.287E+01, 6.061E+00, 0.000E+00/
      DATA (PH1(I, 8, 5, 1),I=1,6) /6.144E+02, 1.593E+02, 4.123E+01, 1.141E+02, 1.287E+00, 0.000E+00/
      DATA (PH1(I, 8, 5, 2),I=1,6) /8.737E+01, 7.942E+00, 1.063E+02, 9.708E+00, 8.183E+00, 0.000E+00/
      DATA (PH1(I, 8, 5, 3),I=1,6) /7.741E+01, 1.937E+01, 1.619E+02, 7.312E+01, 3.648E+00, 3.760E-02/
      DATA (PH1(I, 8, 6, 1),I=1,6) /5.840E+02, 1.620E+02, 3.939E+01, 1.104E+02, 1.289E+00, 0.000E+00/
      DATA (PH1(I, 8, 6, 2),I=1,6) /6.551E+01, 1.594E+01, 1.217E+02, 6.156E+00, 7.271E+00, 0.000E+00/
      DATA (PH1(I, 8, 6, 3),I=1,6) /5.494E+01, 2.067E+01, 2.318E+02, 7.136E+01, 3.618E+00, 5.538E-02/
      DATA (PH1(I, 8, 7, 1),I=1,6) /5.581E+02, 1.690E+02, 3.584E+01, 1.894E+02, 1.185E+00, 0.000E+00/
      DATA (PH1(I, 8, 7, 2),I=1,6) /4.599E+01, 1.759E+01, 1.962E+02, 4.020E+00, 7.999E+00, 0.000E+00/
      DATA (PH1(I, 8, 7, 3),I=1,6) /3.512E+01, 1.745E+01, 5.186E+02, 1.728E+01, 4.995E+00, 2.182E-02/
      DATA (PH1(I, 8, 8, 1),I=1,6) /5.380E+02, 1.774E+02, 3.237E+01, 3.812E+02, 1.083E+00, 0.000E+00/
      DATA (PH1(I, 8, 8, 2),I=1,6) /2.848E+01, 1.994E+01, 2.415E+02, 3.241E+00, 8.037E+00, 0.000E+00/
      DATA (PH1(I, 8, 8, 3),I=1,6) /1.362E+01, 1.391E+01, 1.220E+05, 1.364E+00, 1.140E+01, 4.103E-01/
      DATA (PH1(I, 9, 1, 1),I=1,6) /1.103E+03, 3.485E+01, 6.759E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 9, 2, 1),I=1,6) /9.539E+02, 1.131E+02, 1.039E+02, 2.657E+01, 2.255E+00, 0.000E+00/
      DATA (PH1(I, 9, 3, 1),I=1,6) /8.905E+02, 1.711E+02, 4.890E+01, 6.137E+01, 1.501E+00, 0.000E+00/
      DATA (PH1(I, 9, 3, 2),I=1,6) /1.852E+02, 6.423E+00, 8.497E+01, 6.075E+01, 5.214E+00, 0.000E+00/
      DATA (PH1(I, 9, 4, 1),I=1,6) /8.502E+02, 1.737E+02, 4.668E+01, 5.876E+01, 1.511E+00, 0.000E+00/
      DATA (PH1(I, 9, 4, 2),I=1,6) /1.572E+02, 1.146E+01, 8.453E+01, 2.457E+01, 5.771E+00, 0.000E+00/
      DATA (PH1(I, 9, 5, 1),I=1,6) /8.091E+02, 1.925E+02, 3.680E+01, 7.933E+01, 1.377E+00, 0.000E+00/
      DATA (PH1(I, 9, 5, 2),I=1,6) /1.262E+02, 1.889E+01, 7.149E+01, 1.194E+01, 6.030E+00, 0.000E+00/
      DATA (PH1(I, 9, 5, 3),I=1,6) /1.142E+02, 2.528E+01, 1.419E+02, 7.089E+01, 3.628E+00, 3.418E-01/
      DATA (PH1(I, 9, 6, 1),I=1,6) /7.709E+02, 1.603E+02, 5.302E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I, 9, 6, 2),I=1,6) /9.957E+01, 2.629E+01,3.563E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I, 9, 6, 3),I=1,6) /8.714E+01, 2.292E+01,3.005E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I, 9, 7, 1),I=1,6) /7.392E+02, 2.055E+02,3.144E+01, 1.230E+02, 1.263E+00, 0.000E+00/
      DATA (PH1(I, 9, 7, 2),I=1,6) /7.610E+01, 2.307E+01,7.140E+01, 7.282E+00, 6.543E+00, 0.000E+00/
      DATA (PH1(I, 9, 7, 3),I=1,6) /6.271E+01, 2.744E+01,2.458E+02, 7.329E+01, 3.596E+00, 1.390E-01/
      DATA (PH1(I, 9, 8, 1),I=1,6) /7.122E+02, 1.660E+02,4.798E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I, 9, 8, 2),I=1,6) /5.459E+01, 2.882E+01,2.439E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I, 9, 8, 3),I=1,6) /3.497E+01, 2.658E+01,2.649E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I, 9, 9, 1),I=1,6) /6.940E+02, 2.390E+02,2.295E+01, 1.257E+03, 9.638E-01, 0.000E+00/
      DATA (PH1(I, 9, 9, 2),I=1,6) /3.786E+01, 2.568E+01,1.097E+02, 4.297E+00, 7.303E+00, 0.000E+00/
      DATA (PH1(I, 9, 9, 3),I=1,6) /1.742E+01, 1.658E+01,2.775E+05, 1.242E+00, 1.249E+01, 3.857E-01/
      DATA (PH1(I,10, 1, 1),I=1,6) /1.362E+03, 4.304E+01,5.475E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,10, 2, 1),I=1,6) /1.196E+03, 1.586E+02,6.695E+01, 3.352E+01, 2.002E+00, 0.000E+00/
      DATA (PH1(I,10, 3, 1),I=1,6) /1.125E+03, 2.193E+02,3.719E+01, 8.181E+01, 1.396E+00, 0.000E+00/
      DATA (PH1(I,10, 3, 2),I=1,6) /2.391E+02, 2.859E+01,2.897E+01, 3.836E+01, 3.992E+00, 0.000E+00/
      DATA (PH1(I,10, 4, 1),I=1,6) /1.078E+03, 1.886E+02,5.011E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,10, 4, 2),I=1,6) /2.073E+02, 2.822E+01,4.978E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,10, 5, 1),I=1,6) /1.031E+03, 2.308E+02,3.232E+01, 7.167E+01, 1.411E+00, 0.000E+00/
      DATA (PH1(I,10, 5, 2),I=1,6) /1.719E+02, 2.496E+01,5.575E+01, 1.427E+01, 5.562E+00, 0.000E+00/
      DATA (PH1(I,10, 5, 3),I=1,6) /1.579E+02, 2.740E+01,1.717E+02, 6.423E+01, 3.819E+00, 1.680E-01/
      DATA (PH1(I,10, 6, 1),I=1,6) /9.873E+02, 2.343E+02,3.097E+01, 6.907E+01, 1.416E+00, 0.000E+00/
      DATA (PH1(I,10, 6, 2),I=1,6) /1.415E+02, 2.615E+01,5.421E+01, 1.117E+01, 5.907E+00, 0.000E+00/
      DATA (PH1(I,10, 6, 3),I=1,6) /1.262E+02, 3.316E+01,1.962E+02, 7.397E+01, 3.563E+00, 2.895E-01/
      DATA (PH1(I,10, 7, 1),I=1,6) /9.480E+02, 2.360E+02,3.027E+01, 6.734E+01, 1.422E+00, 0.000E+00/
      DATA (PH1(I,10, 7, 2),I=1,6) /1.132E+02, 3.092E+01,4.190E+01, 1.139E+01, 5.554E+00, 0.000E+00/
      DATA (PH1(I,10, 7, 3),I=1,6) /9.712E+01, 3.494E+01,2.214E+02, 7.203E+01, 3.563E+00, 2.573E-01/
      DATA (PH1(I,10, 8, 1),I=1,6) /9.131E+02, 2.387E+02,2.943E+01, 6.525E+01, 1.424E+00, 0.000E+00/
      DATA (PH1(I,10, 8, 2),I=1,6) /8.721E+01, 3.192E+01,3.797E+01, 1.054E+01, 5.652E+00, 0.000E+00/
      DATA (PH1(I,10, 8, 3),I=1,6) /6.346E+01, 3.485E+01,2.451E+02, 6.937E+01, 3.645E+00, 2.042E-01/
      DATA (PH1(I,10, 9, 1),I=1,6) /8.831E+02, 2.452E+02,2.783E+01, 6.075E+01, 1.420E+00, 0.000E+00/
      DATA (PH1(I,10, 9, 2),I=1,6) /6.374E+01, 3.187E+01,4.025E+01, 8.495E+00, 6.038E+00, 0.000E+00/
      DATA (PH1(I,10, 9, 3),I=1,6) /4.096E+01, 3.428E+01,2.766E+02, 4.179E+01, 4.029E+00, 3.052E-01/
      DATA (PH1(I,10,10, 1),I=1,6) /8.701E+02, 3.144E+02,1.664E+01, 2.042E+05, 8.450E-01, 0.000E+00/
      DATA (PH1(I,10,10, 2),I=1,6) /4.847E+01, 3.204E+01,5.615E+01, 5.808E+00, 6.678E+00, 0.000E+00/
      DATA (PH1(I,10,10, 3),I=1,6) /2.156E+01, 2.000E+01,1.691E+04, 2.442E+00, 1.043E+01, 3.345E-01/
      DATA (PH1(I,11, 1, 1),I=1,6) /1.649E+03, 5.211E+01,4.525E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,11, 2, 1),I=1,6) /1.465E+03, 2.268E+02,3.995E+01, 5.315E+01, 1.678E+00, 0.000E+00/
      DATA (PH1(I,11, 3, 1),I=1,6) /1.386E+03, 2.406E+02,3.850E+01, 5.198E+01, 1.575E+00, 0.000E+00/
      DATA (PH1(I,11, 3, 2),I=1,6) /2.999E+02, 1.758E+01,4.531E+01, 4.549E+01, 4.689E+00, 0.000E+00/
      DATA (PH1(I,11, 4, 1),I=1,6) /1.335E+03, 2.464E+02,3.613E+01, 4.968E+01, 1.579E+00, 0.000E+00/
      DATA (PH1(I,11, 4, 2),I=1,6) /2.642E+02, 2.097E+01,6.127E+01, 2.644E+01, 5.246E+00, 0.000E+00/
      DATA (PH1(I,11, 5, 1),I=1,6) /1.281E+03, 2.777E+02,2.745E+01, 7.512E+01, 1.397E+00, 0.000E+00/
      DATA (PH1(I,11, 5, 2),I=1,6) /2.244E+02, 2.947E+01,4.802E+01, 1.621E+01, 5.395E+00, 0.000E+00/
      DATA (PH1(I,11, 5, 3),I=1,6) /2.085E+02, 3.552E+01,1.374E+02, 7.062E+01, 3.675E+00, 1.613E-01/
      DATA (PH1(I,11, 6, 1),I=1,6) /1.230E+03, 2.806E+02,2.654E+01, 7.237E+01, 1.405E+00, 0.000E+00/
      DATA (PH1(I,11, 6, 2),I=1,6) /1.899E+02, 3.303E+01,4.218E+01, 1.377E+01, 5.458E+00, 0.000E+00/
      DATA (PH1(I,11, 6, 3),I=1,6) /1.722E+02, 4.212E+01,1.650E+02, 9.649E+01, 3.367E+00, 2.182E-01/
      DATA (PH1(I,11, 7, 1),I=1,6) /1.185E+03, 2.360E+02,3.753E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,11, 7, 2),I=1,6) /1.567E+02, 3.868E+01,2.607E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,11, 7, 3),I=1,6) /1.384E+02, 3.523E+01,3.028E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,11, 8, 1),I=1,6) /1.143E+03, 2.841E+02,2.551E+01, 6.977E+01, 1.414E+00, 0.000E+00/
      DATA (PH1(I,11, 8, 2),I=1,6) /1.269E+02, 3.807E+01,3.144E+01, 1.250E+01, 5.405E+00, 0.000E+00/
      DATA (PH1(I,11, 8, 3),I=1,6) /9.892E+01, 4.314E+01,2.232E+02, 8.352E+01, 3.511E+00, 1.644E-01/
      DATA (PH1(I,11, 9, 1),I=1,6) /1.118E+03, 2.426E+02,3.466E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,11, 9, 2),I=1,6) /9.945E+01, 4.130E+01,1.933E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,11, 9, 3),I=1,6) /7.162E+01, 3.957E+01,2.748E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,11,10, 1),I=1,6) /1.097E+03, 3.266E+02,1.889E+01, 2.527E+02, 1.120E+00, 0.000E+00/
      DATA (PH1(I,11,10, 2),I=1,6) /7.347E+01, 3.939E+01,2.629E+01, 1.103E+01, 5.618E+00, 0.000E+00/
      DATA (PH1(I,11,10, 3),I=1,6) /4.729E+01, 3.911E+01,3.551E+02, 2.257E+01, 4.643E+00, 2.900E-01/
      DATA (PH1(I,11,11, 1),I=1,6) /1.079E+03, 4.216E+02,1.119E+01, 5.642E+07, 7.736E-01, 0.000E+00/
      DATA (PH1(I,11,11, 2),I=1,6) /7.084E+01, 4.537E+01,1.142E+01, 2.395E+02, 3.380E+00, 0.000E+00/
      DATA (PH1(I,11,11, 3),I=1,6) /3.814E+01, 3.655E+01,2.486E+02, 3.222E+02, 3.570E+00, 1.465E-01/
      DATA (PH1(I,11,11, 4),I=1,6) /5.139E+00, 5.968E+00,1.460E+00, 2.557E+07, 3.789E+00, 0.000E+00/
      DATA (PH1(I,12, 1, 1),I=1,6) /1.963E+03, 6.203E+01,3.802E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,12, 2, 1),I=1,6) /1.762E+03, 2.042E+02,6.140E+01, 2.778E+01, 2.161E+00, 0.000E+00/
      DATA (PH1(I,12, 3, 1),I=1,6) /1.675E+03, 2.858E+02,3.290E+01, 5.384E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,12, 3, 2),I=1,6) /3.675E+02, 4.041E+01,2.228E+01, 3.904E+01, 3.986E+00, 0.000E+00/
      DATA (PH1(I,12, 4, 1),I=1,6) /1.618E+03, 2.671E+02,3.699E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,12, 4, 2),I=1,6) /3.282E+02, 3.981E+01,3.936E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,12, 5, 1),I=1,6) /1.558E+03, 3.082E+02,2.721E+01, 5.108E+01, 1.539E+00, 0.000E+00/
      DATA (PH1(I,12, 5, 2),I=1,6) /2.839E+02, 3.967E+01,3.755E+01, 1.873E+01, 4.950E+00, 0.000E+00/
      DATA (PH1(I,12, 5, 3),I=1,6) /2.660E+02, 2.461E+01,3.546E+02, 6.160E+01, 4.366E+00, 6.452E-03/
      DATA (PH1(I,12, 6, 1),I=1,6) /1.501E+03, 2.742E+02,3.407E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,12, 6, 2),I=1,6) /2.444E+02, 4.341E+01,2.832E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,12, 6, 3),I=1,6) /2.249E+02, 4.016E+01,2.442E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,12, 7, 1),I=1,6) /1.449E+03, 2.777E+02,3.278E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,12, 7, 2),I=1,6) /2.076E+02, 4.500E+01,2.434E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,12, 7, 3),I=1,6) /1.865E+02, 4.155E+01,2.901E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,12, 8, 1),I=1,6) /1.400E+03, 3.259E+02,2.350E+01, 6.173E+01, 1.464E+00, 0.000E+00/
      DATA (PH1(I,12, 8, 2),I=1,6) /1.735E+02, 4.710E+01,2.448E+01, 1.636E+01, 4.951E+00, 0.000E+00/
      DATA (PH1(I,12, 8, 3),I=1,6) /1.413E+02, 5.114E+01,2.134E+02, 9.563E+01, 3.438E+00, 5.419E-03/
      DATA (PH1(I,12, 9, 1),I=1,6) /1.356E+03, 3.345E+02,2.212E+01, 7.147E+01, 1.413E+00, 0.000E+00/
      DATA (PH1(I,12, 9, 2),I=1,6) /1.411E+02, 4.810E+01,2.118E+01, 1.675E+01, 4.941E+00, 0.000E+00/
      DATA (PH1(I,12, 9, 3),I=1,6) /1.093E+02, 5.137E+01,2.223E+02, 9.301E+01, 3.513E+00, 5.276E-03/
      DATA (PH1(I,12,10, 1),I=1,6) /1.336E+03, 3.607E+02,1.877E+01, 1.401E+02, 1.238E+00, 0.000E+00/
      DATA (PH1(I,12,10, 2),I=1,6) /1.111E+02, 4.804E+01,1.914E+01, 1.613E+01, 5.050E+00, 0.000E+00/
      DATA (PH1(I,12,10, 3),I=1,6) /8.014E+01, 5.319E+01,2.104E+02, 7.813E+01, 3.630E+00, 2.851E-01/
      DATA (PH1(I,12,11, 1),I=1,6) /1.320E+03, 3.709E+02,1.767E+01, 1.346E+02, 1.225E+00, 0.000E+00/
      DATA (PH1(I,12,11, 2),I=1,6) /9.881E+01, 5.142E+01,1.290E+01, 5.775E+01, 3.903E+00, 0.000E+00/
      DATA (PH1(I,12,11, 3),I=1,6) /6.569E+01, 4.940E+01,2.049E+02, 4.112E+03, 2.995E+00, 2.223E-05/
      DATA (PH1(I,12,11, 4),I=1,6) /1.504E+01, 8.139E+00,3.278E+00, 4.341E+07, 3.610E+00, 0.000E+00/
      DATA (PH1(I,12,12, 1),I=1,6) /1.311E+03, 2.711E+02,3.561E+01, 2.374E+01, 1.952E+00, 0.000E+00/
      DATA (PH1(I,12,12, 2),I=1,6) /9.400E+01, 4.587E+01,1.671E+01, 2.389E+01, 4.742E+00, 0.000E+00/
      DATA (PH1(I,12,12, 3),I=1,6) /5.490E+01, 4.937E+01,2.023E+02, 1.079E+04, 2.960E+00, 1.463E-02/
      DATA (PH1(I,12,12, 4),I=1,6) /7.646E+00, 9.393E+00,3.034E+00, 2.625E+07, 3.923E+00, 0.000E+00/
      DATA (PH1(I,13, 1, 1),I=1,6) /2.304E+03, 7.281E+01,3.239E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,13, 2, 1),I=1,6) /2.086E+03, 2.738E+02,4.036E+01, 3.567E+01, 1.915E+00, 0.000E+00/
      DATA (PH1(I,13, 3, 1),I=1,6) /1.992E+03, 3.137E+02,3.272E+01, 4.295E+01, 1.676E+00, 0.000E+00/
      DATA (PH1(I,13, 3, 2),I=1,6) /4.421E+02, 4.434E+01,2.098E+01, 3.913E+01, 4.062E+00, 0.000E+00/
      DATA (PH1(I,13, 4, 1),I=1,6) /1.929E+03, 3.117E+02,3.223E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,13, 4, 2),I=1,6) /3.994E+02, 4.639E+01,3.516E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,13, 5, 1),I=1,6) /1.862E+03, 3.305E+02,2.857E+01, 3.660E+01, 1.712E+00, 0.000E+00/
      DATA (PH1(I,13, 5, 2),I=1,6) /3.502E+02, 4.550E+01,3.393E+01, 1.945E+01, 4.926E+00, 0.000E+00/
      DATA (PH1(I,13, 5, 3),I=1,6) /3.301E+02, 4.084E+01,1.706E+02, 6.480E+01, 3.927E+00, 1.386E-01/
      DATA (PH1(I,13, 6, 1),I=1,6) /1.799E+03, 3.193E+02,2.989E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,13, 6, 2),I=1,6) /3.065E+02, 5.016E+01,2.604E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,13, 6, 3),I=1,6) /2.846E+02, 4.719E+01,2.229E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,13, 7, 1),I=1,6) /1.739E+03, 3.396E+02,2.653E+01, 3.478E+01, 1.721E+00, 0.000E+00/
      DATA (PH1(I,13, 7, 2),I=1,6) /2.662E+02, 5.726E+01,2.274E+01, 1.983E+01, 4.594E+00, 0.000E+00/
      DATA (PH1(I,13, 7, 3),I=1,6) /2.414E+02, 5.746E+01,1.941E+02, 8.747E+01, 3.460E+00, 1.206E-01/
      DATA (PH1(I,13, 8, 1),I=1,6) /1.688E+03, 3.268E+02,2.786E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,13, 8, 2),I=1,6) /2.268E+02, 5.340E+01,1.982E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,13, 8, 3),I=1,6) /1.905E+02, 5.032E+01,2.888E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,13, 9, 1),I=1,6) /1.634E+03, 3.432E+02,2.572E+01, 3.376E+01, 1.730E+00, 0.000E+00/
      DATA (PH1(I,13, 9, 2),I=1,6) /1.903E+02, 5.651E+01,1.865E+01, 1.904E+01, 4.763E+00, 0.000E+00/
      DATA (PH1(I,13, 9, 3),I=1,6) /1.538E+02, 6.049E+01,2.132E+02, 8.633E+01, 3.522E+00, 1.041E-01/
      DATA (PH1(I,13,10, 1),I=1,6) /1.604E+03, 3.482E+02,2.494E+01, 3.260E+01, 1.732E+00, 0.000E+00/
      DATA (PH1(I,13,10, 2),I=1,6) /1.558E+02, 5.735E+01,1.619E+01, 1.989E+01, 4.749E+00, 0.000E+00/
      DATA (PH1(I,13,10, 3),I=1,6) /1.200E+02, 6.033E+01,2.186E+02, 8.353E+01, 3.609E+00, 9.152E-02/
      DATA (PH1(I,13,11, 1),I=1,6) /1.583E+03, 3.200E+02,2.994E+01, 2.757E+01, 1.876E+00, 0.000E+00/
      DATA (PH1(I,13,11, 2),I=1,6) /1.407E+02, 5.629E+01,1.398E+01, 3.550E+01, 4.267E+00, 0.000E+00/
      DATA (PH1(I,13,11, 3),I=1,6) /1.026E+02, 6.016E+01,1.990E+02, 5.082E+02, 3.110E+00, 2.016E-02/
      DATA (PH1(I,13,11, 4),I=1,6) /2.845E+01, 1.027E+01,4.915E+00, 1.990E+06, 3.477E+00, 0.000E+00/
      DATA (PH1(I,13,12, 1),I=1,6) /1.571E+03, 3.049E+02,3.295E+01, 2.905E+01, 1.898E+00, 0.000E+00/
      DATA (PH1(I,13,12, 2),I=1,6) /1.281E+02, 5.171E+01,1.642E+01, 2.275E+01, 4.810E+00, 0.000E+00/
      DATA (PH1(I,13,12, 3),I=1,6) /8.997E+01, 6.154E+01,1.842E+02, 2.404E+03, 2.920E+00, 7.839E-03/
      DATA (PH1(I,13,12, 4),I=1,6) /1.883E+01, 1.012E+01,6.324E+00, 2.195E+02, 4.481E+00, 0.000E+00/
      DATA (PH1(I,13,13, 1),I=1,6) /1.567E+03, 3.670E+02,2.206E+01, 4.405E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,13,13, 2),I=1,6) /1.256E+02, 5.594E+01,1.425E+01, 3.094E+01, 4.399E+00, 0.000E+00/
      DATA (PH1(I,13,13, 3),I=1,6) /8.040E+01, 6.445E+01,1.735E+02, 1.131E+04, 2.762E+00, 2.337E-02/
      DATA (PH1(I,13,13, 4),I=1,6) /1.133E+01, 1.204E+01,5.384E+00, 4.341E+02, 4.088E+00, 0.000E+00/
      DATA (PH1(I,13,13, 5),I=1,6) /5.986E+00, 1.860E+01,1.828E+02, 2.797E+00, 1.084E+01, 3.076E-01/
      DATA (PH1(I,14, 1, 1),I=1,6) /2.673E+03, 8.447E+01,2.793E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,14, 2, 1),I=1,6) /2.438E+03, 2.752E+02,4.754E+01, 2.848E+01, 2.135E+00, 0.000E+00/
      DATA (PH1(I,14, 3, 1),I=1,6) /2.336E+03, 3.150E+02,3.880E+01, 3.034E+01, 1.925E+00, 0.000E+00/
      DATA (PH1(I,14, 3, 2),I=1,6) /5.235E+02, 4.847E+01,1.966E+01, 3.851E+01, 4.148E+00, 0.000E+00/
      DATA (PH1(I,14, 4, 1),I=1,6) /2.268E+03, 3.599E+02,2.832E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 4, 2),I=1,6) /4.761E+02, 5.350E+01,3.155E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 5, 1),I=1,6) /2.194E+03, 3.356E+02,3.343E+01, 2.487E+01, 1.982E+00, 0.000E+00/
      DATA (PH1(I,14, 5, 2),I=1,6) /4.234E+02, 6.098E+01,2.624E+01, 2.213E+01, 4.520E+00, 0.000E+00/
      DATA (PH1(I,14, 5, 3),I=1,6) /4.014E+02, 3.950E+01,2.223E+02, 6.521E+01, 4.118E+00, 4.760E-05/
      DATA (PH1(I,14, 6, 1),I=1,6) /2.125E+03, 3.679E+02,2.641E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 6, 2),I=1,6) /3.756E+02, 5.745E+01,2.394E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 6, 3),I=1,6) /3.511E+02, 5.486E+01,2.031E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,14, 7, 1),I=1,6) /2.058E+03, 3.442E+02,3.135E+01, 2.260E+01, 2.020E+00, 0.000E+00/
      DATA (PH1(I,14, 7, 2),I=1,6) /3.310E+02, 6.593E+01,2.049E+01, 2.121E+01, 4.515E+00, 0.000E+00/
      DATA (PH1(I,14, 7, 3),I=1,6) /3.032E+02, 5.911E+01,2.326E+02, 6.806E+01, 3.718E+00, 2.813E-05/
      DATA (PH1(I,14, 8, 1),I=1,6) /2.001E+03, 3.759E+02,2.474E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 8, 2),I=1,6) /2.872E+02, 6.086E+01,1.857E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 8, 3),I=1,6) /2.465E+02, 5.786E+01,2.777E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,14, 9, 1),I=1,6) /1.946E+03, 3.799E+02,2.398E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 9, 2),I=1,6) /2.468E+02, 6.237E+01,1.648E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 9, 3),I=1,6) /2.051E+02, 6.007E+01,2.789E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,14,10, 1),I=1,6) /1.887E+03, 3.460E+02,3.100E+01, 1.979E+01, 2.094E+00, 0.000E+00/
      DATA (PH1(I,14,10, 2),I=1,6) /2.076E+02, 6.653E+01,1.463E+01, 2.200E+01, 4.611E+00, 0.000E+00/
      DATA (PH1(I,14,10, 3),I=1,6) /1.668E+02, 7.021E+01,2.109E+02, 8.199E+01, 3.588E+00, 6.329E-05/
      DATA (PH1(I,14,11, 1),I=1,6) /1.868E+03, 3.131E+02,3.784E+01, 2.017E+01, 2.186E+00, 0.000E+00/
      DATA (PH1(I,14,11, 2),I=1,6) /1.899E+02, 5.977E+01,1.486E+01, 2.786E+01, 4.565E+00, 0.000E+00/
      DATA (PH1(I,14,11, 3),I=1,6) /1.466E+02, 6.482E+01,2.326E+02, 1.188E+02, 3.547E+00, 8.417E-04/
      DATA (PH1(I,14,11, 4),I=1,6) /4.514E+01, 1.288E+01,6.083E+00, 1.356E+06, 3.353E+00, 0.000E+00/
      DATA (PH1(I,14,12, 1),I=1,6) /1.852E+03, 3.201E+02,3.591E+01, 2.174E+01, 2.123E+00, 0.000E+00/
      DATA (PH1(I,14,12, 2),I=1,6) /1.744E+02, 5.843E+01,1.569E+01, 2.296E+01, 4.803E+00, 0.000E+00/
      DATA (PH1(I,14,12, 3),I=1,6) /1.311E+02, 6.652E+01,2.197E+02, 1.169E+02, 3.529E+00, 8.658E-04/
      DATA (PH1(I,14,12, 4),I=1,6) /3.349E+01, 1.182E+01,8.666E+00, 3.532E+02, 4.208E+00, 0.000E+00/
      DATA (PH1(I,14,13, 1),I=1,6) /1.848E+03, 4.580E+02,1.628E+01, 7.706E+01, 1.385E+00, 0.000E+00/
      DATA (PH1(I,14,13, 2),I=1,6) /1.619E+02, 6.738E+01,1.263E+01, 3.623E+01, 4.172E+00, 0.000E+00/
      DATA (PH1(I,14,13, 3),I=1,6) /1.186E+02, 7.154E+01,1.832E+02, 3.537E+02, 3.133E+00, 2.870E-04/
      DATA (PH1(I,14,13, 4),I=1,6) /2.240E+01, 1.364E+01,8.454E+00, 7.489E+01, 4.676E+00, 0.000E+00/
      DATA (PH1(I,14,13, 5),I=1,6) /1.635E+01, 2.123E+01,6.975E+01, 4.907E+00, 9.525E+00, 3.169E-01/
      DATA (PH1(I,14,14, 1),I=1,6) /1.846E+03, 5.322E+02,1.184E+01, 2.580E+02, 1.102E+00, 0.000E+00/
      DATA (PH1(I,14,14, 2),I=1,6) /1.560E+02, 7.017E+01,1.166E+01, 4.742E+01, 3.933E+00, 0.000E+00/
      DATA (PH1(I,14,14, 3),I=1,6) /1.060E+02, 7.808E+01,1.532E+02, 5.765E+06, 2.639E+00, 2.774E-04/
      DATA (PH1(I,14,14, 4),I=1,6) /1.517E+01, 1.413E+01,1.166E+01, 2.288E+01, 5.334E+00, 0.000E+00/
      DATA (PH1(I,14,14, 5),I=1,6) /8.152E+00, 2.212E+01,1.845E+02, 3.849E+00, 9.721E+00, 2.921E-01/
      DATA (PH1(I,15, 1, 1),I=1,6) /3.070E+03, 9.701E+01,2.433E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,15, 2, 1),I=1,6) /2.817E+03, 3.381E+02,3.646E+01, 3.254E+01, 2.002E+00, 0.000E+00/
      DATA (PH1(I,15, 3, 1),I=1,6) /2.707E+03, 3.711E+02,3.220E+01, 3.257E+01, 1.869E+00, 0.000E+00/
      DATA (PH1(I,15, 3, 2),I=1,6) /6.119E+02, 6.401E+01,1.519E+01, 3.741E+01, 3.995E+00, 0.000E+00/
      DATA (PH1(I,15, 4, 1),I=1,6) /2.633E+03, 3.790E+02,3.053E+01, 3.068E+01, 1.885E+00, 0.000E+00/
      DATA (PH1(I,15, 4, 2),I=1,6) /5.604E+02, 5.712E+01,2.980E+01, 2.988E+01, 4.450E+00, 0.000E+00/
      DATA (PH1(I,15, 5, 1),I=1,6) /2.553E+03, 3.860E+02,2.914E+01, 2.656E+01, 1.945E+00, 0.000E+00/
      DATA (PH1(I,15, 5, 2),I=1,6) /5.035E+02, 7.791E+01,2.079E+01, 2.573E+01, 4.191E+00, 0.000E+00/
      DATA (PH1(I,15, 5, 3),I=1,6) /4.796E+02, 7.305E+01,7.867E+01, 7.330E+01, 3.527E+00, 1.286E-03/
      DATA (PH1(I,15, 6, 1),I=1,6) /2.477E+03, 3.824E+02,2.955E+01, 2.402E+01, 2.010E+00, 0.000E+00/
      DATA (PH1(I,15, 6, 2),I=1,6) /4.522E+02, 7.664E+01,1.966E+01, 2.403E+01, 4.317E+00, 0.000E+00/
      DATA (PH1(I,15, 6, 3),I=1,6) /4.245E+02, 6.502E+01,1.792E+02, 6.744E+01, 3.727E+00, 1.364E-03/
      DATA (PH1(I,15, 7, 1),I=1,6) /2.407E+03, 4.243E+02,2.277E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,15, 7, 2),I=1,6) /4.018E+02, 6.712E+01,1.952E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,15, 7, 3),I=1,6) /3.717E+02, 6.437E+01,2.384E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,15, 8, 1),I=1,6) /2.337E+03, 3.811E+02,2.964E+01, 2.122E+01, 2.087E+00, 0.000E+00/
      DATA (PH1(I,15, 8, 2),I=1,6) /3.552E+02, 7.217E+01,1.771E+01, 2.129E+01, 4.617E+00, 0.000E+00/
      DATA (PH1(I,15, 8, 3),I=1,6) /3.094E+02, 5.587E+01,3.696E+02, 4.738E+01, 4.244E+00, 1.094E-03/
      DATA (PH1(I,15, 9, 1),I=1,6) /2.280E+03, 4.327E+02,2.146E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,15, 9, 2),I=1,6) /3.098E+02, 7.044E+01,1.553E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,15, 9, 3),I=1,6) /2.632E+02, 6.819E+01,2.712E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,15,10, 1),I=1,6) /2.214E+03, 3.828E+02,2.932E+01, 1.991E+01, 2.123E+00, 0.000E+00/
      DATA (PH1(I,15,10, 2),I=1,6) /2.663E+02, 6.951E+01,1.560E+01, 1.845E+01, 4.938E+00, 0.000E+00/
      DATA (PH1(I,15,10, 3),I=1,6) /2.204E+02, 6.893E+01,2.980E+02, 3.497E+01, 4.307E+00, 1.008E-03/
      DATA (PH1(I,15,11, 1),I=1,6) /2.192E+03, 3.423E+02,3.671E+01, 2.010E+01, 2.232E+00, 0.000E+00/
      DATA (PH1(I,15,11, 2),I=1,6) /2.461E+02, 6.270E+01,1.535E+01, 2.486E+01, 4.787E+00, 0.000E+00/
      DATA (PH1(I,15,11, 3),I=1,6) /1.977E+02, 6.186E+01,3.368E+02, 4.700E+01, 4.256E+00, 1.118E-03/
      DATA (PH1(I,15,11, 4),I=1,6) /6.503E+01, 1.526E+00,3.688E+01, 1.111E+03, 4.890E+00, 0.000E+00/
      DATA (PH1(I,15,12, 1),I=1,6) /2.173E+03, 3.969E+02,2.640E+01, 2.945E+01, 1.901E+00, 0.000E+00/
      DATA (PH1(I,15,12, 2),I=1,6) /2.280E+02, 6.564E+01,1.491E+01, 2.333E+01, 4.787E+00, 0.000E+00/
      DATA (PH1(I,15,12, 3),I=1,6) /1.795E+02, 6.946E+01,2.661E+02, 5.760E+01, 3.967E+00, 4.482E-03/
      DATA (PH1(I,15,12, 4),I=1,6) /5.144E+01, 1.151E+01,1.153E+01, 2.303E+02, 4.447E+00, 0.000E+00/
      DATA (PH1(I,15,13, 1),I=1,6) /2.157E+03, 5.025E+02,1.577E+01, 6.163E+01, 1.467E+00, 0.000E+00/
      DATA (PH1(I,15,13, 2),I=1,6) /2.124E+02, 7.708E+01,1.206E+01, 3.388E+01, 4.192E+00, 0.000E+00/
      DATA (PH1(I,15,13, 3),I=1,6) /1.639E+02, 7.943E+01,1.949E+02, 1.392E+02, 3.391E+00, 5.039E-03/
      DATA (PH1(I,15,13, 4),I=1,6) /3.828E+01, 1.832E+01,8.260E+00, 7.058E+02, 3.682E+00, 0.000E+00/
      DATA (PH1(I,15,13, 5),I=1,6) /3.020E+01, 2.401E+01,2.851E+01, 9.795E+00, 8.210E+00, 3.290E-01/
      DATA (PH1(I,15,14, 1),I=1,6) /2.155E+03, 4.445E+02,1.990E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,15,14, 2),I=1,6) /1.987E+02, 7.299E+01,1.328E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,15,14, 3),I=1,6) /1.495E+02, 7.517E+01,2.347E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,15,14, 4),I=1,6) /2.709E+01, 1.666E+01,7.930E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,15,14, 5),I=1,6) /1.973E+01, 2.511E+01,2.564E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,15,15, 1),I=1,6) /2.154E+03, 6.472E+02,9.167E+00, 2.562E+02, 1.063E+00, 0.000E+00/
      DATA (PH1(I,15,15, 2),I=1,6) /1.940E+02, 8.632E+01,9.931E+00, 6.594E+01, 3.617E+00, 0.000E+00/
      DATA (PH1(I,15,15, 3),I=1,6) /1.400E+02, 8.812E+01,1.512E+02, 2.230E+03, 2.795E+00, 3.422E-04/
      DATA (PH1(I,15,15, 4),I=1,6) /2.017E+01, 1.658E+01,1.125E+01, 2.613E+01, 5.205E+00, 0.000E+00/
      DATA (PH1(I,15,15, 5),I=1,6) /1.049E+01, 2.580E+01,9.925E+01, 6.712E+00, 8.516E+00, 2.765E-01/
      DATA (PH1(I,16, 1, 1),I=1,6) /3.494E+03, 1.104E+02,2.139E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,16, 2, 1),I=1,6) /3.224E+03, 4.390E+02,2.453E+01, 4.405E+01, 1.765E+00, 0.000E+00/
      DATA (PH1(I,16, 3, 1),I=1,6) /3.107E+03, 4.667E+02,2.293E+01, 4.459E+01, 1.668E+00, 0.000E+00/
      DATA (PH1(I,16, 3, 2),I=1,6) /7.072E+02, 4.414E+01,2.178E+01, 4.365E+01, 4.480E+00, 0.000E+00/
      DATA (PH1(I,16, 4, 1),I=1,6) /3.029E+03, 4.669E+02,2.233E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 4, 2),I=1,6) /6.517E+02, 6.930E+01,2.571E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 5, 1),I=1,6) /2.941E+03, 4.828E+02,2.087E+01, 3.742E+01, 1.720E+00, 0.000E+00/
      DATA (PH1(I,16, 5, 2),I=1,6) /5.906E+02, 6.791E+01,2.451E+01, 2.306E+01, 4.707E+00, 0.000E+00/
      DATA (PH1(I,16, 5, 3),I=1,6) /5.647E+02, 5.542E+01,1.641E+02, 6.822E+01, 4.000E+00, 4.183E-02/
      DATA (PH1(I,16, 6, 1),I=1,6) /2.859E+03, 4.758E+02,2.101E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 6, 2),I=1,6) /5.346E+02, 7.360E+01,2.030E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 6, 3),I=1,6) /5.048E+02, 7.213E+01,1.689E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,16, 7, 1),I=1,6) /2.782E+03, 4.802E+02,2.041E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 7, 2),I=1,6) /4.804E+02, 7.555E+01,1.812E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 7, 3),I=1,6) /4.471E+02, 7.326E+01,2.213E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,16, 8, 1),I=1,6) /2.705E+03, 4.937E+02,1.946E+01, 3.568E+01, 1.737E+00, 0.000E+00/
      DATA (PH1(I,16, 8, 2),I=1,6) /4.296E+02, 8.031E+01,1.659E+01, 2.170E+01, 4.615E+00, 0.000E+00/
      DATA (PH1(I,16, 8, 3),I=1,6) /3.791E+02, 6.911E+01,2.970E+02, 5.519E+01, 4.011E+00, 2.178E-02/
      DATA (PH1(I,16, 9, 1),I=1,6) /2.641E+03, 4.891E+02,1.931E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 9, 2),I=1,6) /3.797E+02, 7.905E+01,1.462E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 9, 3),I=1,6) /3.282E+02, 7.695E+01,2.615E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,16,10, 1),I=1,6) /2.569E+03, 4.952E+02,1.916E+01, 3.555E+01, 1.742E+00, 0.000E+00/
      DATA (PH1(I,16,10, 2),I=1,6) /3.321E+02, 8.439E+01,1.291E+01, 2.335E+01, 4.545E+00, 0.000E+00/
      DATA (PH1(I,16,10, 3),I=1,6) /2.809E+02, 7.813E+01,2.751E+02, 4.495E+01, 4.110E+00, 2.108E-02/
      DATA (PH1(I,16,11, 1),I=1,6) /2.544E+03, 4.681E+02,2.148E+01, 4.092E+01, 1.738E+00, 0.000E+00/
      DATA (PH1(I,16,11, 2),I=1,6) /3.094E+02, 6.780E+01,1.495E+01, 2.519E+01, 4.826E+00, 0.000E+00/
      DATA (PH1(I,16,11, 3),I=1,6) /2.557E+02, 5.747E+01,4.712E+02, 3.610E+01, 4.742E+00, 2.480E-02/
      DATA (PH1(I,16,11, 4),I=1,6) /8.805E+01, 1.413E+01,9.139E+00, 1.656E+03, 3.626E+00, 0.000E+00/
      DATA (PH1(I,16,12, 1),I=1,6) /2.522E+03, 4.813E+02,2.029E+01, 3.854E+01, 1.736E+00, 0.000E+00/
      DATA (PH1(I,16,12, 2),I=1,6) /2.888E+02, 7.355E+01,1.407E+01, 2.397E+01, 4.754E+00, 0.000E+00/
      DATA (PH1(I,16,12, 3),I=1,6) /2.350E+02, 7.411E+01,2.919E+02, 4.864E+01, 4.142E+00, 2.785E-02/
      DATA (PH1(I,16,12, 4),I=1,6) /7.268E+01, 1.271E+01,1.363E+01, 2.570E+02, 4.361E+00, 0.000E+00/
      DATA (PH1(I,16,13, 1),I=1,6) /2.502E+03, 5.523E+02,1.504E+01, 5.527E+01, 1.517E+00, 0.000E+00/
      DATA (PH1(I,16,13, 2),I=1,6) /2.703E+02, 8.665E+01,1.155E+01, 3.202E+01, 4.227E+00, 0.000E+00/
      DATA (PH1(I,16,13, 3),I=1,6) /2.164E+02, 8.744E+01,2.037E+02, 9.310E+01, 3.565E+00, 9.497E-03/
      DATA (PH1(I,16,13, 4),I=1,6) /5.750E+01, 2.143E+01,9.080E+00, 7.224E+02, 3.618E+00, 0.000E+00/
      DATA (PH1(I,16,13, 5),I=1,6) /4.731E+01, 3.325E+01,1.099E+01, 3.639E+01, 5.977E+00, 3.761E-01/
      DATA (PH1(I,16,14, 1),I=1,6) /2.486E+03, 6.434E+02,1.082E+01, 1.410E+02, 1.221E+00, 0.000E+00/
      DATA (PH1(I,16,14, 2),I=1,6) /2.536E+02, 9.233E+01,1.048E+01, 4.007E+01, 3.968E+00, 0.000E+00/
      DATA (PH1(I,16,14, 3),I=1,6) /1.995E+02, 9.246E+01,1.780E+02, 1.498E+02, 3.319E+00, 2.142E-02/
      DATA (PH1(I,16,14, 4),I=1,6) /4.415E+01, 2.175E+01,8.208E+00, 7.941E+02, 3.591E+00, 0.000E+00/
      DATA (PH1(I,16,14, 5),I=1,6) /3.483E+01, 2.934E+01,3.702E+01, 1.445E+01, 7.321E+00, 2.989E-01/
      DATA (PH1(I,16,15, 1),I=1,6) /2.478E+03, 7.306E+02,8.303E+00, 7.483E+02, 9.844E-01, 0.000E+00/
      DATA (PH1(I,16,15, 2),I=1,6) /2.387E+02, 9.694E+01,9.659E+00, 5.188E+01, 3.738E+00, 0.000E+00/
      DATA (PH1(I,16,15, 3),I=1,6) /1.846E+02, 9.058E+01,1.896E+02, 7.538E+01, 3.635E+00, 2.934E-01/
      DATA (PH1(I,16,15, 4),I=1,6) /3.190E+01, 2.047E+01,7.824E+00, 2.396E+02, 3.950E+00, 0.000E+00/
      DATA (PH1(I,16,15, 5),I=1,6) /2.333E+01, 2.890E+01,1.054E+02, 7.474E+00, 8.421E+00, 2.840E-01/
      DATA (PH1(I,16,16, 1),I=1,6) /2.477E+03, 8.114E+02,6.649E+00, 3.734E+03, 8.646E-01, 0.000E+00/
      DATA (PH1(I,16,16, 2),I=1,6) /2.350E+02, 1.047E+02,8.520E+00, 9.469E+01, 3.346E+00, 0.000E+00/
      DATA (PH1(I,16,16, 3),I=1,6) /1.700E+02, 9.152E+01,1.883E+02, 7.193E+01, 3.633E+00, 2.485E-01/
      DATA (PH1(I,16,16, 4),I=1,6) /2.130E+01, 1.916E+01,1.003E+01, 3.296E+01, 5.038E+00, 0.000E+00/
      DATA (PH1(I,16,16, 5),I=1,6) /1.036E+01, 2.975E+01,5.644E+01, 1.321E+01, 7.513E+00, 2.621E-01/
      DATA (PH1(I,17, 1, 1),I=1,6) /3.946E+03, 1.247E+02,1.894E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,17, 2, 1),I=1,6) /3.659E+03, 4.174E+02,3.143E+01, 3.162E+01, 2.037E+00, 0.000E+00/
      DATA (PH1(I,17, 3, 1),I=1,6) /3.534E+03, 4.432E+02,2.956E+01, 3.058E+01, 1.953E+00, 0.000E+00/
      DATA (PH1(I,17, 3, 2),I=1,6) /8.094E+02, 7.444E+01,1.356E+01, 3.797E+01, 4.103E+00, 0.000E+00/
      DATA (PH1(I,17, 4, 1),I=1,6) /3.448E+03, 4.473E+02,2.878E+01, 2.772E+01, 1.999E+00, 0.000E+00/
      DATA (PH1(I,17, 4, 2),I=1,6) /7.498E+02, 7.310E+01,2.449E+01, 3.044E+01, 4.421E+00, 0.000E+00/
      DATA (PH1(I,17, 5, 1),I=1,6) /3.356E+03, 4.633E+02,2.651E+01, 2.523E+01, 2.025E+00, 0.000E+00/
      DATA (PH1(I,17, 5, 2),I=1,6) /6.846E+02, 1.049E+02,1.622E+01, 2.683E+01, 4.067E+00, 0.000E+00/
      DATA (PH1(I,17, 5, 3),I=1,6) /6.567E+02, 3.905E+01,3.845E+02, 8.984E+01, 4.272E+00, 3.017E-01/
      DATA (PH1(I,17, 6, 1),I=1,6) /3.266E+03, 4.610E+02,2.666E+01, 2.298E+01, 2.085E+00, 0.000E+00/
      DATA (PH1(I,17, 6, 2),I=1,6) /6.249E+02, 9.701E+01,1.646E+01, 2.531E+01, 4.266E+00, 0.000E+00/
      DATA (PH1(I,17, 6, 3),I=1,6) /5.920E+02, 7.621E+01,1.831E+02, 6.513E+01, 3.852E+00, 3.987E-01/
      DATA (PH1(I,17, 7, 1),I=1,6) /3.184E+03, 5.397E+02,1.840E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17, 7, 2),I=1,6) /5.660E+02, 8.450E+01,1.684E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17, 7, 3),I=1,6) /5.293E+02, 8.280E+01,2.053E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,17, 8, 1),I=1,6) /3.100E+03, 4.593E+02,2.668E+01, 2.099E+01, 2.148E+00, 0.000E+00/
      DATA (PH1(I,17, 8, 2),I=1,6) /5.109E+02, 9.134E+01,1.502E+01, 2.321E+01, 4.515E+00, 0.000E+00/
      DATA (PH1(I,17, 8, 3),I=1,6) /4.556E+02, 8.470E+01,2.378E+02, 5.655E+01, 3.889E+00, 3.329E-01/
      DATA (PH1(I,17, 9, 1),I=1,6) /3.030E+03, 5.491E+02,1.746E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17, 9, 2),I=1,6) /4.565E+02, 8.818E+01,1.375E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17, 9, 3),I=1,6) /4.001E+02, 8.635E+01,2.505E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,17,10, 1),I=1,6) /2.951E+03, 4.603E+02,2.645E+01, 1.991E+01, 2.182E+00, 0.000E+00/
      DATA (PH1(I,17,10, 2),I=1,6) /4.048E+02, 8.397E+01,1.392E+01, 2.002E+01, 4.904E+00, 0.000E+00/
      DATA (PH1(I,17,10, 3),I=1,6) /3.483E+02, 8.434E+01,2.880E+02, 4.019E+01, 4.242E+00, 2.672E-01/
      DATA (PH1(I,17,11, 1),I=1,6) /2.923E+03, 4.414E+02,2.836E+01, 2.402E+01, 2.121E+00, 0.000E+00/
      DATA (PH1(I,17,11, 2),I=1,6) /3.797E+02, 6.106E+01,1.676E+01, 2.195E+01, 5.347E+00, 0.000E+00/
      DATA (PH1(I,17,11, 3),I=1,6) /3.207E+02, 6.471E+01,4.512E+02, 3.950E+01, 4.637E+00, 2.356E-01/
      DATA (PH1(I,17,11, 4),I=1,6) /1.142E+02, 1.534E+00,6.782E+01, 1.497E+03, 4.756E+00, 0.000E+00/
      DATA (PH1(I,17,12, 1),I=1,6) /2.898E+03, 5.111E+02,2.068E+01, 3.167E+01, 1.860E+00, 0.000E+00/
      DATA (PH1(I,17,12, 2),I=1,6) /3.566E+02, 7.233E+01,1.520E+01, 2.042E+01, 5.148E+00, 0.000E+00/
      DATA (PH1(I,17,12, 3),I=1,6) /2.974E+02, 7.682E+01,3.345E+02, 3.803E+01, 4.431E+00, 2.396E-01/
      DATA (PH1(I,17,12, 4),I=1,6) /9.703E+01, 8.603E+00,2.080E+01, 2.738E+02, 4.689E+00, 0.000E+00/
      DATA (PH1(I,17,13, 1),I=1,6) /2.875E+03, 6.130E+02,1.391E+01, 5.524E+01, 1.527E+00, 0.000E+00/
      DATA (PH1(I,17,13, 2),I=1,6) /3.354E+02, 8.204E+01,1.391E+01, 2.025E+01, 4.934E+00, 0.000E+00/
      DATA (PH1(I,17,13, 3),I=1,6) /2.760E+02, 8.644E+01,2.701E+02, 3.890E+01, 4.239E+00, 2.462E-01/
      DATA (PH1(I,17,13, 4),I=1,6) /7.997E+01, 1.825E+01,1.148E+01, 1.976E+02, 4.245E+00, 0.000E+00/
      DATA (PH1(I,17,13, 5),I=1,6) /6.782E+01, 3.654E+01,1.195E+01, 3.759E+01, 6.009E+00, 4.036E-01/
      DATA (PH1(I,17,14, 1),I=1,6) /2.851E+03, 5.520E+02,1.713E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17,14, 2),I=1,6) /3.151E+02, 9.076E+01,1.198E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17,14, 3),I=1,6) /2.552E+02, 9.562E+01,2.144E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,17,14, 4),I=1,6) /6.470E+01, 2.019E+01,9.833E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,17,14, 5),I=1,6) /5.347E+01, 3.505E+01,2.148E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,17,15, 1),I=1,6) /2.838E+03, 7.650E+02,8.653E+00, 2.714E+02, 1.105E+00, 0.000E+00/
      DATA (PH1(I,17,15, 2),I=1,6) /2.979E+02, 9.151E+01,1.265E+01, 2.088E+01, 4.703E+00, 0.000E+00/
      DATA (PH1(I,17,15, 3),I=1,6) /2.382E+02, 9.580E+01,2.222E+02, 4.214E+01, 4.019E+00, 2.237E-01/
      DATA (PH1(I,17,15, 4),I=1,6) /5.019E+01, 2.307E+01,8.302E+00, 2.450E+02, 3.950E+00, 0.000E+00/
      DATA (PH1(I,17,15, 5),I=1,6) /3.961E+01, 3.460E+01,3.044E+01, 3.383E+01, 6.261E+00, 2.707E-01/
      DATA (PH1(I,17,16, 1),I=1,6) /2.832E+03, 5.619E+02,1.634E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17,16, 2),I=1,6) /2.837E+02, 9.168E+01,1.167E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17,16, 3),I=1,6) /2.236E+02, 9.695E+01,2.052E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,17,16, 4),I=1,6) /3.686E+01, 2.191E+01,7.493E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,17,16, 5),I=1,6) /2.381E+01, 3.339E+01,4.993E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,17,17, 1),I=1,6) /2.830E+03, 9.700E+02,5.255E+00, 1.856E+06, 7.888E-01, 0.000E+00/
      DATA (PH1(I,17,17, 2),I=1,6) /2.780E+02, 1.092E+02,1.059E+01, 2.491E+01, 4.205E+00, 0.000E+00/
      DATA (PH1(I,17,17, 3),I=1,6) /2.090E+02, 8.004E+01,3.053E+02, 3.498E+01, 4.457E+00, 2.017E-01/
      DATA (PH1(I,17,17, 4),I=1,6) /2.531E+01, 2.231E+01,6.628E+00, 1.843E+02, 4.196E+00, 0.000E+00/
      DATA (PH1(I,17,17, 5),I=1,6) /1.297E+01, 3.398E+01,4.539E+01, 2.232E+01, 6.896E+00, 2.479E-01/
      DATA (PH1(I,18, 1, 1),I=1,6) /4.426E+03, 1.399E+02,1.690E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,18, 2, 1),I=1,6) /4.121E+03, 4.468E+02,3.108E+01, 3.039E+01, 2.092E+00, 0.000E+00/
      DATA (PH1(I,18, 3, 1),I=1,6) /3.988E+03, 4.680E+02,3.003E+01, 2.854E+01, 2.037E+00, 0.000E+00/
      DATA (PH1(I,18, 3, 2),I=1,6) /9.180E+02, 5.440E+01,1.838E+01, 4.656E+01, 4.439E+00, 0.000E+00/
      DATA (PH1(I,18, 4, 1),I=1,6) /3.898E+03, 4.756E+02,2.883E+01, 2.615E+01, 2.074E+00, 0.000E+00/
      DATA (PH1(I,18, 4, 2),I=1,6) /8.548E+02, 8.162E+01,2.235E+01, 3.057E+01, 4.418E+00, 0.000E+00/
      DATA (PH1(I,18, 5, 1),I=1,6) /3.798E+03, 4.749E+02,2.874E+01, 2.235E+01, 2.171E+00, 0.000E+00/
      DATA (PH1(I,18, 5, 2),I=1,6) /7.856E+02, 1.086E+02,1.606E+01, 2.688E+01, 4.173E+00, 0.000E+00/
      DATA (PH1(I,18, 5, 3),I=1,6) /7.558E+02, 8.270E+01,1.006E+02, 7.201E+01, 3.789E+00, 5.509E-03/
      DATA (PH1(I,18, 6, 1),I=1,6) /3.702E+03, 4.731E+02,2.888E+01, 2.042E+01, 2.234E+00, 0.000E+00/
      DATA (PH1(I,18, 6, 2),I=1,6) /7.217E+02, 1.081E+02,1.515E+01, 2.572E+01, 4.251E+00, 0.000E+00/
      DATA (PH1(I,18, 6, 3),I=1,6) /6.861E+02, 5.722E+01,3.784E+02, 7.663E+01, 4.151E+00, 1.194E-02/
      DATA (PH1(I,18, 7, 1),I=1,6) /3.613E+03, 6.028E+02,1.666E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18, 7, 2),I=1,6) /6.584E+02, 9.398E+01,1.567E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18, 7, 3),I=1,6) /6.183E+02, 9.297E+01,1.904E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,18, 8, 1),I=1,6) /3.523E+03, 4.789E+02,2.796E+01, 1.917E+01, 2.271E+00, 0.000E+00/
      DATA (PH1(I,18, 8, 2),I=1,6) /5.992E+02, 9.884E+01,1.431E+01, 2.351E+01, 4.545E+00, 0.000E+00/
      DATA (PH1(I,18, 8, 3),I=1,6) /5.390E+02, 7.562E+01,3.466E+02, 5.329E+01, 4.195E+00, 1.392E-02/
      DATA (PH1(I,18, 9, 1),I=1,6) /3.446E+03, 6.126E+02,1.585E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18, 9, 2),I=1,6) /5.403E+02, 9.783E+01,1.294E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18, 9, 3),I=1,6) /4.787E+02, 9.639E+01,2.389E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,18,10, 1),I=1,6) /3.361E+03, 4.679E+02,2.931E+01, 1.744E+01, 2.362E+00, 0.000E+00/
      DATA (PH1(I,18,10, 2),I=1,6) /4.845E+02, 9.273E+01,1.306E+01, 2.086E+01, 4.862E+00, 0.000E+00/
      DATA (PH1(I,18,10, 3),I=1,6) /4.225E+02, 9.230E+01,2.855E+02, 4.508E+01, 4.165E+00, 8.883E-03/
      DATA (PH1(I,18,11, 1),I=1,6) /3.331E+03, 4.486E+02,3.143E+01, 2.008E+01, 2.315E+00, 0.000E+00/
      DATA (PH1(I,18,11, 2),I=1,6) /4.570E+02, 6.529E+01,1.602E+01, 2.336E+01, 5.323E+00, 0.000E+00/
      DATA (PH1(I,18,11, 3),I=1,6) /3.925E+02, 7.368E+01,4.198E+02, 4.419E+01, 4.492E+00, 7.712E-03/
      DATA (PH1(I,18,11, 4),I=1,6) /1.435E+02, 3.884E+00,3.295E+01, 7.082E+02, 4.645E+00, 0.000E+00/
      DATA (PH1(I,18,12, 1),I=1,6) /3.303E+03, 6.199E+02,1.559E+01, 4.503E+01, 1.664E+00, 0.000E+00/
      DATA (PH1(I,18,12, 2),I=1,6) /4.314E+02, 8.000E+01,1.424E+01, 2.136E+01, 5.094E+00, 0.000E+00/
      DATA (PH1(I,18,12, 3),I=1,6) /3.667E+02, 8.563E+01,3.204E+02, 4.230E+01, 4.329E+00, 7.258E-03/
      DATA (PH1(I,18,12, 4),I=1,6) /1.243E+02, 8.091E+00,2.596E+01, 3.212E+02, 4.685E+00, 0.000E+00/
      DATA (PH1(I,18,13, 1),I=1,6) /3.277E+03, 6.427E+02,1.444E+01, 4.208E+01, 1.659E+00, 0.000E+00/
      DATA (PH1(I,18,13, 2),I=1,6) /4.076E+02, 8.974E+01,1.313E+01, 2.094E+01, 4.919E+00, 0.000E+00/
      DATA (PH1(I,18,13, 3),I=1,6) /3.426E+02, 9.233E+01,2.786E+02, 4.220E+01, 4.227E+00, 7.408E-03/
      DATA (PH1(I,18,13, 4),I=1,6) /1.056E+02, 1.839E+01,1.304E+01, 2.002E+02, 4.303E+00, 0.000E+00/
      DATA (PH1(I,18,13, 5),I=1,6) /9.101E+01, 3.564E+01,1.418E+01, 3.945E+01, 6.194E+00, 7.398E-02/
      DATA (PH1(I,18,14, 1),I=1,6) /3.253E+03, 6.819E+02,1.269E+01, 4.958E+01, 1.562E+00, 0.000E+00/
      DATA (PH1(I,18,14, 2),I=1,6) /3.852E+02, 9.393E+01,1.265E+01, 2.082E+01, 4.848E+00, 0.000E+00/
      DATA (PH1(I,18,14, 3),I=1,6) /3.200E+02, 9.835E+01,2.474E+02, 4.284E+01, 4.125E+00, 7.283E-03/
      DATA (PH1(I,18,14, 4),I=1,6) /8.828E+01, 2.249E+01,1.048E+01, 2.022E+02, 4.127E+00, 0.000E+00/
      DATA (PH1(I,18,14, 5),I=1,6) /7.502E+01, 3.732E+01,2.476E+01, 3.876E+01, 6.142E+00, 1.882E-01/
      DATA (PH1(I,18,15, 1),I=1,6) /3.228E+03, 6.171E+02,1.548E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18,15, 2),I=1,6) /3.642E+02, 1.010E+02,1.120E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18,15, 3),I=1,6) /2.987E+02, 1.073E+02,2.007E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,18,15, 4),I=1,6) /7.174E+01, 2.269E+01,9.306E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,18,15, 5),I=1,6) /5.981E+01, 3.958E+01,3.058E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,18,16, 1),I=1,6) /3.216E+03, 8.408E+02,8.071E+00, 1.847E+02, 1.160E+00, 0.000E+00/
      DATA (PH1(I,18,16, 2),I=1,6) /3.455E+02, 1.031E+02,1.164E+01, 2.134E+01, 4.657E+00, 0.000E+00/
      DATA (PH1(I,18,16, 3),I=1,6) /2.801E+02, 1.025E+02,2.281E+02, 4.380E+01, 4.046E+00, 7.167E-03/
      DATA (PH1(I,18,16, 4),I=1,6) /5.637E+01, 2.612E+01,7.899E+00, 2.318E+02, 3.961E+00, 0.000E+00/
      DATA (PH1(I,18,16, 5),I=1,6) /4.074E+01, 3.862E+01,3.973E+01, 3.208E+01, 6.347E+00, 2.643E-01/
      DATA (PH1(I,18,17, 1),I=1,6) /3.208E+03, 6.260E+02,1.489E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18,17, 2),I=1,6) /3.317E+02, 1.018E+02,1.096E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18,17, 3),I=1,6) /2.662E+02, 1.085E+02,1.937E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,18,17, 4),I=1,6) /4.198E+01, 2.473E+01,7.104E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,18,17, 5),I=1,6) /2.763E+01, 3.780E+01,6.006E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,18,18, 1),I=1,6) /3.203E+03, 1.135E+03,4.280E+00, 3.285E+07, 7.631E-01, 0.000E+00/
      DATA (PH1(I,18,18, 2),I=1,6) /3.260E+02, 1.302E+02,9.185E+00, 2.693E+01, 4.021E+00, 0.000E+00/
      DATA (PH1(I,18,18, 3),I=1,6) /2.492E+02, 1.647E+02,8.372E+01, 5.452E+01, 3.328E+00, 6.270E-01/
      DATA (PH1(I,18,18, 4),I=1,6) /2.892E+01, 2.525E+01,6.394E+00, 1.700E+02, 4.223E+00, 0.000E+00/
      DATA (PH1(I,18,18, 5),I=1,6) /1.576E+01, 3.854E+01,4.872E+01, 2.640E+01, 6.662E+00, 2.355E-01/
      DATA (PH1(I,19, 1, 1),I=1,6) /4.934E+03, 1.559E+02,1.517E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,19, 2, 1),I=1,6) /4.611E+03, 5.177E+02,2.589E+01, 3.282E+01, 2.018E+00, 0.000E+00/
      DATA (PH1(I,19, 3, 1),I=1,6) /4.471E+03, 5.376E+02,2.534E+01, 3.129E+01, 1.966E+00, 0.000E+00/
      DATA (PH1(I,19, 3, 2),I=1,6) /1.035E+03, 1.155E+02,8.875E+00, 3.732E+01, 3.844E+00, 0.000E+00/
      DATA (PH1(I,19, 4, 1),I=1,6) /4.375E+03, 5.442E+02,2.451E+01, 2.858E+01, 2.005E+00, 0.000E+00/
      DATA (PH1(I,19, 4, 2),I=1,6) /9.680E+02, 8.979E+01,2.067E+01, 3.141E+01, 4.401E+00, 0.000E+00/
      DATA (PH1(I,19, 5, 1),I=1,6) /4.269E+03, 5.279E+02,2.596E+01, 2.351E+01, 2.144E+00, 0.000E+00/
      DATA (PH1(I,19, 5, 2),I=1,6) /8.935E+02, 1.199E+02,1.483E+01, 2.724E+01, 4.169E+00, 0.000E+00/
      DATA (PH1(I,19, 5, 3),I=1,6) /8.611E+02, 9.596E+01,7.963E+01, 6.368E+01, 3.884E+00, 1.487E+00/
      DATA (PH1(I,19, 6, 1),I=1,6) /4.166E+03, 5.297E+02,2.562E+01, 2.213E+01, 2.182E+00, 0.000E+00/
      DATA (PH1(I,19, 6, 2),I=1,6) /8.255E+02, 1.179E+02,1.420E+01, 2.609E+01, 4.262E+00, 0.000E+00/
      DATA (PH1(I,19, 6, 3),I=1,6) /7.867E+02, 1.054E+02,1.244E+02, 6.541E+01, 3.764E+00, 9.280E-01/
      DATA (PH1(I,19, 7, 1),I=1,6) /4.070E+03, 6.694E+02,1.515E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19, 7, 2),I=1,6) /7.578E+02, 1.040E+02,1.460E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19, 7, 3),I=1,6) /7.147E+02, 1.038E+02,1.768E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,19, 8, 1),I=1,6) /3.974E+03, 5.204E+02,2.646E+01, 1.958E+01, 2.282E+00, 0.000E+00/
      DATA (PH1(I,19, 8, 2),I=1,6) /6.945E+02, 1.112E+02,1.306E+01, 2.458E+01, 4.472E+00, 0.000E+00/
      DATA (PH1(I,19, 8, 3),I=1,6) /6.295E+02, 1.376E+02,1.164E+02, 5.816E+01, 3.609E+00, 8.872E-01/
      DATA (PH1(I,19, 9, 1),I=1,6) /3.890E+03, 6.796E+02,1.446E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19, 9, 2),I=1,6) /6.310E+02, 1.080E+02,1.218E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19, 9, 3),I=1,6) /5.647E+02, 1.071E+02,2.272E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,19,10, 1),I=1,6) /3.799E+03, 5.158E+02,2.685E+01, 1.837E+01, 2.338E+00, 0.000E+00/
      DATA (PH1(I,19,10, 2),I=1,6) /5.712E+02, 1.100E+02,1.119E+01, 2.456E+01, 4.565E+00, 0.000E+00/
      DATA (PH1(I,19,10, 3),I=1,6) /5.038E+02, 1.558E+02,1.104E+02, 5.635E+01, 3.554E+00, 7.774E-01/
      DATA (PH1(I,19,11, 1),I=1,6) /3.766E+03, 5.034E+02,2.785E+01, 2.112E+01, 2.275E+00, 0.000E+00/
      DATA (PH1(I,19,11, 2),I=1,6) /5.413E+02, 8.084E+01,1.395E+01, 2.604E+01, 4.988E+00, 0.000E+00/
      DATA (PH1(I,19,11, 3),I=1,6) /4.713E+02, 1.737E+02,8.545E+01, 4.865E+01, 3.547E+00, 8.892E-01/
      DATA (PH1(I,19,11, 4),I=1,6) /1.758E+02, 1.352E+01,1.538E+01, 7.752E+02, 3.884E+00, 0.000E+00/
      DATA (PH1(I,19,12, 1),I=1,6) /3.735E+03, 6.944E+02,1.389E+01, 4.826E+01, 1.635E+00, 0.000E+00/
      DATA (PH1(I,19,12, 2),I=1,6) /5.133E+02, 1.060E+02,1.113E+01, 2.784E+01, 4.499E+00, 0.000E+00/
      DATA (PH1(I,19,12, 3),I=1,6) /4.430E+02, 1.414E+02,1.332E+02, 5.593E+01, 3.674E+00, 7.309E-01/
      DATA (PH1(I,19,12, 4),I=1,6) /1.547E+02, 1.930E+01,1.619E+01, 3.223E+02, 4.079E+00, 0.000E+00/
      DATA (PH1(I,19,13, 1),I=1,6) /3.706E+03, 7.193E+02,1.288E+01, 4.554E+01, 1.628E+00, 0.000E+00/
      DATA (PH1(I,19,13, 2),I=1,6) /4.868E+02, 1.155E+02,1.017E+01, 2.932E+01, 4.325E+00, 0.000E+00/
      DATA (PH1(I,19,13, 3),I=1,6) /4.162E+02, 1.360E+02,1.449E+02, 5.618E+01, 3.706E+00, 6.504E-01/
      DATA (PH1(I,19,13, 4),I=1,6) /1.344E+02, 3.470E+01,9.467E+00, 1.032E+03, 3.365E+00, 0.000E+00/
      DATA (PH1(I,19,13, 5),I=1,6) /1.176E+02, 4.576E+01,1.153E+01, 6.308E+01, 5.520E+00, 2.470E-01/
      DATA (PH1(I,19,14, 1),I=1,6) /3.679E+03, 8.147E+02,9.835E+00, 9.408E+01, 1.355E+00, 0.000E+00/
      DATA (PH1(I,19,14, 2),I=1,6) /4.618E+02, 1.184E+02,9.849E+00, 2.972E+01, 4.280E+00, 0.000E+00/
      DATA (PH1(I,19,14, 3),I=1,6) /3.909E+02, 2.054E+02,5.562E+01, 1.347E+04, 2.349E+00, 5.284E-01/
      DATA (PH1(I,19,14, 4),I=1,6) /1.152E+02, 3.602E+01,8.591E+00, 2.262E+03, 3.242E+00, 0.000E+00/
      DATA (PH1(I,19,14, 5),I=1,6) /9.944E+01, 5.453E+01,1.638E+01, 9.114E+01, 5.035E+00, 4.138E-01/
      DATA (PH1(I,19,15, 1),I=1,6) /3.651E+03, 6.817E+02,1.428E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19,15, 2),I=1,6) /4.374E+02, 1.113E+02,1.061E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19,15, 3),I=1,6) /3.665E+02, 1.188E+02,1.922E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,19,15, 4),I=1,6) /9.657E+01, 2.466E+01,9.855E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,19,15, 5),I=1,6) /8.266E+01, 4.382E+01,3.214E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,19,16, 1),I=1,6) /3.633E+03, 8.442E+02,9.086E+00, 8.929E+01, 1.343E+00, 0.000E+00/
      DATA (PH1(I,19,16, 2),I=1,6) /4.165E+02, 1.249E+02,9.136E+00, 3.259E+01, 4.130E+00, 0.000E+00/
      DATA (PH1(I,19,16, 3),I=1,6) /3.452E+02, 1.876E+02,6.770E+01, 1.658E+07, 2.363E+00, 4.546E-01/
      DATA (PH1(I,19,16, 4),I=1,6) /7.925E+01, 3.296E+01,7.619E+00, 1.323E+03, 3.423E+00, 0.000E+00/
      DATA (PH1(I,19,16, 5),I=1,6) /6.091E+01, 4.791E+01,3.344E+01, 5.426E+01, 5.626E+00, 3.140E-01/
      DATA (PH1(I,19,17, 1),I=1,6) /3.623E+03, 6.897E+02,1.383E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19,17, 2),I=1,6) /3.990E+02, 1.121E+02,1.039E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19,17, 3),I=1,6) /3.279E+02, 1.199E+02,1.862E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,19,17, 4),I=1,6) /6.269E+01, 2.670E+01,7.455E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,19,17, 5),I=1,6) /4.581E+01, 4.055E+01,6.784E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,19,18, 1),I=1,6) /3.617E+03, 1.103E+03,5.173E+00, 2.564E+03, 8.917E-01, 0.000E+00/
      DATA (PH1(I,19,18, 2),I=1,6) /3.867E+02, 1.488E+02,7.085E+00, 6.561E+01, 3.449E+00, 0.000E+00/
      DATA (PH1(I,19,18, 3),I=1,6) /3.067E+02, 2.286E+02,4.199E+01, 8.731E+06, 2.239E+00, 5.560E-01/
      DATA (PH1(I,19,18, 4),I=1,6) /4.728E+01, 2.985E+01,6.625E+00, 3.001E+02, 3.891E+00, 0.000E+00/
      DATA (PH1(I,19,18, 5),I=1,6) /3.163E+01, 4.195E+01,7.392E+01, 1.697E+01, 7.208E+00, 2.531E-01/
      DATA (PH1(I,19,19, 1),I=1,6) /3.614E+03, 1.171E+03,4.540E+00, 6.165E+03, 8.392E-01, 0.000E+00/
      DATA (PH1(I,19,19, 2),I=1,6) /3.843E+02, 1.602E+02,6.389E+00, 1.044E+02, 3.159E+00, 0.000E+00/
      DATA (PH1(I,19,19, 3),I=1,6) /3.014E+02, 2.666E+02,3.107E+01, 7.187E+06, 2.067E+00, 5.274E-01/
      DATA (PH1(I,19,19, 4),I=1,6) /4.080E+01, 2.910E+01,6.377E+00, 2.229E+03, 3.587E+00, 0.000E+00/
      DATA (PH1(I,19,19, 5),I=1,6) /2.466E+01, 4.138E+01,2.614E+01, 2.143E+02, 5.631E+00, 2.437E-01/
      DATA (PH1(I,19,19, 6),I=1,6) /1.000E+00, 1.000E+00,0.000E+00, 1.000E+00, 1.000E+00, 1.000E+00/
      DATA (PH1(I,19,19, 7),I=1,6) /4.341E+00, 3.824E+00,7.363E-01, 2.410E+07, 4.427E+00, 2.049E-04/
      DATA (PH1(I,20, 1, 1),I=1,6) /5.470E+03, 1.729E+02,1.369E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,20, 2, 1),I=1,6) /5.129E+03, 6.297E+02,1.936E+01, 3.921E+01, 1.862E+00, 0.000E+00/
      DATA (PH1(I,20, 3, 1),I=1,6) /4.982E+03, 6.640E+02,1.820E+01, 3.979E+01, 1.777E+00, 0.000E+00/
      DATA (PH1(I,20, 3, 2),I=1,6) /1.157E+03, 9.395E+01,1.115E+01, 3.959E+01, 4.176E+00, 0.000E+00/
      DATA (PH1(I,20, 4, 1),I=1,6) /4.880E+03, 7.235E+02,1.486E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 4, 2),I=1,6) /1.087E+03, 1.072E+02,1.788E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 5, 1),I=1,6) /4.767E+03, 6.862E+02,1.665E+01, 3.443E+01, 1.823E+00, 0.000E+00/
      DATA (PH1(I,20, 5, 2),I=1,6) /1.008E+03, 1.169E+02,1.552E+01, 2.618E+01, 4.389E+00, 0.000E+00/
      DATA (PH1(I,20, 5, 3),I=1,6) /9.745E+02, 4.041E+01,5.215E+02, 1.044E+02, 4.452E+00, 1.164E-05/
      DATA (PH1(I,20, 6, 1),I=1,6) /4.659E+03, 7.342E+02,1.417E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 6, 2),I=1,6) /9.357E+02, 1.122E+02,1.490E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 6, 3),I=1,6) /8.946E+02, 1.144E+02,1.199E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,20, 7, 1),I=1,6) /4.555E+03, 7.396E+02,1.384E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 7, 2),I=1,6) /8.642E+02, 1.145E+02,1.363E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 7, 3),I=1,6) /8.177E+02, 1.152E+02,1.643E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,20, 8, 1),I=1,6) /4.453E+03, 6.989E+02,1.570E+01, 3.218E+01, 1.851E+00, 0.000E+00/
      DATA (PH1(I,20, 8, 2),I=1,6) /7.968E+02, 1.154E+02,1.287E+01, 2.461E+01, 4.561E+00, 0.000E+00/
      DATA (PH1(I,20, 8, 3),I=1,6) /7.267E+02, 8.942E+01,3.336E+02, 5.766E+01, 4.171E+00, 1.104E-05/
      DATA (PH1(I,20, 9, 1),I=1,6) /4.362E+03, 7.503E+02,1.324E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 9, 2),I=1,6) /7.287E+02, 1.187E+02,1.147E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 9, 3),I=1,6) /6.572E+02, 1.184E+02,2.157E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,20,10, 1),I=1,6) /4.265E+03, 7.010E+02,1.547E+01, 3.197E+01, 1.858E+00, 0.000E+00/
      DATA (PH1(I,20,10, 2),I=1,6) /6.649E+02, 9.384E+01,1.413E+01, 1.780E+01, 5.361E+00, 0.000E+00/
      DATA (PH1(I,20,10, 3),I=1,6) /5.919E+02, 7.052E+01,5.909E+02, 4.129E+01, 4.878E+00, 7.117E-06/
      DATA (PH1(I,20,11, 1),I=1,6) /4.229E+03, 6.805E+02,1.642E+01, 3.522E+01, 1.842E+00, 0.000E+00/
      DATA (PH1(I,20,11, 2),I=1,6) /6.326E+02, 8.448E+01,1.359E+01, 2.702E+01, 5.021E+00, 0.000E+00/
      DATA (PH1(I,20,11, 3),I=1,6) /5.569E+02, 6.511E+01,6.616E+02, 4.371E+01, 4.937E+00, 7.881E-06/
      DATA (PH1(I,20,11, 4),I=1,6) /2.113E+02, 1.605E+01,1.437E+01, 6.989E+02, 3.857E+00, 0.000E+00/
      DATA (PH1(I,20,12, 1),I=1,6) /4.198E+03, 7.401E+02,1.369E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20,12, 2),I=1,6) /6.018E+02, 1.203E+02,1.052E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20,12, 3),I=1,6) /5.270E+02, 1.293E+02,1.934E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,20,12, 4),I=1,6) /1.883E+02, 2.552E+01,1.335E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,20,13, 1),I=1,6) /4.163E+03, 7.816E+02,1.217E+01, 4.492E+01, 1.646E+00, 0.000E+00/
      DATA (PH1(I,20,13, 2),I=1,6) /5.732E+02, 1.243E+02,9.859E+00, 2.823E+01, 4.388E+00, 0.000E+00/
      DATA (PH1(I,20,13, 3),I=1,6) /4.967E+02, 1.070E+02,2.788E+02, 4.768E+01, 4.200E+00, 4.591E-06/
      DATA (PH1(I,20,13, 4),I=1,6) /1.664E+02, 3.989E+01,9.257E+00, 9.884E+02, 3.326E+00, 0.000E+00/
      DATA (PH1(I,20,13, 5),I=1,6) /1.472E+02, 5.202E+01,1.127E+01, 6.821E+01, 5.408E+00, 2.204E-01/
      DATA (PH1(I,20,14, 1),I=1,6) /4.133E+03, 8.261E+02,1.079E+01, 5.751E+01, 1.531E+00, 0.000E+00/
      DATA (PH1(I,20,14, 2),I=1,6) /5.455E+02, 1.268E+02,9.589E+00, 2.852E+01, 4.353E+00, 0.000E+00/
      DATA (PH1(I,20,14, 3),I=1,6) /4.687E+02, 1.208E+02,2.181E+02, 5.816E+01, 3.907E+00, 4.346E-06/
      DATA (PH1(I,20,14, 4),I=1,6) /1.452E+02, 4.100E+01,8.521E+00, 1.977E+03, 3.218E+00, 0.000E+00/
      DATA (PH1(I,20,14, 5),I=1,6) /1.272E+02, 6.097E+01,1.648E+01, 9.236E+01, 5.007E+00, 4.261E-01/
      DATA (PH1(I,20,15, 1),I=1,6) /4.105E+03, 8.505E+02,1.012E+01, 6.416E+01, 1.482E+00, 0.000E+00/
      DATA (PH1(I,20,15, 2),I=1,6) /5.193E+02, 1.297E+02,9.293E+00, 2.927E+01, 4.300E+00, 0.000E+00/
      DATA (PH1(I,20,15, 3),I=1,6) /4.423E+02, 1.303E+02,1.855E+02, 6.993E+01, 3.707E+00, 1.400E-05/
      DATA (PH1(I,20,15, 4),I=1,6) /1.249E+02, 3.994E+01,8.068E+00, 2.341E+03, 3.235E+00, 0.000E+00/
      DATA (PH1(I,20,15, 5),I=1,6) /1.088E+02, 5.749E+01,2.483E+01, 7.819E+01, 5.203E+00, 3.800E-01/
      DATA (PH1(I,20,16, 1),I=1,6) /4.078E+03, 8.635E+02,9.791E+00, 6.215E+01, 1.479E+00, 0.000E+00/
      DATA (PH1(I,20,16, 2),I=1,6) /4.948E+02, 1.315E+02,9.126E+00, 2.928E+01, 4.281E+00, 0.000E+00/
      DATA (PH1(I,20,16, 3),I=1,6) /4.175E+02, 1.384E+02,1.622E+02, 8.811E+01, 3.521E+00, 4.384E-04/
      DATA (PH1(I,20,16, 4),I=1,6) /1.054E+02, 3.722E+01,7.757E+00, 1.194E+03, 3.398E+00, 0.000E+00/
      DATA (PH1(I,20,16, 5),I=1,6) /8.451E+01, 5.673E+01,3.044E+01, 7.647E+01, 5.265E+00, 3.557E-01/
      DATA (PH1(I,20,17, 1),I=1,6) /4.063E+03, 9.606E+02,7.782E+00, 1.191E+02, 1.271E+00, 0.000E+00/
      DATA (PH1(I,20,17, 2),I=1,6) /4.719E+02, 1.373E+02,8.622E+00, 3.169E+01, 4.152E+00, 0.000E+00/
      DATA (PH1(I,20,17, 3),I=1,6) /3.944E+02, 1.413E+02,1.542E+02, 9.906E+01, 3.446E+00, 1.107E-03/
      DATA (PH1(I,20,17, 4),I=1,6) /8.680E+01, 3.645E+01,7.134E+00, 1.200E+03, 3.447E+00, 0.000E+00/
      DATA (PH1(I,20,17, 5),I=1,6) /6.727E+01, 5.169E+01,4.224E+01, 4.415E+01, 5.864E+00, 3.052E-01/
      DATA (PH1(I,20,18, 1),I=1,6) /4.053E+03, 1.003E+03,7.079E+00, 1.108E+02, 1.256E+00, 0.000E+00/
      DATA (PH1(I,20,18, 2),I=1,6) /4.542E+02, 1.483E+02,7.770E+00, 3.866E+01, 3.891E+00, 0.000E+00/
      DATA (PH1(I,20,18, 3),I=1,6) /3.731E+02, 1.260E+02,1.945E+02, 6.819E+01, 3.770E+00, 4.791E-04/
      DATA (PH1(I,20,18, 4),I=1,6) /6.920E+01, 3.504E+01,6.596E+00, 7.909E+02, 3.588E+00, 0.000E+00/
      DATA (PH1(I,20,18, 5),I=1,6) /5.091E+01, 4.617E+01,7.930E+01, 1.711E+01, 7.186E+00, 2.658E-01/
      DATA (PH1(I,20,19, 1),I=1,6) /4.047E+03, 7.997E+02,1.168E+01, 3.233E+01, 1.744E+00, 0.000E+00/
      DATA (PH1(I,20,19, 2),I=1,6) /4.445E+02, 1.335E+02,8.939E+00, 2.914E+01, 4.269E+00, 0.000E+00/
      DATA (PH1(I,20,19, 3),I=1,6) /3.638E+02, 1.448E+02,1.446E+02, 1.349E+02, 3.300E+00, 3.358E-04/
      DATA (PH1(I,20,19, 4),I=1,6) /6.037E+01, 3.176E+01,6.924E+00, 5.246E+02, 3.771E+00, 0.000E+00/
      DATA (PH1(I,20,19, 5),I=1,6) /4.090E+01, 4.498E+01,7.314E+01, 1.898E+01, 7.152E+00, 2.735E-01/
      DATA (PH1(I,20,19, 6),I=1,6) /1.000E+00, 1.000E+00,0.000E+00, 1.000E+00, 1.000E+00, 1.000E+00/
      DATA (PH1(I,20,19, 7),I=1,6) /1.187E+01, 4.155E+00,2.235E+00, 1.595E+04, 4.313E+00, 3.539E-01/
      DATA (PH1(I,20,20, 1),I=1,6) /4.043E+03, 6.947E+02,1.586E+01, 2.563E+01, 1.966E+00, 0.000E+00/
      DATA (PH1(I,20,20, 2),I=1,6) /4.425E+02, 1.201E+02,1.010E+01, 2.468E+01, 4.592E+00, 0.000E+00/
      DATA (PH1(I,20,20, 3),I=1,6) /3.523E+02, 1.529E+02,1.282E+02, 2.217E+02, 3.087E+00, 3.343E-03/
      DATA (PH1(I,20,20, 4),I=1,6) /4.830E+01, 3.012E+01,7.227E+00, 1.736E+02, 4.165E+00, 0.000E+00/
      DATA (PH1(I,20,20, 5),I=1,6) /3.443E+01, 4.487E+01,9.017E+01, 1.465E+01, 7.498E+00, 2.754E-01/
      DATA (PH1(I,20,20, 6),I=1,6) /1.000E+00, 1.000E+00,0.000E+00, 1.000E+00, 1.000E+00, 1.000E+00/
      DATA (PH1(I,20,20, 7),I=1,6) /6.113E+00, 7.366E+00,2.373E+00, 2.082E+02, 4.841E+00, 5.841E-04/
      DATA (PH1(I,21, 1, 1),I=1,6) /6.034E+03, 1.907E+02,1.241E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,21, 2, 1),I=1,6) /5.675E+03, 6.267E+02,2.182E+01, 3.452E+01, 1.998E+00, 0.000E+00/
      DATA (PH1(I,21, 3, 1),I=1,6) /5.520E+03, 6.439E+02,2.169E+01, 3.290E+01, 1.960E+00, 0.000E+00/
      DATA (PH1(I,21, 3, 2),I=1,6) /1.288E+03, 1.197E+02,8.820E+00, 3.718E+01, 4.054E+00, 0.000E+00/
      DATA (PH1(I,21, 4, 1),I=1,6) /5.413E+03, 6.607E+02,2.038E+01, 3.109E+01, 1.971E+00, 0.000E+00/
      DATA (PH1(I,21, 4, 2),I=1,6) /1.213E+03, 1.109E+02,1.725E+01, 3.167E+01, 4.371E+00, 0.000E+00/
      DATA (PH1(I,21, 5, 1),I=1,6) /5.294E+03, 6.684E+02,1.971E+01, 2.792E+01, 2.022E+00, 0.000E+00/
      DATA (PH1(I,21, 5, 2),I=1,6) /1.130E+03, 1.386E+02,1.325E+01, 2.785E+01, 4.218E+00, 0.000E+00/
      DATA (PH1(I,21, 5, 3),I=1,6) /1.094E+03, 1.435E+02,4.807E+01, 8.557E+01, 3.452E+00, 4.096E-02/
      DATA (PH1(I,21, 6, 1),I=1,6) /5.178E+03, 6.752E+02,1.918E+01, 2.616E+01, 2.052E+00, 0.000E+00/
      DATA (PH1(I,21, 6, 2),I=1,6) /1.054E+03, 1.367E+02,1.270E+01, 2.670E+01, 4.303E+00, 0.000E+00/
      DATA (PH1(I,21, 6, 3),I=1,6) /1.009E+03, 1.269E+02,1.155E+02, 8.134E+01, 3.618E+00, 1.480E-02/
      DATA (PH1(I,21, 7, 1),I=1,6) /5.067E+03, 8.133E+02,1.269E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21, 7, 2),I=1,6) /9.775E+02, 1.256E+02,1.274E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21, 7, 3),I=1,6) /9.275E+02, 1.273E+02,1.529E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,21, 8, 1),I=1,6) /4.960E+03, 6.739E+02,1.909E+01, 2.432E+01, 2.101E+00, 0.000E+00/
      DATA (PH1(I,21, 8, 2),I=1,6) /9.062E+02, 1.274E+02,1.195E+01, 2.504E+01, 4.532E+00, 0.000E+00/
      DATA (PH1(I,21, 8, 3),I=1,6) /8.308E+02, 7.899E+01,4.780E+02, 6.650E+01, 4.304E+00, 9.238E-03/
      DATA (PH1(I,21, 9, 1),I=1,6) /4.861E+03, 8.244E+02,1.216E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21, 9, 2),I=1,6) /8.333E+02, 1.300E+02,1.081E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21, 9, 3),I=1,6) /7.567E+02, 1.304E+02,2.045E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,21,10, 1),I=1,6) /4.759E+03, 6.781E+02,1.869E+01, 2.393E+01, 2.110E+00, 0.000E+00/
      DATA (PH1(I,21,10, 2),I=1,6) /7.657E+02, 1.156E+02,1.130E+01, 2.296E+01, 4.853E+00, 0.000E+00/
      DATA (PH1(I,21,10, 3),I=1,6) /6.874E+02, 8.313E+01,5.075E+02, 4.739E+01, 4.637E+00, 3.716E-03/
      DATA (PH1(I,21,11, 1),I=1,6) /4.720E+03, 6.682E+02,1.915E+01, 2.701E+01, 2.060E+00, 0.000E+00/
      DATA (PH1(I,21,11, 2),I=1,6) /7.309E+02, 7.807E+01,1.420E+01, 2.694E+01, 5.297E+00, 0.000E+00/
      DATA (PH1(I,21,11, 3),I=1,6) /6.505E+02, 8.091E+01,5.240E+02, 4.759E+01, 4.674E+00, 3.660E-03/
      DATA (PH1(I,21,11, 4),I=1,6) /2.498E+02, 3.414E+00,5.976E+01, 1.109E+03, 4.577E+00, 0.000E+00/
      DATA (PH1(I,21,12, 1),I=1,6) /4.684E+03, 7.696E+02,1.415E+01, 3.494E+01, 1.825E+00, 0.000E+00/
      DATA (PH1(I,21,12, 2),I=1,6) /6.982E+02, 1.019E+02,1.215E+01, 2.361E+01, 5.033E+00, 0.000E+00/
      DATA (PH1(I,21,12, 3),I=1,6) /6.173E+02, 1.007E+02,3.568E+02, 4.400E+01, 4.463E+00, 3.640E-03/
      DATA (PH1(I,21,12, 4),I=1,6) /2.251E+02, 1.104E+01,2.980E+01, 3.451E+02, 4.545E+00, 0.000E+00/
      DATA (PH1(I,21,13, 1),I=1,6) /4.649E+03, 8.283E+02,1.206E+01, 4.105E+01, 1.707E+00, 0.000E+00/
      DATA (PH1(I,21,13, 2),I=1,6) /6.666E+02, 1.131E+02,1.124E+01, 2.303E+01, 4.893E+00, 0.000E+00/
      DATA (PH1(I,21,13, 3),I=1,6) /5.852E+02, 1.050E+02,3.290E+02, 4.379E+01, 4.411E+00, 3.675E-03/
      DATA (PH1(I,21,13, 4),I=1,6) /2.015E+02, 1.499E+01,2.155E+01, 2.761E+02, 4.455E+00, 0.000E+00/
      DATA (PH1(I,21,13, 5),I=1,6) /1.800E+02, 4.763E+01,1.441E+01, 4.454E+01, 6.085E+00, 5.576E-03/
      DATA (PH1(I,21,14, 1),I=1,6) /4.615E+03, 9.428E+02,9.117E+00, 7.840E+01, 1.428E+00, 0.000E+00/
      DATA (PH1(I,21,14, 2),I=1,6) /6.362E+02, 1.182E+02,1.081E+01, 2.278E+01, 4.836E+00, 0.000E+00/
      DATA (PH1(I,21,14, 3),I=1,6) /5.544E+02, 1.106E+02,2.980E+02, 4.355E+01, 4.348E+00, 3.645E-03/
      DATA (PH1(I,21,14, 4),I=1,6) /1.784E+02, 2.031E+01,1.547E+01, 2.264E+02, 4.346E+00, 0.000E+00/
      DATA (PH1(I,21,14, 5),I=1,6) /1.581E+02, 4.941E+01,2.617E+01, 4.344E+01, 6.056E+00, 6.050E-03/
      DATA (PH1(I,21,15, 1),I=1,6) /4.582E+03, 8.245E+02,1.210E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21,15, 2),I=1,6) /6.052E+02, 1.335E+02,9.523E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21,15, 3),I=1,6) /5.236E+02, 1.434E+02,1.765E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,21,15, 4),I=1,6) /1.557E+02, 2.891E+01,1.052E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,21,15, 5),I=1,6) /1.380E+02, 5.322E+01,3.348E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,21,16, 1),I=1,6) /4.554E+03, 9.555E+02,8.840E+00, 7.637E+01, 1.426E+00, 0.000E+00/
      DATA (PH1(I,21,16, 2),I=1,6) /5.803E+02, 1.249E+02,1.026E+01, 2.249E+01, 4.764E+00, 0.000E+00/
      DATA (PH1(I,21,16, 3),I=1,6) /4.977E+02, 1.266E+02,2.303E+02, 4.517E+01, 4.139E+00, 3.408E-03/
      DATA (PH1(I,21,16, 4),I=1,6) /1.348E+02, 2.899E+01,9.855E+00, 1.921E+02, 4.175E+00, 0.000E+00/
      DATA (PH1(I,21,16, 5),I=1,6) /1.107E+02, 5.034E+01,4.505E+01, 3.882E+01, 6.159E+00, 6.694E-03/
      DATA (PH1(I,21,17, 1),I=1,6) /4.531E+03, 8.307E+02,1.185E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21,17, 2),I=1,6) /5.552E+02, 1.344E+02,9.326E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21,17, 3),I=1,6) /4.730E+02, 1.444E+02,1.721E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,21,17, 4),I=1,6) /1.141E+02, 3.036E+01,8.550E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,21,17, 5),I=1,6) /9.187E+01, 5.439E+01,4.585E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,21,18, 1),I=1,6) /4.517E+03, 9.957E+02,8.077E+00, 7.124E+01, 1.415E+00, 0.000E+00/
      DATA (PH1(I,21,18, 2),I=1,6) /5.306E+02, 1.371E+02,9.401E+00, 2.320E+01, 4.579E+00, 0.000E+00/
      DATA (PH1(I,21,18, 3),I=1,6) /4.476E+02, 1.422E+02,1.820E+02, 5.428E+01, 3.854E+00, 3.269E-03/
      DATA (PH1(I,21,18, 4),I=1,6) /9.451E+01, 3.451E+01,7.316E+00, 1.881E+02, 4.071E+00, 0.000E+00/
      DATA (PH1(I,21,18, 5),I=1,6) /7.349E+01, 5.468E+01,5.377E+01, 3.321E+01, 6.239E+00, 2.941E-01/
      DATA (PH1(I,21,19, 1),I=1,6) /4.508E+03, 1.030E+03,7.512E+00, 6.545E+01, 1.410E+00, 0.000E+00/
      DATA (PH1(I,21,19, 2),I=1,6) /5.148E+02, 1.443E+02,8.939E+00, 2.406E+01, 4.456E+00, 0.000E+00/
      DATA (PH1(I,21,19, 3),I=1,6) /4.285E+02, 1.348E+02,2.014E+02, 5.215E+01, 3.952E+00, 3.368E-03/
      DATA (PH1(I,21,19, 4),I=1,6) /7.762E+01, 3.513E+01,6.785E+00, 1.898E+02, 4.068E+00, 0.000E+00/
      DATA (PH1(I,21,19, 5),I=1,6) /5.591E+01, 5.169E+01,5.465E+01, 2.788E+01, 6.613E+00, 2.567E-01/
      DATA (PH1(I,21,19, 6),I=1,6) /2.476E+01, 2.183E+01,8.459E+01, 6.603E+01, 6.452E+00, 3.890E-01/
      DATA (PH1(I,21,20, 1),I=1,6) /4.500E+03, 7.669E+02,1.437E+01, 2.629E+01, 1.953E+00, 0.000E+00/
      DATA (PH1(I,21,20, 2),I=1,6) /5.054E+02, 1.274E+02,1.003E+01, 2.215E+01, 4.755E+00, 0.000E+00/
      DATA (PH1(I,21,20, 3),I=1,6) /4.190E+02, 1.472E+02,1.711E+02, 5.517E+01, 3.796E+00, 3.289E-03/
      DATA (PH1(I,21,20, 4),I=1,6) /6.848E+01, 3.246E+01,6.970E+00, 1.757E+02, 4.198E+00, 0.000E+00/
      DATA (PH1(I,21,20, 5),I=1,6) /4.678E+01, 4.990E+01,5.795E+01, 2.476E+01, 6.870E+00, 2.658E-01/
      DATA (PH1(I,21,20, 6),I=1,6) /1.444E+01, 1.611E+01,1.269E+02, 5.229E+01, 7.300E+00, 3.103E-01/
      DATA (PH1(I,21,20, 7),I=1,6) /1.280E+01, 3.837E+00,1.909E+00, 5.825E+02, 5.018E+00, 1.415E-02/
      DATA (PH1(I,21,21, 1),I=1,6) /4.494E+03, 7.367E+02,1.554E+01, 2.940E+01, 1.937E+00, 0.000E+00/
      DATA (PH1(I,21,21, 2),I=1,6) /5.032E+02, 1.180E+02,1.065E+01, 2.172E+01, 4.911E+00, 0.000E+00/
      DATA (PH1(I,21,21, 3),I=1,6) /4.054E+02, 1.593E+02,1.473E+02, 6.007E+01, 3.635E+00, 3.208E-03/
      DATA (PH1(I,21,21, 4),I=1,6) /5.640E+01, 3.284E+01,6.849E+00, 1.709E+02, 4.207E+00, 0.000E+00/
      DATA (PH1(I,21,21, 5),I=1,6) /3.360E+01, 4.936E+01,6.176E+01, 2.186E+01, 7.081E+00, 2.671E-01/
      DATA (PH1(I,21,21, 6),I=1,6) /8.010E+00, 9.733E+00,2.488E+02, 2.066E+01, 1.022E+01, 3.104E-01/
      DATA (PH1(I,21,21, 7),I=1,6) /7.342E+00, 8.324E+00,2.252E+00, 4.118E+02, 4.588E+00, 1.441E-04/
      DATA (PH1(I,22, 1, 1),I=1,6) /6.626E+03, 2.094E+02,1.131E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,22, 2, 1),I=1,6) /6.249E+03, 8.656E+02,1.234E+01, 5.893E+01, 1.622E+00, 0.000E+00/
      DATA (PH1(I,22, 3, 1),I=1,6) /6.087E+03, 8.854E+02,1.240E+01, 8.225E+01, 1.477E+00, 0.000E+00/
      DATA (PH1(I,22, 3, 2),I=1,6) /1.425E+03, 1.025E+02,1.038E+01, 4.131E+01, 4.265E+00, 0.000E+00/
      DATA (PH1(I,22, 4, 1),I=1,6) /5.974E+03, 7.415E+02,1.773E+01, 3.385E+01, 1.915E+00, 0.000E+00/
      DATA (PH1(I,22, 4, 2),I=1,6) /1.346E+03, 1.388E+02,1.408E+01, 3.133E+01, 4.198E+00, 0.000E+00/
      DATA (PH1(I,22, 5, 1),I=1,6) /5.848E+03, 7.545E+02,1.694E+01, 3.060E+01, 1.955E+00, 0.000E+00/
      DATA (PH1(I,22, 5, 2),I=1,6) /1.260E+03, 1.433E+02,1.296E+01, 2.842E+01, 4.278E+00, 0.000E+00/
      DATA (PH1(I,22, 5, 3),I=1,6) /1.221E+03, 1.686E+02,3.856E+01, 9.198E+01, 3.355E+00, 3.372E-05/
      DATA (PH1(I,22, 6, 1),I=1,6) /5.726E+03, 7.617E+02,1.648E+01, 2.934E+01, 1.972E+00, 0.000E+00/
      DATA (PH1(I,22, 6, 2),I=1,6) /1.179E+03, 1.419E+02,1.240E+01, 2.693E+01, 4.373E+00, 0.000E+00/
      DATA (PH1(I,22, 6, 3),I=1,6) /1.131E+03, 1.291E+02,1.268E+02, 8.686E+01, 3.655E+00, 1.321E-05/
      DATA (PH1(I,22, 7, 1),I=1,6) /5.607E+03, 8.905E+02,1.167E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22, 7, 2),I=1,6) /1.098E+03, 1.372E+02,1.192E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22, 7, 3),I=1,6) /1.044E+03, 1.401E+02,1.425E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,22, 8, 1),I=1,6) /5.495E+03, 7.690E+02,1.599E+01, 2.789E+01, 1.997E+00, 0.000E+00/
      DATA (PH1(I,22, 8, 2),I=1,6) /1.022E+03, 1.341E+02,1.154E+01, 2.533E+01, 4.576E+00, 0.000E+00/
      DATA (PH1(I,22, 8, 3),I=1,6) /9.419E+02, 1.084E+02,2.972E+02, 6.609E+01, 4.056E+00, 1.465E-05/
      DATA (PH1(I,22, 9, 1),I=1,6) /5.387E+03, 9.021E+02,1.121E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22, 9, 2),I=1,6) /9.448E+02, 1.417E+02,1.019E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22, 9, 3),I=1,6) /8.631E+02, 1.430E+02,1.937E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,22,10, 1),I=1,6) /5.281E+03, 7.712E+02,1.576E+01, 2.765E+01, 2.005E+00, 0.000E+00/
      DATA (PH1(I,22,10, 2),I=1,6) /8.734E+02, 1.250E+02,1.070E+01, 2.364E+01, 4.830E+00, 0.000E+00/
      DATA (PH1(I,22,10, 3),I=1,6) /7.878E+02, 9.852E+01,4.256E+02, 4.955E+01, 4.485E+00, 1.750E-05/
      DATA (PH1(I,22,11, 1),I=1,6) /5.239E+03, 7.394E+02,1.717E+01, 2.884E+01, 2.020E+00, 0.000E+00/
      DATA (PH1(I,22,11, 2),I=1,6) /8.363E+02, 9.240E+01,1.291E+01, 2.687E+01, 5.160E+00, 0.000E+00/
      DATA (PH1(I,22,11, 3),I=1,6) /7.511E+02, 9.404E+01,4.550E+02, 5.003E+01, 4.543E+00, 1.753E-05/
      DATA (PH1(I,22,11, 4),I=1,6) /2.915E+02, 3.359E+00,6.875E+01, 1.244E+03, 4.559E+00, 0.000E+00/
      DATA (PH1(I,22,12, 1),I=1,6) /5.200E+03, 9.318E+02,1.045E+01, 5.288E+01, 1.604E+00, 0.000E+00/
      DATA (PH1(I,22,12, 2),I=1,6) /8.013E+02, 1.259E+02,1.049E+01, 2.445E+01, 4.779E+00, 0.000E+00/
      DATA (PH1(I,22,12, 3),I=1,6) /7.156E+02, 1.122E+02,3.306E+02, 4.701E+01, 4.376E+00, 1.686E-05/
      DATA (PH1(I,22,12, 4),I=1,6) /2.650E+02, 9.778E+00,3.726E+01, 4.225E+02, 4.552E+00, 0.000E+00/
      DATA (PH1(I,22,13, 1),I=1,6) /5.162E+03, 9.469E+02,1.008E+01, 5.127E+01, 1.603E+00, 0.000E+00/
      DATA (PH1(I,22,13, 2),I=1,6) /7.670E+02, 1.254E+02,1.046E+01, 2.368E+01, 4.830E+00, 0.000E+00/
      DATA (PH1(I,22,13, 3),I=1,6) /6.807E+02, 1.133E+02,3.222E+02, 4.654E+01, 4.373E+00, 1.686E-05/
      DATA (PH1(I,22,13, 4),I=1,6) /2.397E+02, 1.712E+01,2.117E+01, 2.738E+02, 4.407E+00, 0.000E+00/
      DATA (PH1(I,22,13, 5),I=1,6) /2.159E+02, 5.265E+01,1.405E+01, 4.648E+01, 6.025E+00, 8.018E-03/
      DATA (PH1(I,22,14, 1),I=1,6) /5.126E+03, 9.454E+02,1.011E+01, 5.151E+01, 1.603E+00, 0.000E+00/
      DATA (PH1(I,22,14, 2),I=1,6) /7.342E+02, 1.245E+02,1.045E+01, 2.319E+01, 4.873E+00, 0.000E+00/
      DATA (PH1(I,22,14, 3),I=1,6) /6.473E+02, 1.160E+02,3.065E+02, 4.592E+01, 4.355E+00, 1.680E-05/
      DATA (PH1(I,22,14, 4),I=1,6) /2.148E+02, 1.736E+01,1.941E+01, 2.691E+02, 4.420E+00, 0.000E+00/
      DATA (PH1(I,22,14, 5),I=1,6) /1.921E+02, 5.337E+01,2.611E+01, 4.514E+01, 6.042E+00, 7.441E-03/
      DATA (PH1(I,22,15, 1),I=1,6) /5.089E+03, 9.027E+02,1.114E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22,15, 2),I=1,6) /6.999E+02, 1.454E+02,9.019E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22,15, 3),I=1,6) /6.129E+02, 1.565E+02,1.694E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,22,15, 4),I=1,6) /1.900E+02, 3.119E+01,1.069E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,22,15, 5),I=1,6) /1.704E+02, 5.839E+01,3.348E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,22,16, 1),I=1,6) /5.059E+03, 9.613E+02,9.746E+00, 5.327E+01, 1.581E+00, 0.000E+00/
      DATA (PH1(I,22,16, 2),I=1,6) /6.729E+02, 1.313E+02,9.948E+00, 2.269E+01, 4.817E+00, 0.000E+00/
      DATA (PH1(I,22,16, 3),I=1,6) /5.853E+02, 1.296E+02,2.477E+02, 4.576E+01, 4.216E+00, 1.795E-05/
      DATA (PH1(I,22,16, 4),I=1,6) /1.673E+02, 2.786E+01,1.099E+01, 1.974E+02, 4.259E+00, 0.000E+00/
      DATA (PH1(I,22,16, 5),I=1,6) /1.408E+02, 5.484E+01,4.527E+01, 3.999E+01, 6.137E+00, 9.213E-03/
      DATA (PH1(I,22,17, 1),I=1,6) /5.027E+03, 9.080E+02,1.096E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22,17, 2),I=1,6) /6.441E+02, 1.464E+02,8.835E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22,17, 3),I=1,6) /5.563E+02, 1.574E+02,1.656E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,22,17, 4),I=1,6) /1.445E+02, 3.258E+01,8.817E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,22,17, 5),I=1,6) /1.195E+02, 5.940E+01,4.665E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,22,18, 1),I=1,6) /5.011E+03, 1.031E+03,8.364E+00, 6.468E+01, 1.478E+00, 0.000E+00/
      DATA (PH1(I,22,18, 2),I=1,6) /6.178E+02, 1.446E+02,9.107E+00, 2.319E+01, 4.635E+00, 0.000E+00/
      DATA (PH1(I,22,18, 3),I=1,6) /5.296E+02, 1.466E+02,1.952E+02, 4.874E+01, 4.006E+00, 1.732E-05/
      DATA (PH1(I,22,18, 4),I=1,6) /1.231E+02, 3.545E+01,7.729E+00, 1.792E+02, 4.142E+00, 0.000E+00/
      DATA (PH1(I,22,18, 5),I=1,6) /9.930E+01, 5.902E+01,5.520E+01, 3.434E+01, 6.226E+00, 2.921E-01/
      DATA (PH1(I,22,19, 1),I=1,6) /4.998E+03, 1.048E+03,8.085E+00, 6.252E+01, 1.476E+00, 0.000E+00/
      DATA (PH1(I,22,19, 2),I=1,6) /5.949E+02, 1.452E+02,9.082E+00, 2.289E+01, 4.645E+00, 0.000E+00/
      DATA (PH1(I,22,19, 3),I=1,6) /5.066E+02, 1.497E+02,1.878E+02, 4.955E+01, 3.967E+00, 1.717E-05/
      DATA (PH1(I,22,19, 4),I=1,6) /1.035E+02, 3.700E+01,6.988E+00, 1.778E+02, 4.126E+00, 0.000E+00/
      DATA (PH1(I,22,19, 5),I=1,6) /7.917E+01, 5.872E+01,5.237E+01, 3.131E+01, 6.379E+00, 2.839E-01/
      DATA (PH1(I,22,19, 6),I=1,6) /4.327E+01, 2.624E+01,9.279E+01, 6.851E+01, 6.300E+00, 3.297E-01/
      DATA (PH1(I,22,20, 1),I=1,6) /4.988E+03, 1.086E+03,7.494E+00, 5.795E+01, 1.469E+00, 0.000E+00/
      DATA (PH1(I,22,20, 2),I=1,6) /5.792E+02, 1.542E+02,8.585E+00, 2.370E+01, 4.508E+00, 0.000E+00/
      DATA (PH1(I,22,20, 3),I=1,6) /4.870E+02, 1.508E+02,1.834E+02, 5.660E+01, 3.868E+00, 1.691E-05/
      DATA (PH1(I,22,20, 4),I=1,6) /8.603E+01, 3.786E+01,6.489E+00, 1.768E+02, 4.124E+00, 0.000E+00/
      DATA (PH1(I,22,20, 5),I=1,6) /6.201E+01, 5.643E+01,5.284E+01, 2.686E+01, 6.696E+00, 2.571E-01/
      DATA (PH1(I,22,20, 6),I=1,6) /2.749E+01, 2.399E+01,1.731E+02, 5.474E+01, 6.742E+00, 3.650E-01/
      DATA (PH1(I,22,21, 1),I=1,6) /4.980E+03, 6.950E+02,1.983E+01, 1.966E+01, 2.290E+00, 0.000E+00/
      DATA (PH1(I,22,21, 2),I=1,6) /5.704E+02, 1.312E+02,9.886E+00, 2.185E+01, 4.872E+00, 0.000E+00/
      DATA (PH1(I,22,21, 3),I=1,6) /4.772E+02, 1.630E+02,1.584E+02, 5.914E+01, 3.739E+00, 1.643E-05/
      DATA (PH1(I,22,21, 4),I=1,6) /7.654E+01, 3.518E+01,6.620E+00, 1.617E+02, 4.260E+00, 0.000E+00/
      DATA (PH1(I,22,21, 5),I=1,6) /5.253E+01, 5.494E+01,5.561E+01, 2.416E+01, 6.917E+00, 2.642E-01/
      DATA (PH1(I,22,21, 6),I=1,6) /1.613E+01, 1.771E+01,2.597E+02, 3.909E+01, 7.825E+00, 2.982E-01/
      DATA (PH1(I,22,21, 7),I=1,6) /1.358E+01, 5.011E+00,1.639E+00, 5.095E+02, 4.958E+00, 1.633E-03/
      DATA (PH1(I,22,22, 1),I=1,6) /4.972E+03, 6.826E+02,2.029E+01, 2.415E+01, 2.187E+00, 0.000E+00/
      DATA (PH1(I,22,22, 2),I=1,6) /5.690E+02, 1.177E+02,1.068E+01, 2.144E+01, 5.085E+00, 0.000E+00/
      DATA (PH1(I,22,22, 3),I=1,6) /4.640E+02, 1.814E+02,1.291E+02, 6.629E+01, 3.531E+00, 5.519E-05/
      DATA (PH1(I,22,22, 4),I=1,6) /6.500E+01, 3.562E+01,6.517E+00, 1.572E+02, 4.268E+00, 0.000E+00/
      DATA (PH1(I,22,22, 5),I=1,6) /4.000E+01, 5.424E+01,5.924E+01, 2.151E+01, 7.123E+00, 2.645E-01/
      DATA (PH1(I,22,22, 6),I=1,6) /9.940E+00, 1.102E+01,5.478E+02, 1.610E+01, 1.096E+01, 2.972E-01/
      DATA (PH1(I,22,22, 7),I=1,6) /6.820E+00, 9.184E+00,2.167E+00, 4.297E+02, 4.552E+00, 3.612E-03/
      DATA (PH1(I,23, 1, 1),I=1,6) /7.246E+03, 2.290E+02,1.035E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,23, 2, 1),I=1,6) /6.852E+03, 7.669E+02,1.757E+01, 3.884E+01, 1.928E+00, 0.000E+00/
      DATA (PH1(I,23, 3, 1),I=1,6) /6.682E+03, 7.973E+02,1.694E+01, 3.805E+01, 1.874E+00, 0.000E+00/
      DATA (PH1(I,23, 3, 2),I=1,6) /1.570E+03, 4.820E+01,1.999E+01, 7.258E+01, 4.724E+00, 0.000E+00/
      DATA (PH1(I,23, 4, 1),I=1,6) /6.563E+03, 8.146E+02,1.606E+01, 3.587E+01, 1.887E+00, 0.000E+00/
      DATA (PH1(I,23, 4, 2),I=1,6) /1.487E+03, 1.341E+02,1.463E+01, 3.226E+01, 4.336E+00, 0.000E+00/
      DATA (PH1(I,23, 5, 1),I=1,6) /6.431E+03, 8.301E+02,1.528E+01, 3.284E+01, 1.918E+00, 0.000E+00/
      DATA (PH1(I,23, 5, 2),I=1,6) /1.396E+03, 1.597E+02,1.178E+01, 2.849E+01, 4.246E+00, 0.000E+00/
      DATA (PH1(I,23, 5, 3),I=1,6) /1.355E+03, 1.329E+02,7.374E+01, 1.017E+02, 3.596E+00, 2.681E-02/
      DATA (PH1(I,23, 6, 1),I=1,6) /6.302E+03, 8.379E+02,1.487E+01, 3.151E+01, 1.935E+00, 0.000E+00/
      DATA (PH1(I,23, 6, 2),I=1,6) /1.311E+03, 1.437E+02,1.233E+01, 2.743E+01, 4.461E+00, 0.000E+00/
      DATA (PH1(I,23, 6, 3),I=1,6) /1.260E+03, 1.761E+02,7.417E+01, 9.459E+01, 3.402E+00, 1.814E-02/
      DATA (PH1(I,23, 7, 1),I=1,6) /6.174E+03, 9.714E+02,1.078E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23, 7, 2),I=1,6) /1.225E+03, 1.493E+02,1.118E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23, 7, 3),I=1,6) /1.168E+03, 1.535E+02,1.330E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,23, 8, 1),I=1,6) /6.058E+03, 8.454E+02,1.444E+01, 3.020E+01, 1.955E+00, 0.000E+00/
      DATA (PH1(I,23, 8, 2),I=1,6) /1.146E+03, 1.434E+02,1.098E+01, 2.603E+01, 4.573E+00, 0.000E+00/
      DATA (PH1(I,23, 8, 3),I=1,6) /1.060E+03, 1.412E+02,1.992E+02, 6.905E+01, 3.834E+00, 1.640E-02/
      DATA (PH1(I,23, 9, 1),I=1,6) /5.941E+03, 9.834E+02,1.036E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23, 9, 2),I=1,6) /1.063E+03, 1.540E+02,9.623E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23, 9, 3),I=1,6) /9.758E+02, 1.562E+02,1.835E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,23,10, 1),I=1,6) /5.831E+03, 8.491E+02,1.420E+01, 2.950E+01, 1.968E+00, 0.000E+00/
      DATA (PH1(I,23,10, 2),I=1,6) /9.883E+02, 1.272E+02,1.058E+01, 2.434E+01, 4.902E+00, 0.000E+00/
      DATA (PH1(I,23,10, 3),I=1,6) /8.960E+02, 1.198E+02,3.365E+02, 5.275E+01, 4.289E+00, 1.591E-02/
      DATA (PH1(I,23,11, 1),I=1,6) /5.787E+03, 8.254E+02,1.504E+01, 3.096E+01, 1.969E+00, 0.000E+00/
      DATA (PH1(I,23,11, 2),I=1,6) /9.488E+02, 9.010E+01,1.292E+01, 2.883E+01, 5.244E+00, 0.000E+00/
      DATA (PH1(I,23,11, 3),I=1,6) /8.589E+02, 1.178E+02,3.418E+02, 5.198E+01, 4.329E+00, 1.528E-02/
      DATA (PH1(I,23,11, 4),I=1,6) /3.363E+02, 4.802E+00,5.149E+01, 9.708E+02, 4.521E+00, 0.000E+00/
      DATA (PH1(I,23,12, 1),I=1,6) /5.745E+03, 8.791E+02,1.315E+01, 3.287E+01, 1.891E+00, 0.000E+00/
      DATA (PH1(I,23,12, 2),I=1,6) /9.114E+02, 1.194E+02,1.089E+01, 2.492E+01, 4.980E+00, 0.000E+00/
      DATA (PH1(I,23,12, 3),I=1,6) /8.209E+02, 1.210E+02,3.219E+02, 4.969E+01, 4.340E+00, 1.495E-02/
      DATA (PH1(I,23,12, 4),I=1,6) /3.081E+02, 1.475E+01,2.755E+01, 3.277E+02, 4.457E+00, 0.000E+00/
      DATA (PH1(I,23,13, 1),I=1,6) /5.704E+03, 9.381E+02,1.144E+01, 3.678E+01, 1.793E+00, 0.000E+00/
      DATA (PH1(I,23,13, 2),I=1,6) /8.746E+02, 1.297E+02,1.022E+01, 2.418E+01, 4.886E+00, 0.000E+00/
      DATA (PH1(I,23,13, 3),I=1,6) /7.834E+02, 1.219E+02,3.155E+02, 4.959E+01, 4.332E+00, 1.503E-02/
      DATA (PH1(I,23,13, 4),I=1,6) /2.810E+02, 1.396E+01,2.795E+01, 3.475E+02, 4.463E+00, 0.000E+00/
      DATA (PH1(I,23,13, 5),I=1,6) /2.557E+02, 6.549E+01,1.189E+01, 4.523E+01, 5.865E+00, 5.604E-01/
      DATA (PH1(I,23,14, 1),I=1,6) /5.664E+03, 9.448E+02,1.126E+01, 3.775E+01, 1.777E+00, 0.000E+00/
      DATA (PH1(I,23,14, 2),I=1,6) /8.391E+02, 1.492E+02,8.956E+00, 2.627E+01, 4.576E+00, 0.000E+00/
      DATA (PH1(I,23,14, 3),I=1,6) /7.474E+02, 1.218E+02,3.136E+02, 4.700E+01, 4.379E+00, 1.470E-02/
      DATA (PH1(I,23,14, 4),I=1,6) /2.543E+02, 5.439E+01,8.345E+00, 9.136E+02, 3.280E+00, 0.000E+00/
      DATA (PH1(I,23,14, 5),I=1,6) /2.305E+02, 7.687E+01,1.769E+01, 7.556E+01, 5.183E+00, 4.070E-01/
      DATA (PH1(I,23,15, 1),I=1,6) /5.624E+03, 9.854E+02,1.025E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23,15, 2),I=1,6) /8.017E+02, 1.579E+02,8.542E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23,15, 3),I=1,6) /7.093E+02, 1.701E+02,1.628E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,23,15, 4),I=1,6) /2.274E+02, 3.357E+01,1.077E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,23,15, 5),I=1,6) /2.058E+02, 6.386E+01,3.315E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,23,16, 1),I=1,6) /5.591E+03, 9.640E+02,1.077E+01, 4.005E+01, 1.740E+00, 0.000E+00/
      DATA (PH1(I,23,16, 2),I=1,6) /7.727E+02, 1.513E+02,8.749E+00, 2.570E+01, 4.587E+00, 0.000E+00/
      DATA (PH1(I,23,16, 3),I=1,6) /6.800E+02, 1.400E+02,2.393E+02, 5.044E+01, 4.150E+00, 1.006E-02/
      DATA (PH1(I,23,16, 4),I=1,6) /2.031E+02, 5.013E+01,7.844E+00, 8.057E+02, 3.393E+00, 0.000E+00/
      DATA (PH1(I,23,16, 5),I=1,6) /1.735E+02, 6.959E+01,3.541E+01, 6.440E+01, 5.460E+00, 2.830E-01/
      DATA (PH1(I,23,17, 1),I=1,6) /5.557E+03, 9.898E+02,1.013E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23,17, 2),I=1,6) /7.401E+02, 1.589E+02,8.371E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23,17, 3),I=1,6) /6.468E+02, 1.710E+02,1.594E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,23,17, 4),I=1,6) /1.781E+02, 3.490E+01,9.007E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,23,17, 5),I=1,6) /1.506E+02, 6.472E+01,4.696E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,23,18, 1),I=1,6) /5.531E+03, 9.782E+02,1.044E+01, 3.841E+01, 1.745E+00, 0.000E+00/
      DATA (PH1(I,23,18, 2),I=1,6) /7.123E+02, 1.580E+02,8.337E+00, 2.583E+01, 4.520E+00, 0.000E+00/
      DATA (PH1(I,23,18, 3),I=1,6) /6.188E+02, 1.692E+02,1.626E+02, 6.988E+01, 3.712E+00, 1.180E-02/
      DATA (PH1(I,23,18, 4),I=1,6) /1.549E+02, 4.722E+01,6.956E+00, 5.988E+02, 3.553E+00, 0.000E+00/
      DATA (PH1(I,23,18, 5),I=1,6) /1.281E+02, 6.736E+01,4.931E+01, 4.579E+01, 5.859E+00, 3.123E-01/
      DATA (PH1(I,23,19, 1),I=1,6) /5.516E+03, 9.875E+02,1.024E+01, 3.726E+01, 1.748E+00, 0.000E+00/
      DATA (PH1(I,23,19, 2),I=1,6) /6.859E+02, 1.633E+02,8.009E+00, 2.702E+01, 4.427E+00, 0.000E+00/
      DATA (PH1(I,23,19, 3),I=1,6) /5.923E+02, 1.721E+02,1.571E+02, 6.951E+01, 3.697E+00, 1.225E-02/
      DATA (PH1(I,23,19, 4),I=1,6) /1.329E+02, 4.591E+01,6.536E+00, 5.301E+02, 3.633E+00, 0.000E+00/
      DATA (PH1(I,23,19, 5),I=1,6) /1.059E+02, 6.053E+01,5.762E+01, 2.883E+01, 6.578E+00, 2.283E-01/
      DATA (PH1(I,23,19, 6),I=1,6) /6.528E+01, 3.190E+01,8.329E+01, 1.369E+02, 5.672E+00, 1.237E-01/
      DATA (PH1(I,23,20, 1),I=1,6) /5.504E+03, 9.965E+02,9.948E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23,20, 2),I=1,6) /6.654E+02, 1.599E+02,8.204E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23,20, 3),I=1,6) /5.715E+02, 1.723E+02,1.553E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,23,20, 4),I=1,6) /1.118E+02, 3.793E+01,6.694E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,23,20, 5),I=1,6) /8.555E+01, 6.751E+01,4.255E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,23,20, 6),I=1,6) /4.671E+01, 2.794E+01,2.071E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,23,21, 1),I=1,6) /5.493E+03, 1.056E+03,8.926E+00, 3.290E+01, 1.738E+00, 0.000E+00/
      DATA (PH1(I,23,21, 2),I=1,6) /6.421E+02, 1.795E+02,7.208E+00, 3.031E+01, 4.183E+00, 0.000E+00/
      DATA (PH1(I,23,21, 3),I=1,6) /5.487E+02, 1.685E+02,1.634E+02, 6.966E+01, 3.723E+00, 1.158E-02/
      DATA (PH1(I,23,21, 4),I=1,6) /9.452E+01, 4.339E+01,5.872E+00, 3.580E+02, 3.847E+00, 0.000E+00/
      DATA (PH1(I,23,21, 5),I=1,6) /6.813E+01, 6.101E+01,5.383E+01, 2.353E+01, 6.904E+00, 2.586E-01/
      DATA (PH1(I,23,21, 6),I=1,6) /2.931E+01, 1.913E+01,3.556E+03, 5.763E+00, 1.187E+01, 5.113E-01/
      DATA (PH1(I,23,22, 1),I=1,6) /5.484E+03, 1.467E+03,4.369E+00, 3.858E+02, 1.042E+00, 0.000E+00/
      DATA (PH1(I,23,22, 2),I=1,6) /6.389E+02, 1.378E+02,9.535E+00, 2.195E+01, 4.927E+00, 0.000E+00/
      DATA (PH1(I,23,22, 3),I=1,6) /5.292E+02, 1.302E+02,2.706E+02, 4.385E+01, 4.365E+00, 6.490E-05/
      DATA (PH1(I,23,22, 4),I=1,6) /7.942E+01, 4.123E+01,5.887E+00, 1.634E+02, 4.188E+00, 0.000E+00/
      DATA (PH1(I,23,22, 5),I=1,6) /5.323E+01, 6.049E+01,5.264E+01, 2.216E+01, 7.058E+00, 2.583E-01/
      DATA (PH1(I,23,22, 6),I=1,6) /1.466E+01, 1.451E+01,7.154E+02, 1.749E+01, 1.038E+01, 4.665E-01/
      DATA (PH1(I,23,23, 1),I=1,6) /5.475E+03, 6.550E+02,2.420E+01, 2.343E+01, 2.326E+00, 0.000E+00/
      DATA (PH1(I,23,23, 2),I=1,6) /6.380E+02, 9.688E+01,1.212E+01, 1.905E+01, 5.757E+00, 0.000E+00/
      DATA (PH1(I,23,23, 3),I=1,6) /5.270E+02, 2.290E+02,8.513E+01, 5.037E+02, 2.767E+00, 9.964E-04/
      DATA (PH1(I,23,23, 4),I=1,6) /7.700E+01, 3.668E+01,6.690E+00, 7.759E+01, 4.706E+00, 0.000E+00/
      DATA (PH1(I,23,23, 5),I=1,6) /4.700E+01, 5.937E+01,7.599E+01, 1.488E+01, 7.586E+00, 2.690E-01/
      DATA (PH1(I,23,23, 6),I=1,6) /1.200E+01, 1.146E+01,3.134E+03, 7.037E+00, 1.377E+01, 3.417E-01/
      DATA (PH1(I,23,23, 7),I=1,6) /6.740E+00, 1.002E+01,2.059E+00, 3.914E+02, 4.565E+00, 5.057E-04/
      DATA (PH1(I,24, 1, 1),I=1,6) /7.895E+03, 2.495E+02,9.505E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,24, 2, 1),I=1,6) /7.482E+03, 9.695E+02,1.185E+01, 5.351E+01, 1.696E+00, 0.000E+00/
      DATA (PH1(I,24, 3, 1),I=1,6) /7.306E+03, 9.824E+02,1.199E+01, 5.185E+01, 1.668E+00, 0.000E+00/
      DATA (PH1(I,24, 3, 2),I=1,6) /1.721E+03, 1.480E+02,7.312E+00, 3.880E+01, 4.083E+00, 0.000E+00/
      DATA (PH1(I,24, 4, 1),I=1,6) /7.181E+03, 9.952E+02,1.156E+01, 4.841E+01, 1.687E+00, 0.000E+00/
      DATA (PH1(I,24, 4, 2),I=1,6) /1.634E+03, 1.488E+02,1.333E+01, 3.225E+01, 4.310E+00, 0.000E+00/
      DATA (PH1(I,24, 5, 1),I=1,6) /7.042E+03, 1.003E+03,1.125E+01, 4.457E+01, 1.717E+00, 0.000E+00/
      DATA (PH1(I,24, 5, 2),I=1,6) /1.539E+03, 1.748E+02,1.089E+01, 2.892E+01, 4.221E+00, 0.000E+00/
      DATA (PH1(I,24, 5, 3),I=1,6) /1.497E+03, 2.216E+02,2.692E+01, 1.091E+02, 3.193E+00, 3.540E-03/
      DATA (PH1(I,24, 6, 1),I=1,6) /6.907E+03, 1.004E+03,1.114E+01, 4.238E+01, 1.740E+00, 0.000E+00/
      DATA (PH1(I,24, 6, 2),I=1,6) /1.450E+03, 1.749E+02,1.035E+01, 2.790E+01, 4.278E+00, 0.000E+00/
      DATA (PH1(I,24, 6, 3),I=1,6) /1.396E+03, 2.006E+02,6.290E+01, 1.016E+02, 3.330E+00, 1.634E-03/
      DATA (PH1(I,24, 7, 1),I=1,6) /6.769E+03, 1.056E+03,9.975E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24, 7, 2),I=1,6) /1.359E+03, 1.619E+02,1.050E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24, 7, 3),I=1,6) /1.299E+03, 1.675E+02,1.244E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,24, 8, 1),I=1,6) /6.650E+03, 1.010E+03,1.086E+01, 4.093E+01, 1.755E+00, 0.000E+00/
      DATA (PH1(I,24, 8, 2),I=1,6) /1.276E+03, 1.514E+02,1.053E+01, 2.659E+01, 4.591E+00, 0.000E+00/
      DATA (PH1(I,24, 8, 3),I=1,6) /1.185E+03, 1.297E+02,2.646E+02, 7.449E+01, 3.960E+00, 6.187E-03/
      DATA (PH1(I,24, 9, 1),I=1,6) /6.523E+03, 1.068E+03,9.610E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24, 9, 2),I=1,6) /1.189E+03, 1.668E+02,9.095E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24, 9, 3),I=1,6) /1.097E+03, 1.701E+02,1.738E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,24,10, 1),I=1,6) /6.409E+03, 1.017E+03,1.063E+01, 4.093E+01, 1.756E+00, 0.000E+00/
      DATA (PH1(I,24,10, 2),I=1,6) /1.110E+03, 1.388E+02,9.935E+00, 2.496E+01, 4.858E+00, 0.000E+00/
      DATA (PH1(I,24,10, 3),I=1,6) /1.011E+03, 1.397E+02,2.826E+02, 5.599E+01, 4.153E+00, 5.775E-03/
      DATA (PH1(I,24,11, 1),I=1,6) /6.362E+03, 9.927E+02,1.118E+01, 4.335E+01, 1.752E+00, 0.000E+00/
      DATA (PH1(I,24,11, 2),I=1,6) /1.068E+03, 9.748E+01,1.227E+01, 2.968E+01, 5.204E+00, 0.000E+00/
      DATA (PH1(I,24,11, 3),I=1,6) /9.738E+02, 1.312E+02,3.124E+02, 5.494E+01, 4.251E+00, 6.239E-03/
      DATA (PH1(I,24,11, 4),I=1,6) /3.842E+02, 5.856E+00,4.533E+01, 8.787E+02, 4.495E+00, 0.000E+00/
      DATA (PH1(I,24,12, 1),I=1,6) /6.318E+03, 1.103E+03,8.914E+00, 5.657E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,24,12, 2),I=1,6) /1.029E+03, 1.295E+02,1.028E+01, 2.555E+01, 4.945E+00, 0.000E+00/
      DATA (PH1(I,24,12, 3),I=1,6) /9.334E+02, 1.427E+02,2.652E+02, 5.321E+01, 4.178E+00, 5.905E-03/
      DATA (PH1(I,24,12, 4),I=1,6) /3.548E+02, 1.135E+01,3.894E+01, 4.475E+02, 4.496E+00, 0.000E+00/
      DATA (PH1(I,24,13, 1),I=1,6) /6.274E+03, 1.115E+03,8.708E+00, 5.540E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,24,13, 2),I=1,6) /9.893E+02, 1.454E+02,9.421E+00, 2.472E+01, 4.803E+00, 0.000E+00/
      DATA (PH1(I,24,13, 3),I=1,6) /8.933E+02, 1.424E+02,2.644E+02, 5.317E+01, 4.181E+00, 5.962E-03/
      DATA (PH1(I,24,13, 4),I=1,6) /3.255E+02, 1.164E+01,3.651E+01, 4.413E+02, 4.490E+00, 0.000E+00/
      DATA (PH1(I,24,13, 5),I=1,6) /2.981E+02, 6.675E+01,1.246E+01, 4.963E+01, 5.860E+00, 3.057E-01/
      DATA (PH1(I,24,14, 1),I=1,6) /6.231E+03, 1.111E+03,8.763E+00, 5.576E+01, 1.589E+00, 0.000E+00/
      DATA (PH1(I,24,14, 2),I=1,6) /9.512E+02, 1.447E+02,9.388E+00, 2.428E+01, 4.837E+00, 0.000E+00/
      DATA (PH1(I,24,14, 3),I=1,6) /8.546E+02, 1.426E+02,2.618E+02, 5.270E+01, 4.188E+00, 5.989E-03/
      DATA (PH1(I,24,14, 4),I=1,6) /2.969E+02, 2.539E+01,1.632E+01, 2.419E+02, 4.286E+00, 0.000E+00/
      DATA (PH1(I,24,14, 5),I=1,6) /2.708E+02, 6.544E+01,2.411E+01, 4.849E+01, 5.923E+00, 2.753E-01/
      DATA (PH1(I,24,15, 1),I=1,6) /6.187E+03, 1.073E+03,9.448E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,15, 2),I=1,6) /9.108E+02, 1.710E+02,8.093E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,15, 3),I=1,6) /8.130E+02, 1.843E+02,1.566E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,15, 4),I=1,6) /2.680E+02, 3.606E+01,1.078E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,15, 5),I=1,6) /2.444E+02, 6.965E+01,3.257E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,24,16, 1),I=1,6) /6.152E+03, 1.110E+03,8.774E+00, 5.647E+01, 1.586E+00, 0.000E+00/
      DATA (PH1(I,24,16, 2),I=1,6) /8.797E+02, 1.476E+02,9.146E+00, 2.358E+01, 4.847E+00, 0.000E+00/
      DATA (PH1(I,24,16, 3),I=1,6) /7.819E+02, 1.427E+02,2.573E+02, 5.042E+01, 4.228E+00, 5.960E-03/
      DATA (PH1(I,24,16, 4),I=1,6) /2.419E+02, 2.465E+01,1.419E+01, 2.341E+02, 4.369E+00, 0.000E+00/
      DATA (PH1(I,24,16, 5),I=1,6) /2.093E+02, 6.644E+01,4.311E+01, 4.198E+01, 6.057E+00, 3.135E-01/
      DATA (PH1(I,24,17, 1),I=1,6) /6.114E+03, 1.076E+03,9.366E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,17, 2),I=1,6) /8.432E+02, 1.720E+02,7.934E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,17, 3),I=1,6) /7.444E+02, 1.850E+02,1.537E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,17, 4),I=1,6) /2.148E+02, 3.733E+01,9.131E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,17, 5),I=1,6) /1.847E+02, 7.035E+01,4.687E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,24,18, 1),I=1,6) /6.080E+03, 1.122E+03,8.562E+00, 5.603E+01, 1.582E+00, 0.000E+00/
      DATA (PH1(I,24,18, 2),I=1,6) /8.140E+02, 1.553E+02,8.726E+00, 2.320E+01, 4.789E+00, 0.000E+00/
      DATA (PH1(I,24,18, 3),I=1,6) /7.153E+02, 1.571E+02,2.130E+02, 5.017E+01, 4.114E+00, 5.706E-03/
      DATA (PH1(I,24,18, 4),I=1,6) /1.900E+02, 3.531E+01,8.881E+00, 1.815E+02, 4.264E+00, 0.000E+00/
      DATA (PH1(I,24,18, 5),I=1,6) /1.602E+02, 6.754E+01,5.675E+01, 3.729E+01, 6.191E+00, 2.358E-01/
      DATA (PH1(I,24,19, 1),I=1,6) /6.062E+03, 1.128E+03,8.464E+00, 5.569E+01, 1.580E+00, 0.000E+00/
      DATA (PH1(I,24,19, 2),I=1,6) /7.842E+02, 1.551E+02,8.723E+00, 2.287E+01, 4.810E+00, 0.000E+00/
      DATA (PH1(I,24,19, 3),I=1,6) /6.855E+02, 1.640E+02,1.963E+02, 5.019E+01, 4.060E+00, 5.648E-03/
      DATA (PH1(I,24,19, 4),I=1,6) /1.656E+02, 3.830E+01,7.706E+00, 1.719E+02, 4.245E+00, 0.000E+00/
      DATA (PH1(I,24,19, 5),I=1,6) /1.359E+02, 6.693E+01,5.368E+01, 3.420E+01, 6.344E+00, 2.242E-01/
      DATA (PH1(I,24,19, 6),I=1,6) /9.064E+01, 3.098E+01,1.201E+02, 6.811E+01, 6.392E+00, 3.113E-03/
      DATA (PH1(I,24,20, 1),I=1,6) /6.046E+03, 1.081E+03,9.244E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,20, 2),I=1,6) /7.598E+02, 1.731E+02,7.776E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,20, 3),I=1,6) /6.602E+02, 1.862E+02,1.500E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,20, 4),I=1,6) /1.427E+02, 4.027E+01,6.892E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,20, 5),I=1,6) /1.134E+02, 7.290E+01,4.327E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,24,20, 6),I=1,6) /6.946E+01, 3.105E+01,2.269E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,24,21, 1),I=1,6) /6.033E+03, 1.224E+03,7.106E+00, 7.179E+01, 1.455E+00, 0.000E+00/
      DATA (PH1(I,24,21, 2),I=1,6) /7.331E+02, 1.640E+02,8.333E+00, 2.286E+01, 4.721E+00, 0.000E+00/
      DATA (PH1(I,24,21, 3),I=1,6) /6.344E+02, 1.721E+02,1.795E+02, 5.067E+01, 3.993E+00, 5.538E-03/
      DATA (PH1(I,24,21, 4),I=1,6) /1.219E+02, 4.256E+01,6.252E+00, 1.605E+02, 4.224E+00, 0.000E+00/
      DATA (PH1(I,24,21, 5),I=1,6) /9.275E+01, 6.776E+01,4.856E+01, 2.817E+01, 6.628E+00, 2.673E-01/
      DATA (PH1(I,24,21, 6),I=1,6) /4.916E+01, 3.483E+01,2.224E+02, 5.448E+01, 6.484E+00, 3.436E-01/
      DATA (PH1(I,24,22, 1),I=1,6) /6.021E+03, 1.085E+03,9.165E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,22, 2),I=1,6) /7.161E+02, 1.735E+02,7.720E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,22, 3),I=1,6) /6.164E+02, 1.871E+02,1.480E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,22, 4),I=1,6) /1.032E+02, 4.293E+01,5.797E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,22, 5),I=1,6) /7.436E+01, 6.450E+01,5.388E+01, 2.100E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,24,22, 6),I=1,6) /3.096E+01, 9.947E+00,1.924E+03, 9.000E+00, 1.450E+01, 3.500E-01/
      DATA (PH1(I,24,23, 1),I=1,6) /6.009E+03, 1.497E+03,4.593E+00, 1.465E+02, 1.189E+00, 0.000E+00/
      DATA (PH1(I,24,23, 2),I=1,6) /7.074E+02, 2.013E+02,6.869E+00, 2.572E+01, 4.266E+00, 0.000E+00/
      DATA (PH1(I,24,23, 3),I=1,6) /5.970E+02, 1.477E+02,2.382E+02, 4.633E+01, 4.261E+00, 6.156E-03/
      DATA (PH1(I,24,23, 4),I=1,6) /8.745E+01, 4.462E+01,5.566E+00, 1.542E+02, 4.230E+00, 0.000E+00/
      DATA (PH1(I,24,23, 5),I=1,6) /5.879E+01, 6.599E+01,5.021E+01, 2.175E+01, 7.102E+00, 2.576E-01/
      DATA (PH1(I,24,23, 6),I=1,6) /1.650E+01, 1.521E+01,1.031E+03, 1.415E+01, 1.123E+01, 4.827E-01/
      DATA (PH1(I,24,24, 1),I=1,6) /5.996E+03, 1.035E+03,1.025E+01, 3.343E+01, 1.822E+00, 0.000E+00/
      DATA (PH1(I,24,24, 2),I=1,6) /7.030E+02, 1.588E+02,8.555E+00, 2.258E+01, 4.789E+00, 0.000E+00/
      DATA (PH1(I,24,24, 3),I=1,6) /5.850E+02, 1.865E+02,1.540E+02, 5.560E+01, 3.823E+00, 5.785E-03/
      DATA (PH1(I,24,24, 4),I=1,6) /7.900E+01, 4.234E+01,5.602E+00, 1.356E+02, 4.374E+00, 0.000E+00/
      DATA (PH1(I,24,24, 5),I=1,6) /4.900E+01, 6.567E+01,5.313E+01, 1.981E+01, 7.258E+00, 2.601E-01/
      DATA (PH1(I,24,24, 6),I=1,6) /8.660E+00, 7.244E+00,1.485E+03, 9.671E+00, 1.575E+01, 7.760E-01/
      DATA (PH1(I,24,24, 7),I=1,6) /6.767E+00, 9.636E+00,6.532E-01, 5.232E+02, 4.641E+00, 9.332E-05/
      DATA (PH1(I,25, 1, 1),I=1,6) /8.572E+03, 2.709E+02,8.759E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,25, 2, 1),I=1,6) /8.141E+03, 9.661E+02,1.309E+01, 4.829E+01, 1.794E+00, 0.000E+00/
      DATA (PH1(I,25, 3, 1),I=1,6) /7.957E+03, 9.997E+02,1.266E+01, 4.727E+01, 1.749E+00, 0.000E+00/
      DATA (PH1(I,25, 3, 2),I=1,6) /1.880E+03, 2.004E+02,5.352E+00, 3.934E+01, 3.800E+00, 0.000E+00/
      DATA (PH1(I,25, 4, 1),I=1,6) /7.827E+03, 9.916E+02,1.277E+01, 4.256E+01, 1.799E+00, 0.000E+00/
      DATA (PH1(I,25, 4, 2),I=1,6) /1.788E+03, 1.631E+02,1.229E+01, 3.227E+01, 4.294E+00, 0.000E+00/
      DATA (PH1(I,25, 5, 1),I=1,6) /7.682E+03, 9.946E+02,1.255E+01, 3.824E+01, 1.847E+00, 0.000E+00/
      DATA (PH1(I,25, 5, 2),I=1,6) /1.689E+03, 1.952E+02,9.881E+00, 2.837E+01, 4.199E+00, 0.000E+00/
      DATA (PH1(I,25, 5, 3),I=1,6) /1.644E+03, 2.576E+02,2.159E+01, 1.234E+02, 3.086E+00, 1.723E-02/
      DATA (PH1(I,25, 6, 1),I=1,6) /7.539E+03, 9.973E+02,1.238E+01, 3.663E+01, 1.868E+00, 0.000E+00/
      DATA (PH1(I,25, 6, 2),I=1,6) /1.596E+03, 1.879E+02,9.744E+00, 2.825E+01, 4.277E+00, 0.000E+00/
      DATA (PH1(I,25, 6, 3),I=1,6) /1.539E+03, 2.403E+02,4.744E+01, 1.147E+02, 3.190E+00, 1.237E-02/
      DATA (PH1(I,25, 7, 1),I=1,6) /7.391E+03, 1.144E+03,9.260E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25, 7, 2),I=1,6) /1.500E+03, 1.751E+02,9.871E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25, 7, 3),I=1,6) /1.437E+03, 1.822E+02,1.165E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,25, 8, 1),I=1,6) /7.270E+03, 1.005E+03,1.205E+01, 3.464E+01, 1.893E+00, 0.000E+00/
      DATA (PH1(I,25, 8, 2),I=1,6) /1.414E+03, 1.452E+02,1.087E+01, 2.737E+01, 4.741E+00, 0.000E+00/
      DATA (PH1(I,25, 8, 3),I=1,6) /1.317E+03, 1.889E+02,1.381E+02, 7.915E+01, 3.619E+00, 1.713E-02/
      DATA (PH1(I,25, 9, 1),I=1,6) /7.132E+03, 1.157E+03,8.934E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25, 9, 2),I=1,6) /1.321E+03, 1.802E+02,8.605E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25, 9, 3),I=1,6) /1.224E+03, 1.847E+02,1.647E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,25,10, 1),I=1,6) /7.016E+03, 1.007E+03,1.191E+01, 3.374E+01, 1.907E+00, 0.000E+00/
      DATA (PH1(I,25,10, 2),I=1,6) /1.239E+03, 1.468E+02,9.533E+00, 2.564E+01, 4.859E+00, 0.000E+00/
      DATA (PH1(I,25,10, 3),I=1,6) /1.133E+03, 1.588E+02,2.466E+02, 5.979E+01, 4.044E+00, 1.415E-02/
      DATA (PH1(I,25,11, 1),I=1,6) /6.966E+03, 9.952E+02,1.220E+01, 3.635E+01, 1.885E+00, 0.000E+00/
      DATA (PH1(I,25,11, 2),I=1,6) /1.195E+03, 1.084E+02,1.148E+01, 3.003E+01, 5.139E+00, 0.000E+00/
      DATA (PH1(I,25,11, 3),I=1,6) /1.096E+03, 1.576E+02,2.464E+02, 5.814E+01, 4.081E+00, 1.417E-02/
      DATA (PH1(I,25,11, 4),I=1,6) /4.352E+02, 5.053E+00,5.856E+01, 1.103E+03, 4.490E+00, 0.000E+00/
      DATA (PH1(I,25,12, 1),I=1,6) /6.919E+03, 1.016E+03,1.168E+01, 3.428E+01, 1.894E+00, 0.000E+00/
      DATA (PH1(I,25,12, 2),I=1,6) /1.153E+03, 1.412E+02,9.670E+00, 2.617E+01, 4.900E+00, 0.000E+00/
      DATA (PH1(I,25,12, 3),I=1,6) /1.053E+03, 1.572E+02,2.455E+02, 5.660E+01, 4.109E+00, 1.411E-02/
      DATA (PH1(I,25,12, 4),I=1,6) /4.030E+02, 1.345E+01,3.535E+01, 4.180E+02, 4.462E+00, 0.000E+00/
      DATA (PH1(I,25,13, 1),I=1,6) /6.872E+03, 1.025E+03,1.145E+01, 3.313E+01, 1.902E+00, 0.000E+00/
      DATA (PH1(I,25,13, 2),I=1,6) /1.111E+03, 1.437E+02,9.475E+00, 2.535E+01, 4.916E+00, 0.000E+00/
      DATA (PH1(I,25,13, 3),I=1,6) /1.010E+03, 1.534E+02,2.551E+02, 5.638E+01, 4.140E+00, 1.436E-02/
      DATA (PH1(I,25,13, 4),I=1,6) /3.731E+02, 2.286E+01,2.034E+01, 2.772E+02, 4.326E+00, 0.000E+00/
      DATA (PH1(I,25,13, 5),I=1,6) /3.436E+02, 7.288E+01,1.205E+01, 5.254E+01, 5.787E+00, 1.633E-04/
      DATA (PH1(I,25,14, 1),I=1,6) /6.826E+03, 1.018E+03,1.161E+01, 3.333E+01, 1.905E+00, 0.000E+00/
      DATA (PH1(I,25,14, 2),I=1,6) /1.070E+03, 1.522E+02,9.043E+00, 2.474E+01, 4.856E+00, 0.000E+00/
      DATA (PH1(I,25,14, 3),I=1,6) /9.690E+02, 1.541E+02,2.509E+02, 5.596E+01, 4.141E+00, 1.435E-02/
      DATA (PH1(I,25,14, 4),I=1,6) /3.427E+02, 2.273E+01,1.932E+01, 2.771E+02, 4.340E+00, 0.000E+00/
      DATA (PH1(I,25,14, 5),I=1,6) /3.144E+02, 7.442E+01,2.225E+01, 4.974E+01, 5.829E+00, 4.047E-01/
      DATA (PH1(I,25,15, 1),I=1,6) /6.778E+03, 1.164E+03,8.718E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,15, 2),I=1,6) /1.027E+03, 1.846E+02,7.670E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,15, 3),I=1,6) /9.238E+02, 1.990E+02,1.508E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,15, 4),I=1,6) /3.118E+02, 3.864E+01,1.074E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,15, 5),I=1,6) /2.860E+02, 7.574E+01,3.181E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,16, 1),I=1,6) /6.741E+03, 1.013E+03,1.172E+01, 3.345E+01, 1.909E+00, 0.000E+00/
      DATA (PH1(I,25,16, 2),I=1,6) /9.939E+02, 1.554E+02,8.799E+00, 2.403E+01, 4.866E+00, 0.000E+00/
      DATA (PH1(I,25,16, 3),I=1,6) /8.910E+02, 1.555E+02,2.430E+02, 5.379E+01, 4.166E+00, 1.402E-02/
      DATA (PH1(I,25,16, 4),I=1,6) /2.840E+02, 2.364E+01,1.585E+01, 2.573E+02, 4.393E+00, 0.000E+00/
      DATA (PH1(I,25,16, 5),I=1,6) /2.483E+02, 7.277E+01,4.169E+01, 4.324E+01, 6.012E+00, 3.835E-01/
      DATA (PH1(I,25,17, 1),I=1,6) /6.698E+03, 1.167E+03,8.667E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,17, 2),I=1,6) /9.536E+02, 1.856E+02,7.522E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,17, 3),I=1,6) /8.492E+02, 1.997E+02,1.482E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,17, 4),I=1,6) /2.547E+02, 3.986E+01,9.199E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,17, 5),I=1,6) /2.218E+02, 7.629E+01,4.645E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,18, 1),I=1,6) /6.664E+03, 1.021E+03,1.152E+01, 3.256E+01, 1.914E+00, 0.000E+00/
      DATA (PH1(I,25,18, 2),I=1,6) /9.230E+02, 1.679E+02,8.168E+00, 2.438E+01, 4.732E+00, 0.000E+00/
      DATA (PH1(I,25,18, 3),I=1,6) /8.191E+02, 1.748E+02,1.921E+02, 5.809E+01, 3.973E+00, 1.717E-02/
      DATA (PH1(I,25,18, 4),I=1,6) /2.282E+02, 5.334E+01,7.236E+00, 4.454E+02, 3.624E+00, 0.000E+00/
      DATA (PH1(I,25,18, 5),I=1,6) /1.945E+02, 7.818E+01,4.992E+01, 4.763E+01, 5.832E+00, 3.093E-01/
      DATA (PH1(I,25,19, 1),I=1,6) /6.635E+03, 1.021E+03,1.154E+01, 3.219E+01, 1.919E+00, 0.000E+00/
      DATA (PH1(I,25,19, 2),I=1,6) /8.900E+02, 1.668E+02,8.196E+00, 2.394E+01, 4.767E+00, 0.000E+00/
      DATA (PH1(I,25,19, 3),I=1,6) /7.860E+02, 1.739E+02,1.939E+02, 5.345E+01, 4.041E+00, 1.659E-02/
      DATA (PH1(I,25,19, 4),I=1,6) /2.016E+02, 5.184E+01,6.791E+00, 3.811E+02, 3.719E+00, 0.000E+00/
      DATA (PH1(I,25,19, 5),I=1,6) /1.691E+02, 7.618E+01,4.935E+01, 3.879E+01, 6.116E+00, 3.087E-01/
      DATA (PH1(I,25,19, 6),I=1,6) /1.193E+02, 3.865E+01,9.713E+01, 1.133E+02, 5.798E+00, 1.437E-01/
      DATA (PH1(I,25,20, 1),I=1,6) /6.617E+03, 1.171E+03,8.591E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,20, 2),I=1,6) /8.614E+02, 1.867E+02,7.373E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,20, 3),I=1,6) /7.560E+02, 2.007E+02,1.450E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,20, 4),I=1,6) /1.768E+02, 4.271E+01,7.047E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,20, 5),I=1,6) /1.444E+02, 7.861E+01,4.367E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,20, 6),I=1,6) /9.575E+01, 3.452E+01,2.355E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,25,21, 1),I=1,6) /6.602E+03, 1.045E+03,1.101E+01, 3.201E+01, 1.902E+00, 0.000E+00/
      DATA (PH1(I,25,21, 2),I=1,6) /8.319E+02, 1.782E+02,7.660E+00, 2.515E+01, 4.612E+00, 0.000E+00/
      DATA (PH1(I,25,21, 3),I=1,6) /7.280E+02, 1.909E+02,1.610E+02, 6.064E+01, 3.847E+00, 2.553E-02/
      DATA (PH1(I,25,21, 4),I=1,6) /1.530E+02, 5.019E+01,5.951E+00, 3.117E+02, 3.866E+00, 0.000E+00/
      DATA (PH1(I,25,21, 5),I=1,6) /1.210E+02, 6.892E+01,5.293E+01, 2.594E+01, 6.855E+00, 1.836E-01/
      DATA (PH1(I,25,21, 6),I=1,6) /7.240E+01, 3.798E+01,2.289E+02, 7.678E+01, 6.206E+00, 6.108E-02/
      DATA (PH1(I,25,22, 1),I=1,6) /6.589E+03, 1.174E+03,8.541E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,22, 2),I=1,6) /8.118E+02, 1.872E+02,7.318E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,22, 3),I=1,6) /7.063E+02, 2.015E+02,1.433E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,22, 4),I=1,6) /1.300E+02, 4.531E+01,5.863E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,22, 5),I=1,6) /9.879E+01, 8.115E+01,3.665E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,22, 6),I=1,6) /5.120E+01, 3.615E+01,3.188E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,25,23, 1),I=1,6) /6.576E+03, 1.154E+03,8.979E+00, 2.899E+01, 1.859E+00, 0.000E+00/
      DATA (PH1(I,25,23, 2),I=1,6) /7.947E+02, 2.016E+02,6.777E+00, 2.756E+01, 4.332E+00, 0.000E+00/
      DATA (PH1(I,25,23, 3),I=1,6) /6.822E+02, 1.908E+02,1.609E+02, 6.198E+01, 3.835E+00, 2.549E-02/
      DATA (PH1(I,25,23, 4),I=1,6) /1.120E+02, 4.956E+01,5.285E+00, 2.544E+02, 4.005E+00, 0.000E+00/
      DATA (PH1(I,25,23, 5),I=1,6) /8.062E+01, 7.190E+01,4.707E+01, 2.356E+01, 6.977E+00, 2.565E-01/
      DATA (PH1(I,25,23, 6),I=1,6) /3.367E+01, 1.604E+01,1.991E+05, 2.538E+00, 1.808E+01, 6.536E-01/
      DATA (PH1(I,25,24, 1),I=1,6) /6.564E+03, 8.684E+02,1.627E+01, 2.405E+01, 2.207E+00, 0.000E+00/
      DATA (PH1(I,25,24, 2),I=1,6) /7.849E+02, 1.346E+02,9.823E+00, 1.970E+01, 5.386E+00, 0.000E+00/
      DATA (PH1(I,25,24, 3),I=1,6) /6.714E+02, 2.207E+02,1.185E+02, 9.212E+01, 3.460E+00, 1.974E-01/
      DATA (PH1(I,25,24, 4),I=1,6) /1.015E+02, 4.359E+01,5.735E+00, 1.016E+02, 4.563E+00, 0.000E+00/
      DATA (PH1(I,25,24, 5),I=1,6) /7.009E+01, 7.141E+01,4.984E+01, 2.130E+01, 7.133E+00, 2.616E-01/
      DATA (PH1(I,25,24, 6),I=1,6) /2.058E+01, 1.624E+01,2.592E+05, 2.411E+00, 1.797E+01, 4.774E-01/
      DATA (PH1(I,25,24, 7),I=1,6) /1.564E+01, 9.306E+00,1.189E+00, 1.439E+09, 4.095E+00, 8.740E-04/
      DATA (PH1(I,25,25, 1),I=1,6) /6.550E+03, 8.758E+02,1.592E+01, 3.965E+01, 1.947E+00, 0.000E+00/
      DATA (PH1(I,25,25, 2),I=1,6) /7.816E+02, 8.316E+01,1.156E+01, 2.187E+01, 6.149E+00, 0.000E+00/
      DATA (PH1(I,25,25, 3),I=1,6) /6.554E+02, 2.737E+02,7.498E+01, 3.952E+02, 2.793E+00, 4.661E-02/
      DATA (PH1(I,25,25, 4),I=1,6) /9.460E+01, 4.114E+01,6.172E+00, 5.928E+01, 4.989E+00, 0.000E+00/
      DATA (PH1(I,25,25, 5),I=1,6) /5.940E+01, 7.040E+01,6.697E+01, 1.485E+01, 7.643E+00, 2.659E-01/
      DATA (PH1(I,25,25, 6),I=1,6) /1.430E+01, 1.311E+01,1.668E+04, 4.497E+00, 1.646E+01, 3.881E-01/
      DATA (PH1(I,25,25, 7),I=1,6) /7.434E+00, 1.183E+01,1.549E+00, 2.920E+06, 4.113E+00, 3.256E-02/
      DATA (PH1(I,26, 1, 1),I=1,6) /9.278E+03, 2.932E+02,8.099E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,26, 2, 1),I=1,6) /8.829E+03, 1.057E+03,1.195E+01, 5.769E+01, 1.718E+00, 0.000E+00/
      DATA (PH1(I,26, 3, 1),I=1,6) /8.638E+03, 1.087E+03,1.157E+01, 5.086E+01, 1.722E+00, 0.000E+00/
      DATA (PH1(I,26, 3, 2),I=1,6) /2.046E+03, 1.873E+02,5.833E+00, 3.849E+01, 3.998E+00, 0.000E+00/
      DATA (PH1(I,26, 4, 1),I=1,6) /8.484E+03, 1.215E+03,9.098E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 4, 2),I=1,6) /1.950E+03, 1.799E+02,1.138E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 5, 1),I=1,6) /8.350E+03, 1.066E+03,1.181E+01, 4.116E+01, 1.827E+00, 0.000E+00/
      DATA (PH1(I,26, 5, 2),I=1,6) /1.847E+03, 1.708E+02,1.125E+01, 2.839E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 5, 3),I=1,6) /1.799E+03, 8.070E+01,2.675E+02, 1.191E+02, 4.194E+00, 1.901E-05/
      DATA (PH1(I,26, 6, 1),I=1,6) /8.184E+03, 1.228E+03,8.773E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 6, 2),I=1,6) /1.745E+03, 1.860E+02,9.942E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 6, 3),I=1,6) /1.689E+03, 1.970E+02,7.720E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,26, 7, 1),I=1,6) /8.041E+03, 1.235E+03,8.619E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 7, 2),I=1,6) /1.648E+03, 1.888E+02,9.298E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 7, 3),I=1,6) /1.582E+03, 1.975E+02,1.093E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,26, 8, 1),I=1,6) /7.918E+03, 1.068E+03,1.155E+01, 3.578E+01, 1.895E+00, 0.000E+00/
      DATA (PH1(I,26, 8, 2),I=1,6) /1.559E+03, 1.623E+02,1.000E+01, 2.675E+01, 4.711E+00, 0.000E+00/
      DATA (PH1(I,26, 8, 3),I=1,6) /1.456E+03, 1.338E+02,3.000E+02, 7.010E+01, 4.143E+00, 1.631E-05/
      DATA (PH1(I,26, 9, 1),I=1,6) /7.769E+03, 1.249E+03,8.327E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 9, 2),I=1,6) /1.460E+03, 1.940E+02,8.150E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 9, 3),I=1,6) /1.358E+03, 1.999E+02,1.561E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,26,10, 1),I=1,6) /7.651E+03, 1.067E+03,1.150E+01, 3.412E+01, 1.922E+00, 0.000E+00/
      DATA (PH1(I,26,10, 2),I=1,6) /1.375E+03, 1.534E+02,9.191E+00, 2.691E+01, 4.849E+00, 0.000E+00/
      DATA (PH1(I,26,10, 3),I=1,6) /1.262E+03, 1.412E+02,3.368E+02, 5.569E+01, 4.328E+00, 1.177E-05/
      DATA (PH1(I,26,11, 1),I=1,6) /7.599E+03, 1.061E+03,1.161E+01, 3.713E+01, 1.889E+00, 0.000E+00/
      DATA (PH1(I,26,11, 2),I=1,6) /1.329E+03, 1.363E+02,9.940E+00, 3.026E+01, 4.885E+00, 0.000E+00/
      DATA (PH1(I,26,11, 3),I=1,6) /1.216E+03, 1.396E+02,3.383E+02, 5.459E+01, 4.366E+00, 1.179E-05/
      DATA (PH1(I,26,11, 4),I=1,6) /4.893E+02, 2.873E+01,1.207E+01, 5.150E+02, 3.846E+00, 0.000E+00/
      DATA (PH1(I,26,12, 1),I=1,6) /7.553E+03, 1.258E+03,8.097E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26,12, 2),I=1,6) /1.287E+03, 1.967E+02,7.555E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26,12, 3),I=1,6) /1.181E+03, 2.136E+02,1.496E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,26,12, 4),I=1,6) /4.570E+02, 4.064E+01,1.228E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,26,13, 1),I=1,6) /7.499E+03, 1.112E+03,1.052E+01, 3.562E+01, 1.871E+00, 0.000E+00/
      DATA (PH1(I,26,13, 2),I=1,6) /1.240E+03, 1.724E+02,8.240E+00, 2.693E+01, 4.675E+00, 0.000E+00/
      DATA (PH1(I,26,13, 3),I=1,6) /1.125E+03, 1.387E+02,3.367E+02, 5.279E+01, 4.407E+00, 1.111E-05/
      DATA (PH1(I,26,13, 4),I=1,6) /4.238E+02, 7.077E+01,7.969E+00, 5.650E+02, 3.329E+00, 0.000E+00/
      DATA (PH1(I,26,13, 5),I=1,6) /3.922E+02, 7.993E+01,1.135E+01, 5.807E+01, 5.676E+00, 6.512E-02/
      DATA (PH1(I,26,14, 1),I=1,6) /7.450E+03, 1.103E+03,1.069E+01, 3.585E+01, 1.875E+00, 0.000E+00/
      DATA (PH1(I,26,14, 2),I=1,6) /1.197E+03, 1.722E+02,8.199E+00, 2.626E+01, 4.710E+00, 0.000E+00/
      DATA (PH1(I,26,14, 3),I=1,6) /1.081E+03, 1.392E+02,3.316E+02, 5.237E+01, 4.410E+00, 1.107E-05/
      DATA (PH1(I,26,14, 4),I=1,6) /3.916E+02, 6.712E+01,8.093E+00, 5.886E+02, 3.360E+00, 0.000E+00/
      DATA (PH1(I,26,14, 5),I=1,6) /3.610E+02, 8.203E+01,2.092E+01, 6.128E+01, 5.596E+00, 3.441E-02/
      DATA (PH1(I,26,15, 1),I=1,6) /7.403E+03, 1.096E+03,1.083E+01, 3.612E+01, 1.878E+00, 0.000E+00/
      DATA (PH1(I,26,15, 2),I=1,6) /1.155E+03, 1.720E+02,8.161E+00, 2.570E+01, 4.741E+00, 0.000E+00/
      DATA (PH1(I,26,15, 3),I=1,6) /1.039E+03, 1.383E+02,3.325E+02, 5.090E+01, 4.446E+00, 1.114E-05/
      DATA (PH1(I,26,15, 4),I=1,6) /3.600E+02, 6.329E+01,8.143E+00, 5.637E+02, 3.421E+00, 0.000E+00/
      DATA (PH1(I,26,15, 5),I=1,6) /3.308E+02, 8.150E+01,3.022E+01, 6.039E+01, 5.612E+00, 2.772E-02/
      DATA (PH1(I,26,16, 1),I=1,6) /7.359E+03, 1.095E+03,1.084E+01, 3.598E+01, 1.880E+00, 0.000E+00/
      DATA (PH1(I,26,16, 2),I=1,6) /1.115E+03, 1.714E+02,8.138E+00, 2.516E+01, 4.774E+00, 0.000E+00/
      DATA (PH1(I,26,16, 3),I=1,6) /9.983E+02, 1.632E+02,2.449E+02, 5.452E+01, 4.187E+00, 1.594E-05/
      DATA (PH1(I,26,16, 4),I=1,6) /3.292E+02, 6.060E+01,7.963E+00, 5.214E+02, 3.487E+00, 0.000E+00/
      DATA (PH1(I,26,16, 5),I=1,6) /2.902E+02, 7.972E+01,3.922E+01, 5.574E+01, 5.729E+00, 2.828E-02/
      DATA (PH1(I,26,17, 1),I=1,6) /7.316E+03, 1.098E+03,1.078E+01, 3.596E+01, 1.878E+00, 0.000E+00/
      DATA (PH1(I,26,17, 2),I=1,6) /1.076E+03, 1.710E+02,8.135E+00, 2.440E+01, 4.817E+00, 0.000E+00/
      DATA (PH1(I,26,17, 3),I=1,6) /9.590E+02, 1.647E+02,2.392E+02, 5.276E+01, 4.204E+00, 1.574E-05/
      DATA (PH1(I,26,17, 4),I=1,6) /2.990E+02, 5.854E+01,7.650E+00, 4.679E+02, 3.561E+00, 0.000E+00/
      DATA (PH1(I,26,17, 5),I=1,6) /2.621E+02, 7.914E+01,4.686E+01, 5.125E+01, 5.835E+00, 2.574E-02/
      DATA (PH1(I,26,18, 1),I=1,6) /7.275E+03, 1.102E+03,1.068E+01, 3.576E+01, 1.878E+00, 0.000E+00/
      DATA (PH1(I,26,18, 2),I=1,6) /1.039E+03, 1.706E+02,8.130E+00, 2.387E+01, 4.850E+00, 0.000E+00/
      DATA (PH1(I,26,18, 3),I=1,6) /9.211E+02, 1.715E+02,2.205E+02, 5.298E+01, 4.154E+00, 1.508E-05/
      DATA (PH1(I,26,18, 4),I=1,6) /2.696E+02, 5.701E+01,7.237E+00, 4.092E+02, 3.643E+00, 0.000E+00/
      DATA (PH1(I,26,18, 5),I=1,6) /2.336E+02, 7.704E+01,5.530E+01, 4.380E+01, 6.059E+00, 2.673E-02/
      DATA (PH1(I,26,19, 1),I=1,6) /7.237E+03, 1.103E+03,1.067E+01, 3.584E+01, 1.876E+00, 0.000E+00/
      DATA (PH1(I,26,19, 2),I=1,6) /1.003E+03, 1.694E+02,8.158E+00, 2.345E+01, 4.887E+00, 0.000E+00/
      DATA (PH1(I,26,19, 3),I=1,6) /8.849E+02, 1.788E+02,2.032E+02, 5.293E+01, 4.107E+00, 1.398E-05/
      DATA (PH1(I,26,19, 4),I=1,6) /2.409E+02, 5.450E+01,6.938E+00, 3.605E+02, 3.741E+00, 0.000E+00/
      DATA (PH1(I,26,19, 5),I=1,6) /2.055E+02, 7.556E+01,5.336E+01, 3.855E+01, 6.265E+00, 2.667E-02/
      DATA (PH1(I,26,19, 6),I=1,6) /1.511E+02, 3.652E+01,1.323E+02, 9.277E+01, 6.150E+00, 2.051E-02/
      DATA (PH1(I,26,20, 1),I=1,6) /7.217E+03, 1.110E+03,1.055E+01, 3.563E+01, 1.873E+00, 0.000E+00/
      DATA (PH1(I,26,20, 2),I=1,6) /9.693E+02, 1.694E+02,8.138E+00, 2.316E+01, 4.905E+00, 0.000E+00/
      DATA (PH1(I,26,20, 3),I=1,6) /8.512E+02, 1.840E+02,1.920E+02, 5.238E+01, 4.082E+00, 1.377E-05/
      DATA (PH1(I,26,20, 4),I=1,6) /2.135E+02, 5.365E+01,6.463E+00, 3.169E+02, 3.824E+00, 0.000E+00/
      DATA (PH1(I,26,20, 5),I=1,6) /1.783E+02, 7.433E+01,5.210E+01, 3.306E+01, 6.509E+00, 2.404E-02/
      DATA (PH1(I,26,20, 6),I=1,6) /1.250E+02, 4.332E+01,1.701E+02, 9.791E+01, 5.905E+00, 2.604E-02/
      DATA (PH1(I,26,21, 1),I=1,6) /7.199E+03, 1.113E+03,1.049E+01, 3.537E+01, 1.874E+00, 0.000E+00/
      DATA (PH1(I,26,21, 2),I=1,6) /9.383E+02, 1.803E+02,7.682E+00, 2.415E+01, 4.757E+00, 0.000E+00/
      DATA (PH1(I,26,21, 3),I=1,6) /8.202E+02, 1.952E+02,1.709E+02, 5.690E+01, 3.953E+00, 1.330E-05/
      DATA (PH1(I,26,21, 4),I=1,6) /1.876E+02, 5.277E+01,6.047E+00, 2.797E+02, 3.910E+00, 0.000E+00/
      DATA (PH1(I,26,21, 5),I=1,6) /1.527E+02, 7.260E+01,5.246E+01, 2.751E+01, 6.823E+00, 2.105E-02/
      DATA (PH1(I,26,21, 6),I=1,6) /9.906E+01, 4.487E+01,2.082E+02, 9.579E+01, 5.914E+00, 2.368E-02/
      DATA (PH1(I,26,22, 1),I=1,6) /7.184E+03, 1.124E+03,1.027E+01, 3.534E+01, 1.866E+00, 0.000E+00/
      DATA (PH1(I,26,22, 2),I=1,6) /9.101E+02, 1.818E+02,7.623E+00, 2.401E+01, 4.752E+00, 0.000E+00/
      DATA (PH1(I,26,22, 3),I=1,6) /7.920E+02, 1.996E+02,1.636E+02, 5.719E+01, 3.924E+00, 1.320E-05/
      DATA (PH1(I,26,22, 4),I=1,6) /1.633E+02, 5.218E+01,5.672E+00, 2.482E+02, 3.991E+00, 0.000E+00/
      DATA (PH1(I,26,22, 5),I=1,6) /1.288E+02, 7.273E+01,5.143E+01, 2.428E+01, 7.028E+00, 1.365E-01/
      DATA (PH1(I,26,22, 6),I=1,6) /7.501E+01, 4.369E+01,2.661E+02, 6.273E+01, 6.344E+00, 2.402E-01/
      DATA (PH1(I,26,23, 1),I=1,6) /7.169E+03, 1.155E+03,9.711E+00, 3.329E+01, 1.869E+00, 0.000E+00/
      DATA (PH1(I,26,23, 2),I=1,6) /8.871E+02, 2.059E+02,6.695E+00, 2.772E+01, 4.410E+00, 0.000E+00/
      DATA (PH1(I,26,23, 3),I=1,6) /7.669E+02, 2.032E+02,1.578E+02, 5.933E+01, 3.879E+00, 1.308E-05/
      DATA (PH1(I,26,23, 4),I=1,6) /1.411E+02, 5.256E+01,5.295E+00, 2.332E+02, 4.036E+00, 0.000E+00/
      DATA (PH1(I,26,23, 5),I=1,6) /1.067E+02, 7.625E+01,4.810E+01, 2.286E+01, 7.043E+00, 2.449E-01/
      DATA (PH1(I,26,23, 6),I=1,6) /5.480E+01, 4.154E+01,3.683E+02, 3.529E+01, 7.056E+00, 3.316E-01/
      DATA (PH1(I,26,24, 1),I=1,6) /7.155E+03, 1.245E+03,8.309E+00, 3.170E+01, 1.826E+00, 0.000E+00/
      DATA (PH1(I,26,24, 2),I=1,6) /8.710E+02, 2.158E+02,6.454E+00, 2.722E+01, 4.355E+00, 0.000E+00/
      DATA (PH1(I,26,24, 3),I=1,6) /7.451E+02, 2.023E+02,1.591E+02, 5.990E+01, 3.878E+00, 1.309E-05/
      DATA (PH1(I,26,24, 4),I=1,6) /1.211E+02, 5.277E+01,5.015E+00, 2.131E+02, 4.093E+00, 0.000E+00/
      DATA (PH1(I,26,24, 5),I=1,6) /8.705E+01, 7.777E+01,4.422E+01, 2.336E+01, 7.017E+00, 2.557E-01/
      DATA (PH1(I,26,24, 6),I=1,6) /3.065E+01, 2.670E+01,6.301E+03, 5.385E+00, 1.232E+01, 4.407E-01/
      DATA (PH1(I,26,25, 1),I=1,6) /7.140E+03, 8.931E+02,1.666E+01, 2.381E+01, 2.262E+00, 0.000E+00/
      DATA (PH1(I,26,25, 2),I=1,6) /8.608E+02, 1.431E+02,9.289E+00, 2.026E+01, 5.376E+00, 0.000E+00/
      DATA (PH1(I,26,25, 3),I=1,6) /7.341E+02, 2.281E+02,1.241E+02, 8.058E+01, 3.572E+00, 1.223E-01/
      DATA (PH1(I,26,25, 4),I=1,6) /1.102E+02, 4.663E+01,5.434E+00, 9.271E+01, 4.640E+00, 0.000E+00/
      DATA (PH1(I,26,25, 5),I=1,6) /7.617E+01, 7.750E+01,4.624E+01, 2.155E+01, 7.138E+00, 2.599E-01/
      DATA (PH1(I,26,25, 6),I=1,6) /2.193E+01, 1.933E+01,3.679E+04, 3.564E+00, 1.566E+01, 4.144E-01/
      DATA (PH1(I,26,25, 7),I=1,6) /1.619E+01, 1.014E+01,1.084E+00, 2.562E+04, 4.167E+00, 1.598E-02/
      DATA (PH1(I,26,26, 1),I=1,6) /7.124E+03, 8.044E+02,2.055E+01, 3.633E+01, 2.118E+00, 0.000E+00/
      DATA (PH1(I,26,26, 2),I=1,6) /8.570E+02, 5.727E+01,1.076E+01, 2.785E+01, 6.635E+00, 0.000E+00/
      DATA (PH1(I,26,26, 3),I=1,6) /7.240E+02, 2.948E+02,7.191E+01, 3.219E+02, 2.837E+00, 6.314E-02/
      DATA (PH1(I,26,26, 4),I=1,6) /1.040E+02, 4.334E+01,5.921E+00, 5.293E+01, 5.129E+00, 0.000E+00/
      DATA (PH1(I,26,26, 5),I=1,6) /6.600E+01, 7.630E+01,6.298E+01, 1.479E+01, 7.672E+00, 2.646E-01/
      DATA (PH1(I,26,26, 6),I=1,6) /1.470E+01, 1.407E+01,1.850E+04, 4.458E+00, 1.691E+01, 4.039E-01/
      DATA (PH1(I,26,26, 7),I=1,6) /7.902E+00, 1.277E+01,1.468E+00, 1.116E+05, 4.112E+00, 3.238E-02/
      DATA (PH1(I,27, 1, 1),I=1,6) /1.001E+04, 3.163E+02,7.510E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,27, 2, 1),I=1,6) /9.545E+03, 1.251E+03,9.046E+00, 6.786E+01, 1.613E+00, 0.000E+00/
      DATA (PH1(I,27, 3, 1),I=1,6) /9.347E+03, 1.249E+03,9.384E+00, 6.378E+01, 1.609E+00, 0.000E+00/
      DATA (PH1(I,27, 3, 2),I=1,6) /2.219E+03, 2.352E+02,4.599E+00, 3.873E+01, 3.808E+00, 0.000E+00/
      DATA (PH1(I,27, 4, 1),I=1,6) /9.205E+03, 1.267E+03,9.026E+00, 6.059E+01, 1.619E+00, 0.000E+00/
      DATA (PH1(I,27, 4, 2),I=1,6) /2.119E+03, 1.593E+02,1.253E+01, 3.366E+01, 4.493E+00, 0.000E+00/
      DATA (PH1(I,27, 5, 1),I=1,6) /9.046E+03, 1.275E+03,8.806E+00, 5.513E+01, 1.651E+00, 0.000E+00/
      DATA (PH1(I,27, 5, 2),I=1,6) /2.012E+03, 1.799E+02,1.074E+01, 2.878E+01, 4.518E+00, 0.000E+00/
      DATA (PH1(I,27, 5, 3),I=1,6) /1.961E+03, 1.923E+02,4.987E+01, 9.392E+01, 3.631E+00, 2.239E-05/
      DATA (PH1(I,27, 6, 1),I=1,6) /8.890E+03, 1.277E+03,8.703E+00, 5.267E+01, 1.668E+00, 0.000E+00/
      DATA (PH1(I,27, 6, 2),I=1,6) /1.910E+03, 1.760E+02,1.043E+01, 2.775E+01, 4.612E+00, 0.000E+00/
      DATA (PH1(I,27, 6, 3),I=1,6) /1.846E+03, 2.159E+02,7.348E+01, 9.457E+01, 3.535E+00, 1.182E-05/
      DATA (PH1(I,27, 7, 1),I=1,6) /8.718E+03, 1.330E+03,8.042E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27, 7, 2),I=1,6) /1.803E+03, 2.030E+02,8.771E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27, 7, 3),I=1,6) /1.735E+03, 2.134E+02,1.027E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,27, 8, 1),I=1,6) /8.595E+03, 1.260E+03,8.868E+00, 4.790E+01, 1.717E+00, 0.000E+00/
      DATA (PH1(I,27, 8, 2),I=1,6) /1.711E+03, 1.706E+02,9.597E+00, 2.731E+01, 4.722E+00, 0.000E+00/
      DATA (PH1(I,27, 8, 3),I=1,6) /1.603E+03, 2.448E+02,9.893E+01, 9.263E+01, 3.432E+00, 9.870E-06/
      DATA (PH1(I,27, 9, 1),I=1,6) /8.433E+03, 1.344E+03,7.779E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27, 9, 2),I=1,6) /1.606E+03, 2.084E+02,7.728E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27, 9, 3),I=1,6) /1.505E+03, 2.157E+02,1.481E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,27,10, 1),I=1,6) /8.315E+03, 1.252E+03,8.903E+00, 4.635E+01, 1.738E+00, 0.000E+00/
      DATA (PH1(I,27,10, 2),I=1,6) /1.519E+03, 1.643E+02,8.754E+00, 2.695E+01, 4.851E+00, 0.000E+00/
      DATA (PH1(I,27,10, 3),I=1,6) /1.397E+03, 2.062E+02,1.816E+02, 6.860E+01, 3.824E+00, 9.070E-06/
      DATA (PH1(I,27,11, 1),I=1,6) /8.260E+03, 1.249E+03,8.967E+00, 4.909E+01, 1.719E+00, 0.000E+00/
      DATA (PH1(I,27,11, 2),I=1,6) /1.470E+03, 1.433E+02,9.597E+00, 2.996E+01, 4.928E+00, 0.000E+00/
      DATA (PH1(I,27,11, 3),I=1,6) /1.362E+03, 2.054E+02,1.804E+02, 6.660E+01, 3.855E+00, 8.999E-06/
      DATA (PH1(I,27,11, 4),I=1,6) /5.466E+02, 1.262E+01,2.734E+01, 6.463E+02, 4.280E+00, 0.000E+00/
      DATA (PH1(I,27,12, 1),I=1,6) /8.207E+03, 1.277E+03,8.536E+00, 4.760E+01, 1.714E+00, 0.000E+00/
      DATA (PH1(I,27,12, 2),I=1,6) /1.423E+03, 1.673E+02,8.543E+00, 2.705E+01, 4.820E+00, 0.000E+00/
      DATA (PH1(I,27,12, 3),I=1,6) /1.314E+03, 2.058E+02,1.783E+02, 6.530E+01, 3.868E+00, 8.972E-06/
      DATA (PH1(I,27,12, 4),I=1,6) /5.120E+02, 2.690E+01,2.122E+01, 3.166E+02, 4.210E+00, 0.000E+00/
      DATA (PH1(I,27,13, 1),I=1,6) /8.154E+03, 1.277E+03,8.527E+00, 4.599E+01, 1.727E+00, 0.000E+00/
      DATA (PH1(I,27,13, 2),I=1,6) /1.376E+03, 1.726E+02,8.292E+00, 2.627E+01, 4.812E+00, 0.000E+00/
      DATA (PH1(I,27,13, 3),I=1,6) /1.266E+03, 1.933E+02,2.002E+02, 6.453E+01, 3.945E+00, 9.099E-06/
      DATA (PH1(I,27,13, 4),I=1,6) /4.777E+02, 2.652E+01,1.977E+01, 2.830E+02, 4.297E+00, 0.000E+00/
      DATA (PH1(I,27,13, 5),I=1,6) /4.440E+02, 1.358E+02,6.171E+00, 4.715E+01, 5.247E+00, 8.240E-01/
      DATA (PH1(I,27,14, 1),I=1,6) /8.102E+03, 1.272E+03,8.593E+00, 4.604E+01, 1.730E+00, 0.000E+00/
      DATA (PH1(I,27,14, 2),I=1,6) /1.331E+03, 1.718E+02,8.247E+00, 2.592E+01, 4.840E+00, 0.000E+00/
      DATA (PH1(I,27,14, 3),I=1,6) /1.219E+03, 1.931E+02,1.993E+02, 6.407E+01, 3.953E+00, 9.111E-06/
      DATA (PH1(I,27,14, 4),I=1,6) /4.436E+02, 2.529E+01,1.979E+01, 2.943E+02, 4.315E+00, 0.000E+00/
      DATA (PH1(I,27,14, 5),I=1,6) /4.110E+02, 9.918E+01,1.786E+01, 5.216E+01, 5.592E+00, 5.707E-01/
      DATA (PH1(I,27,15, 1),I=1,6) /8.045E+03, 1.361E+03,7.456E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,15, 2),I=1,6) /1.281E+03, 2.135E+02,6.900E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,15, 3),I=1,6) /1.167E+03, 2.301E+02,1.403E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,15, 4),I=1,6) /4.088E+02, 4.412E+01,1.052E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,15, 5),I=1,6) /3.790E+02, 8.885E+01,2.993E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,16, 1),I=1,6) /8.005E+03, 1.264E+03,8.703E+00, 4.638E+01, 1.732E+00, 0.000E+00/
      DATA (PH1(I,27,16, 2),I=1,6) /1.244E+03, 1.755E+02,8.017E+00, 2.508E+01, 4.854E+00, 0.000E+00/
      DATA (PH1(I,27,16, 3),I=1,6) /1.131E+03, 1.925E+02,1.972E+02, 6.228E+01, 3.980E+00, 9.096E-06/
      DATA (PH1(I,27,16, 4),I=1,6) /3.775E+02, 2.185E+01,1.954E+01, 3.125E+02, 4.423E+00, 0.000E+00/
      DATA (PH1(I,27,16, 5),I=1,6) /3.360E+02, 9.235E+01,3.595E+01, 4.513E+01, 5.839E+00, 5.704E-01/
      DATA (PH1(I,27,17, 1),I=1,6) /7.951E+03, 1.362E+03,7.447E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,17, 2),I=1,6) /1.196E+03, 2.146E+02,6.772E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,17, 3),I=1,6) /1.080E+03, 2.305E+02,1.383E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,17, 4),I=1,6) /3.439E+02, 4.522E+01,9.201E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,17, 5),I=1,6) /3.053E+02, 8.909E+01,4.487E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,18, 1),I=1,6) /7.915E+03, 1.266E+03,8.671E+00, 4.646E+01, 1.731E+00, 0.000E+00/
      DATA (PH1(I,27,18, 2),I=1,6) /1.163E+03, 1.808E+02,7.771E+00, 2.435E+01, 4.849E+00, 0.000E+00/
      DATA (PH1(I,27,18, 3),I=1,6) /1.048E+03, 1.922E+02,1.949E+02, 5.995E+01, 4.013E+00, 9.089E-06/
      DATA (PH1(I,27,18, 4),I=1,6) /3.142E+02, 3.099E+01,1.193E+01, 2.258E+02, 4.398E+00, 0.000E+00/
      DATA (PH1(I,27,18, 5),I=1,6) /2.754E+02, 9.052E+01,5.045E+01, 4.052E+01, 6.001E+00, 4.620E-01/
      DATA (PH1(I,27,19, 1),I=1,6) /7.877E+03, 1.287E+03,8.375E+00, 4.957E+01, 1.697E+00, 0.000E+00/
      DATA (PH1(I,27,19, 2),I=1,6) /1.123E+03, 1.802E+02,7.760E+00, 2.395E+01, 4.878E+00, 0.000E+00/
      DATA (PH1(I,27,19, 3),I=1,6) /1.009E+03, 1.943E+02,1.901E+02, 5.795E+01, 4.029E+00, 8.959E-06/
      DATA (PH1(I,27,19, 4),I=1,6) /2.833E+02, 3.765E+01,9.284E+00, 1.936E+02, 4.355E+00, 0.000E+00/
      DATA (PH1(I,27,19, 5),I=1,6) /2.450E+02, 9.011E+01,4.789E+01, 3.807E+01, 6.102E+00, 4.132E-01/
      DATA (PH1(I,27,19, 6),I=1,6) /1.861E+02, 3.461E+01,1.735E+02, 7.681E+01, 6.518E+00, 1.706E-04/
      DATA (PH1(I,27,20, 1),I=1,6) /7.840E+03, 1.363E+03,7.433E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,20, 2),I=1,6) /1.086E+03, 2.158E+02,6.639E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,20, 3),I=1,6) /9.692E+02, 2.313E+02,1.357E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,20, 4),I=1,6) /2.544E+02, 4.790E+01,7.244E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,20, 5),I=1,6) /2.158E+02, 9.095E+01,4.368E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,20, 6),I=1,6) /1.578E+02, 4.257E+01,2.299E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,27,21, 1),I=1,6) /7.820E+03, 1.295E+03,8.260E+00, 5.104E+01, 1.682E+00, 0.000E+00/
      DATA (PH1(I,27,21, 2),I=1,6) /1.052E+03, 1.797E+02,7.740E+00, 2.344E+01, 4.914E+00, 0.000E+00/
      DATA (PH1(I,27,21, 3),I=1,6) /9.379E+02, 2.015E+02,1.767E+02, 5.582E+01, 4.018E+00, 8.773E-06/
      DATA (PH1(I,27,21, 4),I=1,6) /2.255E+02, 4.430E+01,6.990E+00, 1.659E+02, 4.350E+00, 0.000E+00/
      DATA (PH1(I,27,21, 5),I=1,6) /1.876E+02, 8.695E+01,4.460E+01, 3.324E+01, 6.378E+00, 3.125E-01/
      DATA (PH1(I,27,21, 6),I=1,6) /1.290E+02, 4.433E+01,2.748E+02, 5.829E+01, 6.482E+00, 1.086E-04/
      DATA (PH1(I,27,22, 1),I=1,6) /7.803E+03, 1.364E+03,7.423E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,22, 2),I=1,6) /1.025E+03, 2.163E+02,6.587E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,22, 3),I=1,6) /9.074E+02, 2.319E+02,1.344E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,22, 4),I=1,6) /1.999E+02, 5.038E+01,6.114E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,22, 5),I=1,6) /1.619E+02, 9.318E+01,3.733E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,22, 6),I=1,6) /1.020E+02, 4.352E+01,3.483E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,27,23, 1),I=1,6) /7.788E+03, 1.313E+03,8.026E+00, 5.052E+01, 1.675E+00, 0.000E+00/
      DATA (PH1(I,27,23, 2),I=1,6) /9.918E+02, 1.840E+02,7.606E+00, 2.298E+01, 4.903E+00, 0.000E+00/
      DATA (PH1(I,27,23, 3),I=1,6) /8.775E+02, 2.077E+02,1.666E+02, 5.503E+01, 3.993E+00, 8.505E-06/
      DATA (PH1(I,27,23, 4),I=1,6) /1.740E+02, 5.017E+01,5.619E+00, 1.481E+02, 4.342E+00, 0.000E+00/
      DATA (PH1(I,27,23, 5),I=1,6) /1.366E+02, 8.377E+01,4.275E+01, 2.773E+01, 6.744E+00, 2.457E-01/
      DATA (PH1(I,27,23, 6),I=1,6) /7.950E+01, 4.369E+01,3.792E+02, 4.396E+01, 6.895E+00, 1.079E-04/
      DATA (PH1(I,27,24, 1),I=1,6) /7.773E+03, 1.365E+03,7.414E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,24, 2),I=1,6) /9.708E+02, 2.166E+02,6.563E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,24, 3),I=1,6) /8.556E+02, 2.326E+02,1.332E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,24, 4),I=1,6) /1.495E+02, 5.342E+01,5.159E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,24, 5),I=1,6) /1.128E+02, 9.621E+01,3.175E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,24, 6),I=1,6) /5.127E+01, 4.530E+01,3.755E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,27,25, 1),I=1,6) /7.758E+03, 1.487E+03,6.149E+00, 5.282E+01, 1.571E+00, 0.000E+00/
      DATA (PH1(I,27,25, 2),I=1,6) /9.547E+02, 2.172E+02,6.667E+00, 2.353E+01, 4.602E+00, 0.000E+00/
      DATA (PH1(I,27,25, 3),I=1,6) /8.295E+02, 2.089E+02,1.651E+02, 5.457E+01, 3.992E+00, 8.474E-06/
      DATA (PH1(I,27,25, 4),I=1,6) /1.305E+02, 5.435E+01,4.874E+00, 1.379E+02, 4.336E+00, 0.000E+00/
      DATA (PH1(I,27,25, 5),I=1,6) /9.362E+01, 8.419E+01,4.192E+01, 2.297E+01, 7.043E+00, 2.547E-01/
      DATA (PH1(I,27,25, 6),I=1,6) /3.350E+01, 2.181E+01,1.205E+03, 1.434E+01, 1.114E+01, 1.266E-02/
      DATA (PH1(I,27,26, 1),I=1,6) /7.742E+03, 2.079E+03,3.020E+00, 1.067E+03, 9.653E-01, 0.000E+00/
      DATA (PH1(I,27,26, 2),I=1,6) /9.443E+02, 1.709E+02,8.001E+00, 2.274E+01, 5.044E+00, 0.000E+00/
      DATA (PH1(I,27,26, 3),I=1,6) /8.130E+02, 1.989E+02,1.810E+02, 5.274E+01, 4.078E+00, 1.431E-04/
      DATA (PH1(I,27,26, 4),I=1,6) /1.213E+02, 5.602E+01,4.666E+00, 1.358E+02, 4.323E+00, 0.000E+00/
      DATA (PH1(I,27,26, 5),I=1,6) /7.621E+01, 8.408E+01,4.300E+01, 2.074E+01, 7.214E+00, 2.553E-01/
      DATA (PH1(I,27,26, 6),I=1,6) /1.708E+01, 2.419E+00,1.275E+01, 2.823E+01, 1.843E+01, 4.806E-05/
      DATA (PH1(I,27,27, 1),I=1,6) /7.725E+03, 6.269E+02,3.582E+01, 3.161E+01, 2.476E+00, 0.000E+00/
      DATA (PH1(I,27,27, 2),I=1,6) /9.400E+02, 8.888E+01,1.030E+01, 2.797E+01, 5.913E+00, 0.000E+00/
      DATA (PH1(I,27,27, 3),I=1,6) /8.000E+02, 2.832E+02,9.075E+01, 7.686E+01, 3.416E+00, 4.833E-05/
      DATA (PH1(I,27,27, 4),I=1,6) /1.150E+02, 5.097E+01,4.896E+00, 1.198E+02, 4.513E+00, 0.000E+00/
      DATA (PH1(I,27,27, 5),I=1,6) /7.300E+01, 8.256E+01,4.587E+01, 1.973E+01, 7.331E+00, 2.573E-01/
      DATA (PH1(I,27,27, 6),I=1,6) /1.580E+01, 1.581E+01,4.931E+03, 6.607E+00, 1.532E+01, 3.676E-01/
      DATA (PH1(I,27,27, 7),I=1,6) /7.864E+00, 1.370E+01,1.555E+00, 7.559E+02, 4.337E+00, 3.355E-02/
      DATA (PH1(I,28, 1, 1),I=1,6) /1.078E+04, 3.406E+02,6.983E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,28, 2, 1),I=1,6) /1.029E+04, 1.349E+03,8.382E+00, 7.431E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,28, 3, 1),I=1,6) /1.008E+04, 1.373E+03,8.336E+00, 7.347E+01, 1.558E+00, 0.000E+00/
      DATA (PH1(I,28, 3, 2),I=1,6) /2.399E+03, 1.314E+02,8.142E+00, 4.522E+01, 4.477E+00, 0.000E+00/
      DATA (PH1(I,28, 4, 1),I=1,6) /9.937E+03, 1.382E+03,8.149E+00, 6.782E+01, 1.580E+00, 0.000E+00/
      DATA (PH1(I,28, 4, 2),I=1,6) /2.295E+03, 1.813E+02,1.117E+01, 3.332E+01, 4.426E+00, 0.000E+00/
      DATA (PH1(I,28, 5, 1),I=1,6) /9.771E+03, 1.404E+03,7.788E+00, 6.396E+01, 1.592E+00, 0.000E+00/
      DATA (PH1(I,28, 5, 2),I=1,6) /2.184E+03, 1.649E+02,1.145E+01, 2.997E+01, 4.706E+00, 0.000E+00/
      DATA (PH1(I,28, 5, 3),I=1,6) /2.131E+03, 2.473E+02,3.198E+01, 1.014E+02, 3.435E+00, 2.062E-06/
      DATA (PH1(I,28, 6, 1),I=1,6) /9.609E+03, 1.415E+03,7.600E+00, 6.229E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28, 6, 2),I=1,6) /2.077E+03, 1.726E+02,1.052E+01, 2.738E+01, 4.768E+00, 0.000E+00/
      DATA (PH1(I,28, 6, 3),I=1,6) /2.011E+03, 2.139E+02,8.219E+01, 9.511E+01, 3.615E+00, 6.838E-07/
      DATA (PH1(I,28, 7, 1),I=1,6) /9.423E+03, 1.429E+03,7.520E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28, 7, 2),I=1,6) /1.965E+03, 2.177E+02,8.285E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28, 7, 3),I=1,6) /1.894E+03, 2.300E+02,9.666E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,28, 8, 1),I=1,6) /9.300E+03, 1.401E+03,7.673E+00, 5.681E+01, 1.640E+00, 0.000E+00/
      DATA (PH1(I,28, 8, 2),I=1,6) /1.870E+03, 1.608E+02,9.953E+00, 2.835E+01, 4.878E+00, 0.000E+00/
      DATA (PH1(I,28, 8, 3),I=1,6) /1.756E+03, 2.140E+02,1.434E+02, 7.628E+01, 3.752E+00, 5.528E-07/
      DATA (PH1(I,28, 9, 1),I=1,6) /9.125E+03, 1.443E+03,7.283E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28, 9, 2),I=1,6) /1.760E+03, 2.233E+02,7.335E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28, 9, 3),I=1,6) /1.648E+03, 2.321E+02,1.405E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,28,10, 1),I=1,6) /9.007E+03, 1.400E+03,7.626E+00, 5.528E+01, 1.653E+00, 0.000E+00/
      DATA (PH1(I,28,10, 2),I=1,6) /1.669E+03, 1.441E+02,9.448E+00, 2.779E+01, 5.132E+00, 0.000E+00/
      DATA (PH1(I,28,10, 3),I=1,6) /1.541E+03, 2.358E+02,1.530E+02, 7.451E+01, 3.705E+00, 4.687E-07/
      DATA (PH1(I,28,11, 1),I=1,6) /8.950E+03, 1.383E+03,7.837E+00, 5.839E+01, 1.644E+00, 0.000E+00/
      DATA (PH1(I,28,11, 2),I=1,6) /1.618E+03, 1.444E+02,9.470E+00, 3.104E+01, 4.981E+00, 0.000E+00/
      DATA (PH1(I,28,11, 3),I=1,6) /1.506E+03, 2.334E+02,1.540E+02, 7.184E+01, 3.746E+00, 4.400E-07/
      DATA (PH1(I,28,11, 4),I=1,6) /6.071E+02, 1.032E+01,3.674E+01, 8.252E+02, 4.304E+00, 0.000E+00/
      DATA (PH1(I,28,12, 1),I=1,6) /8.894E+03, 1.442E+03,7.159E+00, 6.122E+01, 1.602E+00, 0.000E+00/
      DATA (PH1(I,28,12, 2),I=1,6) /1.569E+03, 1.824E+02,8.013E+00, 2.738E+01, 4.780E+00, 0.000E+00/
      DATA (PH1(I,28,12, 3),I=1,6) /1.456E+03, 2.327E+02,1.539E+02, 7.045E+01, 3.763E+00, 4.221E-07/
      DATA (PH1(I,28,12, 4),I=1,6) /5.713E+02, 2.577E+01,2.339E+01, 3.455E+02, 4.232E+00, 0.000E+00/
      DATA (PH1(I,28,13, 1),I=1,6) /8.838E+03, 1.452E+03,7.048E+00, 6.036E+01, 1.602E+00, 0.000E+00/
      DATA (PH1(I,28,13, 2),I=1,6) /1.519E+03, 1.853E+02,7.851E+00, 2.679E+01, 4.785E+00, 0.000E+00/
      DATA (PH1(I,28,13, 3),I=1,6) /1.405E+03, 2.180E+02,1.742E+02, 6.951E+01, 3.842E+00, 4.401E-07/
      DATA (PH1(I,28,13, 4),I=1,6) /5.347E+02, 2.505E+01,2.361E+01, 3.622E+02, 4.225E+00, 0.000E+00/
      DATA (PH1(I,28,13, 5),I=1,6) /4.984E+02, 1.711E+02,4.708E+00, 4.790E+01, 5.019E+00, 8.478E-01/
      DATA (PH1(I,28,14, 1),I=1,6) /8.783E+03, 1.448E+03,7.081E+00, 6.137E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,14, 2),I=1,6) /1.471E+03, 1.841E+02,7.825E+00, 2.644E+01, 4.816E+00, 0.000E+00/
      DATA (PH1(I,28,14, 3),I=1,6) /1.356E+03, 2.108E+02,1.846E+02, 6.874E+01, 3.886E+00, 4.288E-07/
      DATA (PH1(I,28,14, 4),I=1,6) /4.988E+02, 2.375E+01,2.383E+01, 3.786E+02, 4.245E+00, 0.000E+00/
      DATA (PH1(I,28,14, 5),I=1,6) /4.637E+02, 1.543E+02,1.053E+01, 4.743E+01, 5.176E+00, 8.373E-01/
      DATA (PH1(I,28,15, 1),I=1,6) /8.720E+03, 1.467E+03,6.912E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,15, 2),I=1,6) /1.419E+03, 2.287E+02,6.549E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,15, 3),I=1,6) /1.299E+03, 2.464E+02,1.356E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,15, 4),I=1,6) /4.620E+02, 4.702E+01,1.037E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,15, 5),I=1,6) /4.302E+02, 9.587E+01,2.890E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,16, 1),I=1,6) /8.680E+03, 1.439E+03,7.182E+00, 6.248E+01, 1.598E+00, 0.000E+00/
      DATA (PH1(I,28,16, 2),I=1,6) /1.379E+03, 1.886E+02,7.589E+00, 2.560E+01, 4.825E+00, 0.000E+00/
      DATA (PH1(I,28,16, 3),I=1,6) /1.262E+03, 2.088E+02,1.852E+02, 6.640E+01, 3.924E+00, 4.306E-07/
      DATA (PH1(I,28,16, 4),I=1,6) /4.291E+02, 2.003E+01,2.274E+01, 3.611E+02, 4.438E+00, 0.000E+00/
      DATA (PH1(I,28,16, 5),I=1,6) /3.840E+02, 1.158E+02,2.906E+01, 4.380E+01, 5.646E+00, 7.430E-01/
      DATA (PH1(I,28,17, 1),I=1,6) /8.620E+03, 1.467E+03,6.916E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,17, 2),I=1,6) /1.328E+03, 2.299E+02,6.430E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,17, 3),I=1,6) /1.207E+03, 2.467E+02,1.337E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,17, 4),I=1,6) /3.933E+02, 4.806E+01,9.150E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,17, 5),I=1,6) /3.521E+02, 9.596E+01,4.382E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,18, 1),I=1,6) /8.584E+03, 1.440E+03,7.163E+00, 6.246E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,18, 2),I=1,6) /1.293E+03, 1.925E+02,7.399E+00, 2.482E+01, 4.836E+00, 0.000E+00/
      DATA (PH1(I,28,18, 3),I=1,6) /1.174E+03, 2.087E+02,1.825E+02, 6.397E+01, 3.955E+00, 4.274E-07/
      DATA (PH1(I,28,18, 4),I=1,6) /3.620E+02, 2.901E+01,1.340E+01, 2.515E+02, 4.429E+00, 0.000E+00/
      DATA (PH1(I,28,18, 5),I=1,6) /3.210E+02, 1.005E+02,4.742E+01, 4.168E+01, 5.919E+00, 5.228E-01/
      DATA (PH1(I,28,19, 1),I=1,6) /8.542E+03, 1.439E+03,7.170E+00, 6.247E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,19, 2),I=1,6) /1.251E+03, 1.931E+02,7.348E+00, 2.450E+01, 4.851E+00, 0.000E+00/
      DATA (PH1(I,28,19, 3),I=1,6) /1.132E+03, 2.109E+02,1.780E+02, 6.201E+01, 3.969E+00, 4.273E-07/
      DATA (PH1(I,28,19, 4),I=1,6) /3.290E+02, 3.448E+01,1.054E+01, 2.150E+02, 4.412E+00, 0.000E+00/
      DATA (PH1(I,28,19, 5),I=1,6) /2.877E+02, 1.033E+02,4.372E+01, 3.878E+01, 5.977E+00, 5.276E-01/
      DATA (PH1(I,28,19, 6),I=1,6) /2.246E+02, 3.839E+01,1.714E+02, 7.972E+01, 6.440E+00, 1.220E-04/
      DATA (PH1(I,28,20, 1),I=1,6) /8.503E+03, 1.439E+03,7.178E+00, 6.253E+01, 1.598E+00, 0.000E+00/
      DATA (PH1(I,28,20, 2),I=1,6) /1.211E+03, 1.935E+02,7.316E+00, 2.416E+01, 4.868E+00, 0.000E+00/
      DATA (PH1(I,28,20, 3),I=1,6) /1.092E+03, 2.189E+02,1.653E+02, 6.160E+01, 3.932E+00, 4.163E-07/
      DATA (PH1(I,28,20, 4),I=1,6) /2.972E+02, 3.985E+01,8.589E+00, 1.892E+02, 4.394E+00, 0.000E+00/
      DATA (PH1(I,28,20, 5),I=1,6) /2.561E+02, 1.023E+02,4.197E+01, 3.659E+01, 6.079E+00, 4.805E-01/
      DATA (PH1(I,28,20, 6),I=1,6) /1.930E+02, 4.212E+01,2.596E+02, 6.923E+01, 6.485E+00, 1.416E-04/
      DATA (PH1(I,28,21, 1),I=1,6) /8.472E+03, 1.466E+03,6.924E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,21, 2),I=1,6) /1.174E+03, 2.314E+02,6.277E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,21, 3),I=1,6) /1.051E+03, 2.476E+02,1.308E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,21, 4),I=1,6) /2.682E+02, 5.179E+01,6.728E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,21, 5),I=1,6) /2.266E+02, 9.852E+01,4.031E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,21, 6),I=1,6) /1.620E+02, 4.735E+01,2.988E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,28,22, 1),I=1,6) /8.452E+03, 1.430E+03,7.274E+00, 6.333E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,22, 2),I=1,6) /1.139E+03, 1.947E+02,7.253E+00, 2.382E+01, 4.880E+00, 0.000E+00/
      DATA (PH1(I,28,22, 3),I=1,6) /1.019E+03, 2.193E+02,1.639E+02, 5.860E+01, 3.969E+00, 4.138E-07/
      DATA (PH1(I,28,22, 4),I=1,6) /2.378E+02, 4.756E+01,6.461E+00, 1.609E+02, 4.380E+00, 0.000E+00/
      DATA (PH1(I,28,22, 5),I=1,6) /1.971E+02, 9.385E+01,4.144E+01, 3.216E+01, 6.437E+00, 3.229E-01/
      DATA (PH1(I,28,22, 6),I=1,6) /1.330E+02, 5.059E+01,3.124E+02, 5.417E+01, 6.505E+00, 1.508E-04/
      DATA (PH1(I,28,23, 1),I=1,6) /8.434E+03, 1.466E+03,6.929E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,23, 2),I=1,6) /1.112E+03, 2.318E+02,6.240E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,23, 3),I=1,6) /9.886E+02, 2.482E+02,1.297E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,23, 4),I=1,6) /2.119E+02, 5.449E+01,5.709E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,23, 5),I=1,6) /1.709E+02, 1.010E+02,3.462E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,23, 6),I=1,6) /1.080E+02, 4.837E+01,3.797E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,28,24, 1),I=1,6) /8.418E+03, 1.463E+03,6.925E+00, 6.343E+01, 1.583E+00, 0.000E+00/
      DATA (PH1(I,28,24, 2),I=1,6) /1.077E+03, 1.956E+02,7.233E+00, 2.315E+01, 4.910E+00, 0.000E+00/
      DATA (PH1(I,28,24, 3),I=1,6) /9.576E+02, 2.215E+02,1.608E+02, 5.699E+01, 3.980E+00, 4.065E-07/
      DATA (PH1(I,28,24, 4),I=1,6) /1.849E+02, 5.389E+01,5.242E+00, 1.434E+02, 4.370E+00, 0.000E+00/
      DATA (PH1(I,28,24, 5),I=1,6) /1.447E+02, 9.036E+01,4.016E+01, 2.699E+01, 6.796E+00, 2.509E-01/
      DATA (PH1(I,28,24, 6),I=1,6) /7.610E+01, 4.694E+01,4.360E+02, 3.874E+01, 7.119E+00, 3.263E-03/
      DATA (PH1(I,28,25, 1),I=1,6) /8.402E+03, 1.466E+03,6.933E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,25, 2),I=1,6) /1.059E+03, 2.320E+02,6.226E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,25, 3),I=1,6) /9.358E+02, 2.488E+02,1.288E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,25, 4),I=1,6) /1.597E+02, 5.774E+01,4.849E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,25, 5),I=1,6) /1.201E+02, 1.043E+02,2.963E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,25, 6),I=1,6) /5.490E+01, 5.023E+01,3.918E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,28,26, 1),I=1,6) /8.386E+03, 1.647E+03,5.367E+00, 6.758E+01, 1.485E+00, 0.000E+00/
      DATA (PH1(I,28,26, 2),I=1,6) /1.040E+03, 2.324E+02,6.314E+00, 2.361E+01, 4.604E+00, 0.000E+00/
      DATA (PH1(I,28,26, 3),I=1,6) /9.084E+02, 2.230E+02,1.588E+02, 5.650E+01, 3.977E+00, 4.038E-07/
      DATA (PH1(I,28,26, 4),I=1,6) /1.401E+02, 5.823E+01,4.586E+00, 1.330E+02, 4.367E+00, 0.000E+00/
      DATA (PH1(I,28,26, 5),I=1,6) /1.003E+02, 9.061E+01,3.980E+01, 2.247E+01, 7.093E+00, 2.543E-01/
      DATA (PH1(I,28,26, 6),I=1,6) /3.532E+01, 2.971E+01,1.965E+03, 9.751E+00, 1.140E+01, 4.019E-01/
      DATA (PH1(I,28,27, 1),I=1,6) /8.368E+03, 2.249E+03,2.776E+00, 1.325E+03, 9.540E-01, 0.000E+00/
      DATA (PH1(I,28,27, 2),I=1,6) /1.029E+03, 1.808E+02,7.632E+00, 2.289E+01, 5.060E+00, 0.000E+00/
      DATA (PH1(I,28,27, 3),I=1,6) /8.903E+02, 2.205E+02,1.625E+02, 5.653E+01, 3.990E+00, 1.246E-04/
      DATA (PH1(I,28,27, 4),I=1,6) /1.315E+02, 6.024E+01,4.394E+00, 1.314E+02, 4.347E+00, 0.000E+00/
      DATA (PH1(I,28,27, 5),I=1,6) /8.232E+01, 9.066E+01,4.075E+01, 2.048E+01, 7.246E+00, 2.545E-01/
      DATA (PH1(I,28,27, 6),I=1,6) /1.817E+01, 1.501E+00,9.258E-01, 4.442E+01, 1.928E+01, 3.708E-02/
      DATA (PH1(I,28,28, 1),I=1,6) /8.348E+03, 7.366E+02,2.836E+01, 3.622E+01, 2.316E+00, 0.000E+00/
      DATA (PH1(I,28,28, 2),I=1,6) /1.024E+03, 1.132E+02,9.424E+00, 2.712E+01, 5.643E+00, 0.000E+00/
      DATA (PH1(I,28,28, 3),I=1,6) /8.760E+02, 3.043E+02,8.611E+01, 7.868E+01, 3.408E+00, 1.680E-05/
      DATA (PH1(I,28,28, 4),I=1,6) /1.250E+02, 5.448E+01,4.611E+00, 1.157E+02, 4.548E+00, 0.000E+00/
      DATA (PH1(I,28,28, 5),I=1,6) /8.200E+01, 8.896E+01,4.351E+01, 1.942E+01, 7.372E+00, 2.566E-01/
      DATA (PH1(I,28,28, 6),I=1,6) /1.700E+01, 6.063E+00,1.186E+03, 6.823E+00, 2.223E+01, 6.227E-03/
      DATA (PH1(I,28,28, 7),I=1,6) /7.637E+00, 1.468E+01,1.437E+00, 7.411E+02, 4.342E+00, 3.908E-02/
      DATA (PH1(I,29, 1, 1),I=1,6) /1.157E+04, 3.656E+02,6.510E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,29, 2, 1),I=1,6) /1.106E+04, 1.427E+03,8.057E+00, 7.912E+01, 1.583E+00, 0.000E+00/
      DATA (PH1(I,29, 3, 1),I=1,6) /1.085E+04, 1.435E+03,8.196E+00, 7.379E+01, 1.577E+00, 0.000E+00/
      DATA (PH1(I,29, 3, 2),I=1,6) /2.585E+03, 1.732E+02,6.336E+00, 4.099E+01, 4.316E+00, 0.000E+00/
      DATA (PH1(I,29, 4, 1),I=1,6) /1.070E+04, 1.476E+03,7.658E+00, 7.281E+01, 1.567E+00, 0.000E+00/
      DATA (PH1(I,29, 4, 2),I=1,6) /2.459E+03, 1.822E+02,1.108E+01, 3.417E+01, 4.489E+00, 0.000E+00/
      DATA (PH1(I,29, 5, 1),I=1,6) /1.053E+04, 1.516E+03,7.156E+00, 7.157E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,29, 5, 2),I=1,6) /2.363E+03, 2.174E+02,9.074E+00, 2.892E+01, 4.445E+00, 0.000E+00/
      DATA (PH1(I,29, 5, 3),I=1,6) /2.298E+03, 1.099E+02,1.898E+02, 1.272E+02, 4.081E+00, 2.821E-08/
      DATA (PH1(I,29, 6, 1),I=1,6) /1.036E+04, 1.537E+03,6.899E+00, 7.045E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,29, 6, 2),I=1,6) /2.253E+03, 2.118E+02,8.876E+00, 2.857E+01, 4.513E+00, 0.000E+00/
      DATA (PH1(I,29, 6, 3),I=1,6) /2.173E+03, 1.503E+02,1.866E+02, 1.055E+02, 3.944E+00, 1.972E-08/
      DATA (PH1(I,29, 7, 1),I=1,6) /1.016E+04, 1.531E+03,7.047E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29, 7, 2),I=1,6) /2.133E+03, 2.330E+02,7.838E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29, 7, 3),I=1,6) /2.045E+03, 2.473E+02,9.111E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,29, 8, 1),I=1,6) /1.004E+04, 1.556E+03,6.650E+00, 6.987E+01, 1.559E+00, 0.000E+00/
      DATA (PH1(I,29, 8, 2),I=1,6) /2.037E+03, 2.003E+02,8.422E+00, 2.803E+01, 4.659E+00, 0.000E+00/
      DATA (PH1(I,29, 8, 3),I=1,6) /1.905E+03, 1.999E+02,1.798E+02, 7.715E+01, 3.888E+00, 1.496E-08/
      DATA (PH1(I,29, 9, 1),I=1,6) /9.844E+03, 1.546E+03,6.833E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29, 9, 2),I=1,6) /1.920E+03, 2.388E+02,6.970E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29, 9, 3),I=1,6) /1.793E+03, 2.493E+02,1.335E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,29,10, 1),I=1,6) /9.729E+03, 1.540E+03,6.748E+00, 6.512E+01, 1.590E+00, 0.000E+00/
      DATA (PH1(I,29,10, 2),I=1,6) /1.827E+03, 1.901E+02,7.838E+00, 2.795E+01, 4.791E+00, 0.000E+00/
      DATA (PH1(I,29,10, 3),I=1,6) /1.690E+03, 2.206E+02,1.906E+02, 6.390E+01, 3.955E+00, 1.558E-08/
      DATA (PH1(I,29,11, 1),I=1,6) /9.668E+03, 1.523E+03,6.915E+00, 6.870E+01, 1.581E+00, 0.000E+00/
      DATA (PH1(I,29,11, 2),I=1,6) /1.773E+03, 1.626E+02,8.723E+00, 3.155E+01, 4.881E+00, 0.000E+00/
      DATA (PH1(I,29,11, 3),I=1,6) /1.657E+03, 2.651E+02,1.307E+02, 7.793E+01, 3.637E+00, 1.894E-08/
      DATA (PH1(I,29,11, 4),I=1,6) /6.706E+02, 1.529E+01,2.453E+01, 6.172E+02, 4.260E+00, 0.000E+00/
      DATA (PH1(I,29,12, 1),I=1,6) /9.610E+03, 1.558E+03,6.567E+00, 6.667E+01, 1.576E+00, 0.000E+00/
      DATA (PH1(I,29,12, 2),I=1,6) /1.722E+03, 1.953E+02,7.587E+00, 2.834E+01, 4.738E+00, 0.000E+00/
      DATA (PH1(I,29,12, 3),I=1,6) /1.604E+03, 2.619E+02,1.331E+02, 7.654E+01, 3.661E+00, 1.879E-08/
      DATA (PH1(I,29,12, 4),I=1,6) /6.330E+02, 2.604E+01,2.426E+01, 3.618E+02, 4.240E+00, 0.000E+00/
      DATA (PH1(I,29,13, 1),I=1,6) /9.550E+03, 1.563E+03,6.515E+00, 6.616E+01, 1.577E+00, 0.000E+00/
      DATA (PH1(I,29,13, 2),I=1,6) /1.670E+03, 1.993E+02,7.415E+00, 2.757E+01, 4.743E+00, 0.000E+00/
      DATA (PH1(I,29,13, 3),I=1,6) /1.551E+03, 2.401E+02,1.577E+02, 7.468E+01, 3.767E+00, 2.341E-08/
      DATA (PH1(I,29,13, 4),I=1,6) /5.949E+02, 3.706E+01,1.653E+01, 2.918E+02, 4.111E+00, 0.000E+00/
      DATA (PH1(I,29,13, 5),I=1,6) /5.570E+02, 2.011E+02,3.929E+00, 5.061E+01, 4.854E+00, 8.398E-01/
      DATA (PH1(I,29,14, 1),I=1,6) /9.493E+03, 1.558E+03,6.561E+00, 6.691E+01, 1.575E+00, 0.000E+00/
      DATA (PH1(I,29,14, 2),I=1,6) /1.619E+03, 1.974E+02,7.417E+00, 2.691E+01, 4.791E+00, 0.000E+00/
      DATA (PH1(I,29,14, 3),I=1,6) /1.499E+03, 2.384E+02,1.588E+02, 7.405E+01, 3.780E+00, 2.405E-08/
      DATA (PH1(I,29,14, 4),I=1,6) /5.572E+02, 3.531E+01,1.668E+01, 3.042E+02, 4.130E+00, 0.000E+00/
      DATA (PH1(I,29,14, 5),I=1,6) /5.200E+02, 1.207E+02,1.569E+01, 5.849E+01, 5.405E+00, 4.555E-01/
      DATA (PH1(I,29,15, 1),I=1,6) /9.423E+03, 1.577E+03,6.419E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,15, 2),I=1,6) /1.563E+03, 2.446E+02,6.221E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,15, 3),I=1,6) /1.439E+03, 2.632E+02,1.311E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,15, 4),I=1,6) /5.184E+02, 5.001E+01,1.019E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,15, 5),I=1,6) /4.840E+02, 1.032E+02,2.783E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,16, 1),I=1,6) /9.384E+03, 1.553E+03,6.609E+00, 6.923E+01, 1.568E+00, 0.000E+00/
      DATA (PH1(I,29,16, 2),I=1,6) /1.522E+03, 1.990E+02,7.271E+00, 2.604E+01, 4.826E+00, 0.000E+00/
      DATA (PH1(I,29,16, 3),I=1,6) /1.400E+03, 2.302E+02,1.675E+02, 7.119E+01, 3.847E+00, 2.441E-08/
      DATA (PH1(I,29,16, 4),I=1,6) /4.838E+02, 3.424E+01,1.441E+01, 2.579E+02, 4.283E+00, 0.000E+00/
      DATA (PH1(I,29,16, 5),I=1,6) /4.350E+02, 1.064E+02,3.399E+01, 5.025E+01, 5.723E+00, 4.731E-01/
      DATA (PH1(I,29,17, 1),I=1,6) /9.317E+03, 1.576E+03,6.433E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,17, 2),I=1,6) /1.467E+03, 2.457E+02,6.110E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,17, 3),I=1,6) /1.340E+03, 2.635E+02,1.295E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,17, 4),I=1,6) /4.458E+02, 5.099E+01,9.072E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,17, 5),I=1,6) /4.010E+02, 1.031E+02,4.266E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,18, 1),I=1,6) /9.282E+03, 1.550E+03,6.633E+00, 6.922E+01, 1.570E+00, 0.000E+00/
      DATA (PH1(I,29,18, 2),I=1,6) /1.431E+03, 1.960E+02,7.274E+00, 2.502E+01, 4.908E+00, 0.000E+00/
      DATA (PH1(I,29,18, 3),I=1,6) /1.307E+03, 2.291E+02,1.666E+02, 6.860E+01, 3.881E+00, 2.421E-08/
      DATA (PH1(I,29,18, 4),I=1,6) /4.129E+02, 3.356E+01,1.254E+01, 2.441E+02, 4.373E+00, 0.000E+00/
      DATA (PH1(I,29,18, 5),I=1,6) /3.688E+02, 1.023E+02,4.888E+01, 4.483E+01, 5.919E+00, 3.996E-01/
      DATA (PH1(I,29,19, 1),I=1,6) /9.237E+03, 1.550E+03,6.630E+00, 6.917E+01, 1.570E+00, 0.000E+00/
      DATA (PH1(I,29,19, 2),I=1,6) /1.386E+03, 1.983E+02,7.172E+00, 2.486E+01, 4.901E+00, 0.000E+00/
      DATA (PH1(I,29,19, 3),I=1,6) /1.262E+03, 2.322E+02,1.615E+02, 6.669E+01, 3.889E+00, 2.403E-08/
      DATA (PH1(I,29,19, 4),I=1,6) /3.778E+02, 3.308E+01,1.152E+01, 2.347E+02, 4.436E+00, 0.000E+00/
      DATA (PH1(I,29,19, 5),I=1,6) /3.335E+02, 1.017E+02,4.646E+01, 4.207E+01, 6.020E+00, 3.632E-01/
      DATA (PH1(I,29,19, 6),I=1,6) /2.661E+02, 4.605E+01,1.458E+02, 7.967E+01, 6.285E+00, 4.503E-03/
      DATA (PH1(I,29,20, 1),I=1,6) /9.181E+03, 1.574E+03,6.454E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,20, 2),I=1,6) /1.339E+03, 2.470E+02,5.993E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,20, 3),I=1,6) /1.211E+03, 2.640E+02,1.274E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,20, 4),I=1,6) /3.447E+02, 5.350E+01,7.318E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,20, 5),I=1,6) /2.999E+02, 1.045E+02,4.286E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,20, 6),I=1,6) /2.320E+02, 5.210E+01,2.093E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,21, 1),I=1,6) /9.155E+03, 1.546E+03,6.671E+00, 6.964E+01, 1.570E+00, 0.000E+00/
      DATA (PH1(I,29,21, 2),I=1,6) /1.302E+03, 1.958E+02,7.186E+00, 2.423E+01, 4.960E+00, 0.000E+00/
      DATA (PH1(I,29,21, 3),I=1,6) /1.178E+03, 2.383E+02,1.525E+02, 6.417E+01, 3.890E+00, 2.517E-08/
      DATA (PH1(I,29,21, 4),I=1,6) /3.113E+02, 4.346E+01,7.773E+00, 1.823E+02, 4.413E+00, 0.000E+00/
      DATA (PH1(I,29,21, 5),I=1,6) /2.672E+02, 9.821E+01,4.314E+01, 3.667E+01, 6.291E+00, 2.759E-01/
      DATA (PH1(I,29,21, 6),I=1,6) /1.990E+02, 5.000E+01,3.041E+02, 6.355E+01, 6.458E+00, 2.040E-03/
      DATA (PH1(I,29,22, 1),I=1,6) /9.127E+03, 1.573E+03,6.468E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,22, 2),I=1,6) /1.267E+03, 2.476E+02,5.944E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,22, 3),I=1,6) /1.137E+03, 2.644E+02,1.263E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,22, 4),I=1,6) /2.824E+02, 5.587E+01,6.260E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,22, 5),I=1,6) /2.375E+02, 1.064E+02,3.728E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,22, 6),I=1,6) /1.670E+02, 5.236E+01,3.470E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,23, 1),I=1,6) /9.107E+03, 1.541E+03,6.731E+00, 7.196E+01, 1.562E+00, 0.000E+00/
      DATA (PH1(I,29,23, 2),I=1,6) /1.229E+03, 1.968E+02,7.138E+00, 2.379E+01, 4.979E+00, 0.000E+00/
      DATA (PH1(I,29,23, 3),I=1,6) /1.104E+03, 2.376E+02,1.529E+02, 6.163E+01, 3.925E+00, 2.511E-08/
      DATA (PH1(I,29,23, 4),I=1,6) /2.504E+02, 5.085E+01,5.992E+00, 1.563E+02, 4.410E+00, 0.000E+00/
      DATA (PH1(I,29,23, 5),I=1,6) /2.067E+02, 1.010E+02,3.860E+01, 3.121E+01, 6.491E+00, 3.331E-01/
      DATA (PH1(I,29,23, 6),I=1,6) /1.390E+02, 5.652E+01,3.452E+02, 5.007E+01, 6.558E+00, 2.171E-03/
      DATA (PH1(I,29,24, 1),I=1,6) /9.088E+03, 1.572E+03,6.482E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,24, 2),I=1,6) /1.203E+03, 2.479E+02,5.918E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,24, 3),I=1,6) /1.073E+03, 2.649E+02,1.255E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,24, 4),I=1,6) /2.243E+02, 5.878E+01,5.342E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,24, 5),I=1,6) /1.801E+02, 1.092E+02,3.217E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,24, 6),I=1,6) /1.030E+02, 5.346E+01,4.010E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,25, 1),I=1,6) /9.070E+03, 1.572E+03,6.440E+00, 6.913E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,29,25, 2),I=1,6) /1.166E+03, 2.072E+02,6.884E+00, 2.337E+01, 4.917E+00, 0.000E+00/
      DATA (PH1(I,29,25, 3),I=1,6) /1.041E+03, 2.368E+02,1.537E+02, 5.951E+01, 3.955E+00, 2.525E-08/
      DATA (PH1(I,29,25, 4),I=1,6) /1.960E+02, 5.776E+01,4.896E+00, 1.391E+02, 4.397E+00, 0.000E+00/
      DATA (PH1(I,29,25, 5),I=1,6) /1.528E+02, 9.650E+01,3.804E+01, 2.618E+01, 6.870E+00, 2.460E-01/
      DATA (PH1(I,29,25, 6),I=1,6) /7.990E+01, 5.875E+01,3.868E+02, 3.774E+01, 6.878E+00, 3.091E-01/
      DATA (PH1(I,29,26, 1),I=1,6) /9.052E+03, 1.571E+03,6.497E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,26, 2),I=1,6) /1.146E+03, 2.480E+02,5.914E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,26, 3),I=1,6) /1.020E+03, 2.654E+02,1.248E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,26, 4),I=1,6) /1.703E+02, 6.225E+01,4.565E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,26, 5),I=1,6) /1.277E+02, 1.127E+02,2.769E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,26, 6),I=1,6) /5.738E+01, 5.539E+01,4.028E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,27, 1),I=1,6) /9.032E+03, 2.025E+03,3.761E+00, 2.259E+02, 1.184E+00, 0.000E+00/
      DATA (PH1(I,29,27, 2),I=1,6) /1.128E+03, 2.508E+02,5.929E+00, 2.384E+01, 4.581E+00, 0.000E+00/
      DATA (PH1(I,29,27, 3),I=1,6) /9.911E+02, 2.369E+02,1.537E+02, 5.847E+01, 3.967E+00, 2.525E-08/
      DATA (PH1(I,29,27, 4),I=1,6) /1.499E+02, 6.247E+01,4.311E+00, 1.291E+02, 4.390E+00, 0.000E+00/
      DATA (PH1(I,29,27, 5),I=1,6) /1.072E+02, 9.745E+01,3.768E+01, 2.211E+01, 7.130E+00, 2.537E-01/
      DATA (PH1(I,29,27, 6),I=1,6) /3.684E+01, 3.541E+01,1.181E+03, 1.390E+01, 1.017E+01, 3.416E-01/
      DATA (PH1(I,29,28, 1),I=1,6) /9.012E+03, 2.378E+03,2.659E+00, 7.489E+02, 9.937E-01, 0.000E+00/
      DATA (PH1(I,29,28, 2),I=1,6) /1.114E+03, 3.054E+02,4.889E+00, 2.675E+01, 4.155E+00, 0.000E+00/
      DATA (PH1(I,29,28, 3),I=1,6) /9.715E+02, 2.328E+02,1.588E+02, 5.776E+01, 3.997E+00, 2.543E-08/
      DATA (PH1(I,29,28, 4),I=1,6) /1.312E+02, 6.472E+01,4.139E+00, 1.278E+02, 4.366E+00, 0.000E+00/
      DATA (PH1(I,29,28, 5),I=1,6) /8.861E+01, 9.750E+01,3.865E+01, 2.022E+01, 7.277E+00, 2.537E-01/
      DATA (PH1(I,29,28, 6),I=1,6) /2.029E+01, 1.590E+01,3.992E+03, 7.098E+00, 1.616E+01, 7.439E-01/
      DATA (PH1(I,29,29, 1),I=1,6) /8.988E+03, 1.788E+03,4.870E+00, 7.645E+01, 1.451E+00, 0.000E+00/
      DATA (PH1(I,29,29, 2),I=1,6) /1.106E+03, 2.502E+02,5.938E+00, 2.402E+01, 4.576E+00, 0.000E+00/
      DATA (PH1(I,29,29, 3),I=1,6) /9.470E+02, 2.826E+02,1.093E+02, 6.688E+01, 3.668E+00, 9.174E-07/
      DATA (PH1(I,29,29, 4),I=1,6) /1.288E+02, 6.198E+01,4.164E+00, 1.158E+02, 4.488E+00, 0.000E+00/
      DATA (PH1(I,29,29, 5),I=1,6) /8.300E+01, 9.617E+01,4.275E+01, 1.747E+01, 7.555E+00, 2.599E-01/
      DATA (PH1(I,29,29, 6),I=1,6) /1.064E+01, 7.279E+00,1.027E+03, 7.988E+00, 2.033E+01, 1.582E+00/
      DATA (PH1(I,29,29, 7),I=1,6) /7.726E+00, 1.436E+01,4.681E-01, 2.383E+03, 4.224E+00, 3.736E-02/
      DATA (PH1(I,30, 1, 1),I=1,6) /1.239E+04, 3.915E+02,6.083E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,30, 2, 1),I=1,6) /1.187E+04, 1.564E+03,7.165E+00, 9.095E+01, 1.535E+00, 0.000E+00/
      DATA (PH1(I,30, 3, 1),I=1,6) /1.165E+04, 1.600E+03,7.032E+00, 9.272E+01, 1.497E+00, 0.000E+00/
      DATA (PH1(I,30, 3, 2),I=1,6) /2.780E+03, 3.219E+02,3.342E+00, 3.904E+01, 3.684E+00, 0.000E+00/
      DATA (PH1(I,30, 4, 1),I=1,6) /1.149E+04, 1.565E+03,7.294E+00, 7.835E+01, 1.557E+00, 0.000E+00/
      DATA (PH1(I,30, 4, 2),I=1,6) /2.647E+03, 1.989E+02,1.025E+01, 3.502E+01, 4.434E+00, 0.000E+00/
      DATA (PH1(I,30, 5, 1),I=1,6) /1.131E+04, 1.582E+03,7.052E+00, 7.246E+01, 1.577E+00, 0.000E+00/
      DATA (PH1(I,30, 5, 2),I=1,6) /2.550E+03, 2.168E+02,9.068E+00, 2.977E+01, 4.514E+00, 0.000E+00/
      DATA (PH1(I,30, 5, 3),I=1,6) /2.479E+03, 1.230E+02,1.663E+02, 1.330E+02, 4.013E+00, 2.250E-08/
      DATA (PH1(I,30, 6, 1),I=1,6) /1.113E+04, 1.607E+03,6.765E+00, 7.173E+01, 1.574E+00, 0.000E+00/
      DATA (PH1(I,30, 6, 2),I=1,6) /2.435E+03, 2.335E+02,8.158E+00, 2.769E+01, 4.506E+00, 0.000E+00/
      DATA (PH1(I,30, 6, 3),I=1,6) /2.363E+03, 1.386E+02,2.396E+02, 1.169E+02, 4.016E+00, 1.579E-08/
      DATA (PH1(I,30, 7, 1),I=1,6) /1.092E+04, 1.637E+03,6.617E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30, 7, 2),I=1,6) /2.309E+03, 2.488E+02,7.424E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30, 7, 3),I=1,6) /2.216E+03, 2.651E+02,8.601E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,30, 8, 1),I=1,6) /1.080E+04, 1.629E+03,6.509E+00, 7.085E+01, 1.574E+00, 0.000E+00/
      DATA (PH1(I,30, 8, 2),I=1,6) /2.210E+03, 2.179E+02,7.856E+00, 2.852E+01, 4.612E+00, 0.000E+00/
      DATA (PH1(I,30, 8, 3),I=1,6) /2.070E+03, 2.034E+02,1.893E+02, 8.099E+01, 3.908E+00, 6.713E-09/
      DATA (PH1(I,30, 9, 1),I=1,6) /1.059E+04, 1.652E+03,6.423E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30, 9, 2),I=1,6) /2.087E+03, 2.548E+02,6.629E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30, 9, 3),I=1,6) /1.953E+03, 2.670E+02,1.269E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,30,10, 1),I=1,6) /1.048E+04, 1.636E+03,6.400E+00, 6.864E+01, 1.584E+00, 0.000E+00/
      DATA (PH1(I,30,10, 2),I=1,6) /1.992E+03, 2.104E+02,7.257E+00, 2.813E+01, 4.728E+00, 0.000E+00/
      DATA (PH1(I,30,10, 3),I=1,6) /1.846E+03, 2.254E+02,1.992E+02, 6.645E+01, 3.978E+00, 6.413E-09/
      DATA (PH1(I,30,11, 1),I=1,6) /1.042E+04, 1.605E+03,6.675E+00, 7.053E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,30,11, 2),I=1,6) /1.936E+03, 1.953E+02,7.654E+00, 3.173E+01, 4.697E+00, 0.000E+00/
      DATA (PH1(I,30,11, 3),I=1,6) /1.815E+03, 2.273E+02,1.928E+02, 6.347E+01, 4.011E+00, 6.017E-09/
      DATA (PH1(I,30,11, 4),I=1,6) /7.374E+02, 1.960E+01,1.939E+01, 5.293E+02, 4.225E+00, 0.000E+00/
      DATA (PH1(I,30,12, 1),I=1,6) /1.035E+04, 1.631E+03,6.432E+00, 6.791E+01, 1.590E+00, 0.000E+00/
      DATA (PH1(I,30,12, 2),I=1,6) /1.882E+03, 2.070E+02,7.256E+00, 2.865E+01, 4.736E+00, 0.000E+00/
      DATA (PH1(I,30,12, 3),I=1,6) /1.760E+03, 2.265E+02,1.930E+02, 6.308E+01, 4.019E+00, 6.018E-09/
      DATA (PH1(I,30,12, 4),I=1,6) /6.980E+02, 3.115E+01,2.108E+01, 3.355E+02, 4.197E+00, 0.000E+00/
      DATA (PH1(I,30,13, 1),I=1,6) /1.029E+04, 1.636E+03,6.381E+00, 6.693E+01, 1.592E+00, 0.000E+00/
      DATA (PH1(I,30,13, 2),I=1,6) /1.828E+03, 2.085E+02,7.155E+00, 2.796E+01, 4.758E+00, 0.000E+00/
      DATA (PH1(I,30,13, 3),I=1,6) /1.704E+03, 2.296E+02,1.871E+02, 6.376E+01, 3.995E+00, 5.868E-09/
      DATA (PH1(I,30,13, 4),I=1,6) /6.583E+02, 3.116E+01,2.073E+01, 3.448E+02, 4.181E+00, 0.000E+00/
      DATA (PH1(I,30,13, 5),I=1,6) /6.190E+02, 2.539E+02,2.886E+00, 5.384E+01, 4.592E+00, 8.255E-01/
      DATA (PH1(I,30,14, 1),I=1,6) /1.023E+04, 1.633E+03,6.399E+00, 6.753E+01, 1.591E+00, 0.000E+00/
      DATA (PH1(I,30,14, 2),I=1,6) /1.775E+03, 2.101E+02,7.058E+00, 2.748E+01, 4.771E+00, 0.000E+00/
      DATA (PH1(I,30,14, 3),I=1,6) /1.650E+03, 2.720E+02,1.333E+02, 8.071E+01, 3.661E+00, 4.757E-09/
      DATA (PH1(I,30,14, 4),I=1,6) /6.187E+02, 2.985E+01,2.083E+01, 3.577E+02, 4.196E+00, 0.000E+00/
      DATA (PH1(I,30,14, 5),I=1,6) /5.790E+02, 1.403E+02,1.361E+01, 6.070E+01, 5.261E+00, 5.034E-01/
      DATA (PH1(I,30,15, 1),I=1,6) /1.015E+04, 1.691E+03,5.972E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,15, 2),I=1,6) /1.716E+03, 2.609E+02,5.912E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,15, 3),I=1,6) /1.586E+03, 2.806E+02,1.269E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,15, 4),I=1,6) /5.780E+02, 5.311E+01,1.000E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,15, 5),I=1,6) /5.420E+02, 1.108E+02,2.676E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,16, 1),I=1,6) /1.012E+04, 1.617E+03,6.536E+00, 6.839E+01, 1.594E+00, 0.000E+00/
      DATA (PH1(I,30,16, 2),I=1,6) /1.673E+03, 2.142E+02,6.875E+00, 2.631E+01, 4.800E+00, 0.000E+00/
      DATA (PH1(I,30,16, 3),I=1,6) /1.546E+03, 2.596E+02,1.443E+02, 7.722E+01, 3.739E+00, 5.496E-09/
      DATA (PH1(I,30,16, 4),I=1,6) /5.416E+02, 3.912E+01,1.395E+01, 2.883E+02, 4.155E+00, 0.000E+00/
      DATA (PH1(I,30,16, 5),I=1,6) /4.900E+02, 1.153E+02,3.254E+01, 5.308E+01, 5.645E+00, 4.056E-01/
      DATA (PH1(I,30,17, 1),I=1,6) /1.004E+04, 1.689E+03,5.993E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,17, 2),I=1,6) /1.613E+03, 2.621E+02,5.810E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,17, 3),I=1,6) /1.481E+03, 2.807E+02,1.254E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,17, 4),I=1,6) /5.015E+02, 5.403E+01,8.971E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,17, 5),I=1,6) /4.540E+02, 1.106E+02,4.142E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,18, 1),I=1,6) /1.001E+04, 1.617E+03,6.532E+00, 6.768E+01, 1.597E+00, 0.000E+00/
      DATA (PH1(I,30,18, 2),I=1,6) /1.576E+03, 2.135E+02,6.802E+00, 2.575E+01, 4.839E+00, 0.000E+00/
      DATA (PH1(I,30,18, 3),I=1,6) /1.448E+03, 2.538E+02,1.487E+02, 7.400E+01, 3.794E+00, 5.551E-09/
      DATA (PH1(I,30,18, 4),I=1,6) /4.671E+02, 3.780E+01,1.186E+01, 2.389E+02, 4.334E+00, 0.000E+00/
      DATA (PH1(I,30,18, 5),I=1,6) /4.197E+02, 1.073E+02,4.863E+01, 4.776E+01, 5.875E+00, 2.432E-01/
      DATA (PH1(I,30,19, 1),I=1,6) /9.960E+03, 1.610E+03,6.593E+00, 6.774E+01, 1.600E+00, 0.000E+00/
      DATA (PH1(I,30,19, 2),I=1,6) /1.528E+03, 2.118E+02,6.804E+00, 2.536E+01, 4.875E+00, 0.000E+00/
      DATA (PH1(I,30,19, 3),I=1,6) /1.399E+03, 2.529E+02,1.489E+02, 7.156E+01, 3.823E+00, 5.330E-09/
      DATA (PH1(I,30,19, 4),I=1,6) /4.299E+02, 3.741E+01,1.098E+01, 2.313E+02, 4.386E+00, 0.000E+00/
      DATA (PH1(I,30,19, 5),I=1,6) /3.828E+02, 1.057E+02,4.664E+01, 4.477E+01, 5.995E+00, 2.234E-01/
      DATA (PH1(I,30,19, 6),I=1,6) /3.108E+02, 5.749E+01,1.122E+02, 7.804E+01, 6.097E+00, 3.731E-01/
      DATA (PH1(I,30,20, 1),I=1,6) /9.915E+03, 1.611E+03,6.591E+00, 6.799E+01, 1.598E+00, 0.000E+00/
      DATA (PH1(I,30,20, 2),I=1,6) /1.482E+03, 2.023E+02,6.987E+00, 2.500E+01, 4.972E+00, 0.000E+00/
      DATA (PH1(I,30,20, 3),I=1,6) /1.353E+03, 2.543E+02,1.465E+02, 6.961E+01, 3.838E+00, 5.027E-09/
      DATA (PH1(I,30,20, 4),I=1,6) /3.939E+02, 3.705E+01,1.008E+01, 2.216E+02, 4.447E+00, 0.000E+00/
      DATA (PH1(I,30,20, 5),I=1,6) /3.467E+02, 1.046E+02,4.454E+01, 4.180E+01, 6.115E+00, 1.997E-01/
      DATA (PH1(I,30,20, 6),I=1,6) /2.740E+02, 5.735E+01,2.033E+02, 7.143E+01, 6.214E+00, 3.888E-01/
      DATA (PH1(I,30,21, 1),I=1,6) /9.854E+03, 1.685E+03,6.036E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,21, 2),I=1,6) /1.436E+03, 2.638E+02,5.673E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,21, 3),I=1,6) /1.302E+03, 2.813E+02,1.231E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,21, 4),I=1,6) /3.609E+02, 5.754E+01,6.797E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,21, 5),I=1,6) /3.128E+02, 1.126E+02,3.961E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,21, 6),I=1,6) /2.380E+02, 5.727E+01,2.777E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,30,22, 1),I=1,6) /9.831E+03, 1.605E+03,6.653E+00, 6.874E+01, 1.597E+00, 0.000E+00/
      DATA (PH1(I,30,22, 2),I=1,6) /1.397E+03, 2.114E+02,6.748E+00, 2.443E+01, 4.934E+00, 0.000E+00/
      DATA (PH1(I,30,22, 3),I=1,6) /1.267E+03, 2.512E+02,1.490E+02, 6.562E+01, 3.896E+00, 5.943E-09/
      DATA (PH1(I,30,22, 4),I=1,6) /3.257E+02, 4.721E+01,7.066E+00, 1.761E+02, 4.431E+00, 0.000E+00/
      DATA (PH1(I,30,22, 5),I=1,6) /2.785E+02, 1.016E+02,4.118E+01, 3.606E+01, 6.396E+00, 1.450E-01/
      DATA (PH1(I,30,22, 6),I=1,6) /2.030E+02, 5.906E+01,3.187E+02, 5.850E+01, 6.430E+00, 2.472E-01/
      DATA (PH1(I,30,23, 1),I=1,6) /9.809E+03, 1.683E+03,6.057E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,23, 2),I=1,6) /1.362E+03, 2.642E+02,5.637E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,23, 3),I=1,6) /1.227E+03, 2.816E+02,1.222E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,23, 4),I=1,6) /2.969E+02, 6.012E+01,5.837E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,23, 5),I=1,6) /2.487E+02, 1.147E+02,3.456E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,23, 6),I=1,6) /1.750E+02, 5.760E+01,3.804E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,30,24, 1),I=1,6) /9.788E+03, 1.697E+03,5.935E+00, 9.066E+01, 1.488E+00, 0.000E+00/
      DATA (PH1(I,30,24, 2),I=1,6) /1.323E+03, 2.094E+02,6.772E+00, 2.401E+01, 4.975E+00, 0.000E+00/
      DATA (PH1(I,30,24, 3),I=1,6) /1.193E+03, 2.522E+02,1.476E+02, 6.376E+01, 3.914E+00, 6.079E-09/
      DATA (PH1(I,30,24, 4),I=1,6) /2.633E+02, 5.492E+01,5.518E+00, 1.512E+02, 4.431E+00, 0.000E+00/
      DATA (PH1(I,30,24, 5),I=1,6) /2.165E+02, 1.096E+02,3.576E+01, 3.034E+01, 6.524E+00, 3.502E-01/
      DATA (PH1(I,30,24, 6),I=1,6) /1.360E+02, 6.530E+01,3.456E+02, 4.608E+01, 6.555E+00, 2.476E-01/
      DATA (PH1(I,30,25, 1),I=1,6) /9.770E+03, 1.682E+03,6.079E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,25, 2),I=1,6) /1.298E+03, 2.645E+02,5.621E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,25, 3),I=1,6) /1.162E+03, 2.821E+02,1.215E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,25, 4),I=1,6) /2.369E+02, 6.326E+01,5.006E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,25, 5),I=1,6) /1.896E+02, 1.177E+02,2.996E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,25, 6),I=1,6) /1.080E+02, 5.878E+01,4.151E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,30,26, 1),I=1,6) /9.751E+03, 1.732E+03,5.666E+00, 8.751E+01, 1.484E+00, 0.000E+00/
      DATA (PH1(I,30,26, 2),I=1,6) /1.259E+03, 2.231E+02,6.480E+00, 2.364E+01, 4.891E+00, 0.000E+00/
      DATA (PH1(I,30,26, 3),I=1,6) /1.129E+03, 2.523E+02,1.473E+02, 6.205E+01, 3.933E+00, 6.078E-09/
      DATA (PH1(I,30,26, 4),I=1,6) /2.075E+02, 6.183E+01,4.578E+00, 1.352E+02, 4.422E+00, 0.000E+00/
      DATA (PH1(I,30,26, 5),I=1,6) /1.612E+02, 1.030E+02,3.598E+01, 2.554E+01, 6.932E+00, 2.426E-01/
      DATA (PH1(I,30,26, 6),I=1,6) /8.260E+01, 6.378E+01,4.172E+02, 3.462E+01, 7.005E+00, 3.108E-01/
      DATA (PH1(I,30,27, 1),I=1,6) /9.733E+03, 1.890E+03,4.708E+00, 1.202E+02, 1.363E+00, 0.000E+00/
      DATA (PH1(I,30,27, 2),I=1,6) /1.239E+03, 2.366E+02,6.210E+00, 2.347E+01, 4.803E+00, 0.000E+00/
      DATA (PH1(I,30,27, 3),I=1,6) /1.101E+03, 2.528E+02,1.467E+02, 6.144E+01, 3.938E+00, 6.071E-09/
      DATA (PH1(I,30,27, 4),I=1,6) /1.826E+02, 6.427E+01,4.285E+00, 1.295E+02, 4.424E+00, 0.000E+00/
      DATA (PH1(I,30,27, 5),I=1,6) /1.365E+02, 1.033E+02,3.582E+01, 2.337E+01, 7.076E+00, 2.470E-01/
      DATA (PH1(I,30,27, 6),I=1,6) /5.940E+01, 5.376E+01,6.007E+02, 2.429E+01, 7.985E+00, 3.026E-01/
      DATA (PH1(I,30,28, 1),I=1,6) /9.713E+03, 1.981E+03,4.228E+00, 1.057E+02, 1.361E+00, 0.000E+00/
      DATA (PH1(I,30,28, 2),I=1,6) /1.222E+03, 2.686E+02,5.600E+00, 2.394E+01, 4.574E+00, 0.000E+00/
      DATA (PH1(I,30,28, 3),I=1,6) /1.077E+03, 2.519E+02,1.478E+02, 6.090E+01, 3.949E+00, 6.087E-09/
      DATA (PH1(I,30,28, 4),I=1,6) /1.601E+02, 6.684E+01,4.057E+00, 1.255E+02, 4.414E+00, 0.000E+00/
      DATA (PH1(I,30,28, 5),I=1,6) /1.143E+02, 1.045E+02,3.568E+01, 2.181E+01, 7.163E+00, 2.533E-01/
      DATA (PH1(I,30,28, 6),I=1,6) /3.972E+01, 3.640E+01,1.498E+03, 1.188E+01, 1.082E+01, 3.519E-01/
      DATA (PH1(I,30,29, 1),I=1,6) /9.691E+03, 1.125E+03,1.382E+01, 2.765E+01, 2.242E+00, 0.000E+00/
      DATA (PH1(I,30,29, 2),I=1,6) /1.209E+03, 2.121E+02,6.713E+00, 2.342E+01, 4.988E+00, 0.000E+00/
      DATA (PH1(I,30,29, 3),I=1,6) /1.065E+03, 2.664E+02,1.329E+02, 6.260E+01, 3.864E+00, 5.913E-09/
      DATA (PH1(I,30,29, 4),I=1,6) /1.479E+02, 6.443E+01,4.075E+00, 1.175E+02, 4.506E+00, 0.000E+00/
      DATA (PH1(I,30,29, 5),I=1,6) /1.021E+02, 1.047E+02,3.618E+01, 2.122E+01, 7.201E+00, 2.534E-01/
      DATA (PH1(I,30,29, 6),I=1,6) /2.694E+01, 2.593E+01,3.408E+03, 7.832E+00, 1.363E+01, 3.562E-01/
      DATA (PH1(I,30,29, 7),I=1,6) /1.796E+01, 1.411E+01,8.168E-01, 8.775E+02, 4.393E+00, 4.021E-04/
      DATA (PH1(I,30,30, 1),I=1,6) /9.667E+03, 8.320E+02,2.586E+01, 4.497E+01, 2.215E+00, 0.000E+00/
      DATA (PH1(I,30,30, 2),I=1,6) /1.203E+03, 9.755E+01,9.077E+00, 3.219E+01, 5.888E+00, 0.000E+00/
      DATA (PH1(I,30,30, 3),I=1,6) /1.037E+03, 3.486E+02,7.784E+01, 8.298E+01, 3.391E+00, 1.125E-08/
      DATA (PH1(I,30,30, 4),I=1,6) /1.450E+02, 6.195E+01,4.094E+00, 1.086E+02, 4.614E+00, 0.000E+00/
      DATA (PH1(I,30,30, 5),I=1,6) /9.700E+01, 1.025E+02,3.930E+01, 1.876E+01, 7.456E+00, 2.559E-01/
      DATA (PH1(I,30,30, 6),I=1,6) /1.730E+01, 1.818E+01,1.017E+04, 5.288E+00, 1.736E+01, 4.667E-01/
      DATA (PH1(I,30,30, 7),I=1,6) /9.394E+00, 1.673E+01,1.236E+00, 1.029E+03, 4.259E+00, 3.962E-02/
      DATA (PH2(I, 1, 1),I=1,7) /4.298E-01, 5.475E+04,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 2, 1),I=1,7) /1.720E+00, 1.369E+04,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 2, 2),I=1,7) /1.361E+01, 9.492E+02,1.469E+00, 3.188E+00, 2.039E+00, 4.434E-01, 2.136E+00/
      DATA (PH2(I, 3, 1),I=1,7) /3.871E+00, 6.083E+03,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 3, 2),I=1,7) /2.006E+01, 3.201E+02,7.391E+00, 2.916E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 3, 3),I=1,7) /3.107E+00, 6.245E+01,1.501E+01, 4.895E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 1),I=1,7) /6.879E+00, 3.422E+03,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 2),I=1,7) /1.760E+01, 5.458E+02,1.719E+01, 3.157E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 3),I=1,7) /1.181E+00, 2.678E+02,5.645E+00, 1.170E+01, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 4),I=1,7) /9.539E+00, 2.932E+05,4.301E-01, 1.052E+01, 3.655E-01, 8.278E-04, 1.269E-02/
      DATA (PH2(I, 5, 1),I=1,7) /1.075E+01, 2.190E+03,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 5, 2),I=1,7) /3.336E+01, 2.846E+02,2.163E+01, 2.624E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 5, 3),I=1,7) /1.041E+00, 5.393E+01,1.767E+01, 9.540E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 5, 4),I=1,7) /2.869E+00, 1.859E+04,1.783E+00, 1.618E+01, 3.503E+00, 4.960E-03, 3.400E-02/
      DATA (PH2(I, 5, 5),I=1,7) /5.213E-01, 5.466E+00,8.618E+00, 1.728E+01, 1.887E+01, 1.319E+01, 4.556E+00/
      DATA (PH2(I, 6, 1),I=1,7) /1.548E+01, 1.521E+03,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 6, 2),I=1,7) /4.624E+01, 2.344E+02,2.183E+01, 2.581E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 6, 3),I=1,7) /3.506E+00, 1.068E+02,1.436E+01, 7.457E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 6, 4),I=1,7) /4.614E+00, 1.539E+04,1.737E+00, 1.593E+01, 5.922E+00, 4.378E-03, 2.528E-02/
      DATA (PH2(I, 6, 5),I=1,7) /4.058E-01, 8.709E+00,1.261E+02, 8.578E+00, 2.093E+00, 4.929E+01, 3.234E+00/
      DATA (PH2(I, 6, 6),I=1,7) /2.144E+00, 5.027E+02,6.216E+01, 5.101E+00, 9.157E-02, 1.133E+00, 1.607E+00/
      DATA (PH2(I, 7, 1),I=1,7) /2.108E+01, 1.117E+03,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 7, 2),I=1,7) /6.943E+01, 1.519E+02,2.627E+01, 2.315E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 7, 3),I=1,7) /4.471E+00, 8.376E+01,3.297E+01, 6.003E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 7, 4),I=1,7) /5.494E+00, 1.690E+04,1.714E+00, 1.706E+01, 7.904E+00, 6.415E-03, 1.937E-02/
      DATA (PH2(I, 7, 5),I=1,7) /2.420E-01, 9.375E-01,2.788E+02, 9.156E+00, 1.850E+00, 1.877E+02, 3.999E+00/
      DATA (PH2(I, 7, 6),I=1,7) /6.128E-02, 1.944E+00,8.163E+02, 8.773E+00, 1.043E+01, 4.280E+02, 2.030E+01/
      DATA (PH2(I, 7, 7),I=1,7) /4.034E+00, 8.235E+02,8.033E+01, 3.928E+00, 9.097E-02, 8.598E-01, 2.325E+00/
      DATA (PH2(I, 8, 1),I=1,7) /2.754E+01, 8.554E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 8, 2),I=1,7) /8.709E+01, 1.329E+02,2.535E+01, 2.336E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 8, 3),I=1,7) /7.824E+00, 6.864E+01,3.210E+01, 5.495E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 8, 4),I=1,7) /2.854E+00, 1.642E+04,1.792E+00, 2.647E+01, 2.836E+01, 3.036E-02, 5.554E-02/
      DATA (PH2(I, 8, 5),I=1,7) /2.044E-01, 8.659E-01,4.931E+02, 8.785E+00, 3.143E+00, 3.328E+02, 4.285E+01/
      DATA (PH2(I, 8, 6),I=1,7) /1.723E-01, 6.753E+02,3.852E+02, 6.822E+00, 1.191E-01, 3.839E-03, 4.569E-01/
      DATA (PH2(I, 8, 7),I=1,7) /1.386E+00, 5.967E+01,3.175E+01, 8.943E+00, 1.934E-02, 2.131E+01, 1.503E-02/
      DATA (PH2(I, 8, 8),I=1,7) /1.240E+00, 1.745E+03,3.784E+00, 1.764E+01, 7.589E-02, 8.698E+00, 1.271E-01/
      DATA (PH2(I, 9, 1),I=1,7) /3.485E+01, 6.759E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 9, 2),I=1,7) /1.131E+02, 1.039E+02,2.657E+01, 2.255E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 9, 3),I=1,7) /2.563E+00, 6.930E+01,7.547E+01, 6.448E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 9, 4),I=1,7) /4.008E+00, 1.157E+04,1.848E+00, 2.446E+01, 2.411E+01, 2.071E-02, 3.998E-02/
      DATA (PH2(I, 9, 5),I=1,7) /7.286E-01, 4.690E-01,1.400E+02, 9.718E+00, 2.570E-01, 1.506E+02, 2.574E-01/
      DATA (PH2(I, 9, 6),I=1,7) /7.744E-01, 3.165E+00,1.099E+02, 9.203E+00, 6.812E+00, 9.531E+01, 9.781E+00/
      DATA (PH2(I, 9, 7),I=1,7) /2.542E+00, 1.541E+02,5.742E+01, 6.614E+00, 1.115E+00, 1.641E+01, 5.124E+00/
      DATA (PH2(I, 9, 8),I=1,7) /1.763E+00, 8.013E+01,1.667E+01, 1.050E+01, 5.103E-01, 1.715E+01, 7.724E-01/
      DATA (PH2(I, 9, 9),I=1,7) /1.297E+01, 3.803E+03,2.587E+00, 7.275E+00, 2.170E-03, 1.701E-04, 1.345E-02/
      DATA (PH2(I,10, 1),I=1,7) /4.304E+01, 5.475E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,10, 2),I=1,7) /1.586E+02, 6.695E+01,3.352E+01, 2.002E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,10, 3),I=1,7) /1.003E+01, 5.631E+01,3.628E+01, 5.585E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,10, 4),I=1,7) /4.888E+00, 1.198E+04,1.788E+00, 2.550E+01, 2.811E+01, 2.536E-02, 4.417E-02/
      DATA (PH2(I,10, 5),I=1,7) /1.499E+00, 9.854E-01,1.350E+02, 8.836E+00, 1.656E+00, 1.042E+02, 1.435E+00/
      DATA (PH2(I,10, 6),I=1,7) /1.248E+00, 2.430E+00,1.066E+02, 8.999E+00, 6.855E-01, 9.169E+01, 3.702E-01/
      DATA (PH2(I,10, 7),I=1,7) /5.566E+00, 1.685E+03,6.409E+02, 3.056E+00, 8.290E-03, 5.149E+00, 6.687E+00/
      DATA (PH2(I,10, 8),I=1,7) /7.753E-01, 5.708E+00,6.725E+01, 1.005E+01, 4.633E-01, 7.654E+01, 2.023E+00/
      DATA (PH2(I,10, 9),I=1,7) /1.247E+01, 1.583E+03,3.935E+00, 7.810E+00, 6.558E-02, 1.520E+00, 1.084E-01/
      DATA (PH2(I,10,10),I=1,7) /4.870E+00, 4.287E+03,5.798E+00, 8.355E+00, 2.434E-01, 4.236E-02, 5.873E+00/
      DATA (PH2(I,11, 1),I=1,7) /5.211E+01, 4.525E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,11, 2),I=1,7) /2.268E+02, 3.995E+01,5.315E+01, 1.678E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,11, 3),I=1,7) /1.391E+01, 4.729E+01,3.889E+01, 5.265E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,11, 4),I=1,7) /1.535E+02, 7.215E+03,3.886E-01, 8.476E+00, 9.121E-01, 1.667E-01, 1.766E-02/
      DATA (PH2(I,11, 5),I=1,7) /2.096E+00, 1.609E+00,2.473E+02, 7.681E+00, 1.895E+00, 9.940E+01, 3.278E+00/
      DATA (PH2(I,11, 6),I=1,7) /4.846E+01, 7.101E+01,3.945E+01, 2.832E+00, 1.285E-02, 9.603E-04, 6.378E-03/
      DATA (PH2(I,11, 7),I=1,7) /5.408E+00, 2.346E+01,2.913E+01, 8.260E+00, 9.275E-01, 2.204E+01, 7.577E-01/
      DATA (PH2(I,11, 8),I=1,7) /6.690E-01, 2.330E+00,1.205E+02, 9.714E+00, 7.365E-01, 1.383E+02, 4.260E+00/
      DATA (PH2(I,11, 9),I=1,7) /1.069E+01, 1.885E+03,3.613E+00, 9.803E+00, 8.579E-02, 3.725E+00, 2.279E-01/
      DATA (PH2(I,11,10),I=1,7) /8.203E+00, 1.040E+03,8.259E+00, 7.362E+00, 2.328E+00, 3.375E+00, 4.010E+00/
      DATA (PH2(I,11,11),I=1,7) /6.139E+00, 1.601E+00,6.148E+03, 3.839E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 1),I=1,7) /6.203E+01, 3.802E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 2),I=1,7) /2.042E+02, 6.140E+01,2.778E+01, 2.161E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 3),I=1,7) /1.452E+01, 4.427E+01,3.826E+01, 5.460E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 4),I=1,7) /3.482E+01, 9.008E+02,1.823E+00, 1.444E+01, 2.751E+00, 5.444E+00, 7.918E-02/
      DATA (PH2(I,12, 5),I=1,7) /4.884E-01, 6.344E-02,5.085E+02, 9.385E+00, 6.666E-01, 5.348E+02, 3.997E-03/
      DATA (PH2(I,12, 6),I=1,7) /3.570E+00, 3.104E+00,6.060E+01, 8.857E+00, 1.422E+00, 5.452E+01, 2.078E+00/
      DATA (PH2(I,12, 7),I=1,7) /1.711E+00, 2.185E+00,9.350E+01, 9.202E+00, 6.325E-01, 1.007E+02, 1.729E+00/
      DATA (PH2(I,12, 8),I=1,7) /9.762E-01, 1.728E+00,9.184E+01, 1.006E+01, 8.090E-01, 1.276E+02, 3.979E+00/
      DATA (PH2(I,12, 9),I=1,7) /2.912E+01, 1.394E+03,2.895E+00, 6.487E+00, 4.326E-02, 9.402E-01, 1.135E-01/
      DATA (PH2(I,12,10),I=1,7) /1.086E+01, 5.377E+02,9.779E+00, 7.117E+00, 2.604E+00, 4.860E+00, 3.722E+00/
      DATA (PH2(I,12,11),I=1,7) /8.139E+00, 3.278E+00,4.341E+07, 3.610E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12,12),I=1,7) /1.197E+01, 1.372E+08,2.228E-01, 1.574E+01, 2.805E-01, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 1),I=1,7) /7.281E+01, 3.239E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 2),I=1,7) /2.738E+02, 4.036E+01,3.567E+01, 1.915E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 3),I=1,7) /2.355E+01, 3.388E+01,3.432E+01, 5.085E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 4),I=1,7) /8.044E+00, 1.774E+04,1.653E+00, 2.655E+01, 2.953E+01, 2.538E-02, 1.203E-02/
      DATA (PH2(I,13, 5),I=1,7) /1.842E+00, 4.982E-01,2.568E+02, 8.406E+00, 6.945E-01, 1.719E+02, 6.595E+00/
      DATA (PH2(I,13, 6),I=1,7) /4.866E-01, 2.350E-01,7.216E+02, 8.659E+00, 2.773E-01, 5.704E+02, 1.580E-01/
      DATA (PH2(I,13, 7),I=1,7) /2.636E+00, 1.889E+02,1.338E+02, 6.204E+00, 1.836E+00, 3.552E+01, 8.223E-03/
      DATA (PH2(I,13, 8),I=1,7) /3.483E-01, 1.962E-02,1.856E+01, 2.084E+01, 8.839E+00, 5.675E-02, 2.768E-01/
      DATA (PH2(I,13, 9),I=1,7) /2.414E+01, 2.925E+02,6.973E+00, 6.724E+00, 1.000E-01, 3.495E+00, 2.701E-01/
      DATA (PH2(I,13,10),I=1,7) /3.130E+00, 1.513E+01,1.674E+01, 1.180E+01, 5.342E+00, 3.994E+01, 4.803E+00/
      DATA (PH2(I,13,11),I=1,7) /1.027E+01, 4.915E+00,1.990E+06, 3.477E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13,12),I=1,7) /2.048E-01, 6.948E-02,5.675E+02, 9.049E+00, 4.615E-01, 9.149E+01, 6.565E-01/
      DATA (PH2(I,13,13),I=1,7) /1.381E+01, 7.195E+00,1.621E+03, 3.642E+00, 3.166E-01, 2.041E-01, 4.753E-01/
      DATA (PH2(I,14, 1),I=1,7) /8.447E+01, 2.793E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14, 2),I=1,7) /2.752E+02, 4.754E+01,2.848E+01, 2.135E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14, 3),I=1,7) /3.560E+01, 2.539E+01,3.307E+01, 4.728E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14, 4),I=1,7) /1.205E+01, 1.992E+04,1.582E+00, 2.425E+01, 2.392E+01, 1.990E-02, 1.007E-02/
      DATA (PH2(I,14, 5),I=1,7) /8.787E-01, 1.950E-01,7.461E+02, 8.302E+00, 4.489E-01, 4.528E+02, 1.015E+00/
      DATA (PH2(I,14, 6),I=1,7) /3.343E-01, 1.465E-01,1.404E+03, 8.503E+00, 1.646E+00, 1.036E+03, 2.936E-01/
      DATA (PH2(I,14, 7),I=1,7) /7.655E-01, 3.477E-01,3.733E+02, 8.986E+00, 1.476E-03, 3.850E+02, 8.999E-02/
      DATA (PH2(I,14, 8),I=1,7) /3.277E-01, 6.680E-02,4.132E+01, 1.606E+01, 3.280E+00, 1.149E-02, 6.396E-01/
      DATA (PH2(I,14, 9),I=1,7) /6.305E+01, 7.293E+01,1.558E+02, 2.400E+00, 2.989E-03, 1.115E+00, 8.051E-02/
      DATA (PH2(I,14,10),I=1,7) /7.761E-01, 8.863E-01,1.541E+02, 9.980E+00, 1.303E+00, 2.009E+02, 4.537E+00/
      DATA (PH2(I,14,11),I=1,7) /1.288E+01, 6.083E+00,1.356E+06, 3.353E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14,12),I=1,7) /1.659E-01, 5.790E-04,1.474E+02, 1.336E+01, 8.626E-01, 9.613E+01, 6.442E-01/
      DATA (PH2(I,14,13),I=1,7) /2.556E+00, 4.140E+00,1.337E+01, 1.191E+01, 1.570E+00, 6.634E+00, 1.272E-01/
      DATA (PH2(I,14,14),I=1,7) /2.317E+01, 2.506E+01,2.057E+01, 3.546E+00, 2.837E-01, 1.672E-05, 4.207E-01/
      DATA (PH2(I,16, 1),I=1,7) /1.104E+02, 2.139E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16, 2),I=1,7) /4.390E+02, 2.453E+01,4.405E+01, 1.765E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16, 3),I=1,7) /3.310E+01, 2.555E+01,3.821E+01, 5.037E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16, 4),I=1,7) /1.474E+01, 2.294E+04,1.529E+00, 2.568E+01, 2.738E+01, 2.203E-02, 1.073E-02/
      DATA (PH2(I,16, 5),I=1,7) /2.443E+00, 3.490E-01,5.411E+02, 7.769E+00, 7.033E-01, 2.279E+02, 1.172E+00/
      DATA (PH2(I,16, 6),I=1,7) /6.485E+00, 1.275E+01,6.583E+01, 7.692E+00, 1.678E+00, 3.426E+01, 1.370E-01/
      DATA (PH2(I,16, 7),I=1,7) /1.040E+01, 5.364E+01,3.641E+01, 7.090E+00, 2.310E+00, 1.775E+01, 1.663E+00/
      DATA (PH2(I,16, 8),I=1,7) /1.526E-01, 9.646E+03,1.438E+03, 5.977E+00, 1.492E+00, 1.615E-03, 4.049E-01/
      DATA (PH2(I,16, 9),I=1,7) /1.462E+01, 3.161E+01,1.611E+01, 8.642E+00, 1.153E-03, 1.869E+01, 3.037E-01/
      DATA (PH2(I,16,10),I=1,7) /3.757E-01, 5.703E-01,1.460E+02, 1.135E+01, 1.503E+00, 2.222E+02, 4.606E+00/
      DATA (PH2(I,16,11),I=1,7) /1.413E+01, 9.139E+00,1.656E+03, 3.626E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16,12),I=1,7) /1.713E-01, 5.072E-04,1.986E+02, 1.307E+01, 7.880E-01, 9.424E+01, 6.265E-01/
      DATA (PH2(I,16,13),I=1,7) /2.173E+00, 2.606E+00,6.641E+01, 8.655E+00, 1.863E+00, 1.975E+01, 3.361E+00/
      DATA (PH2(I,16,14),I=1,7) /2.027E+00, 6.666E+00,5.454E+01, 8.611E+00, 4.109E+00, 1.568E+01, 9.421E+00/
      DATA (PH2(I,16,15),I=1,7) /8.787E+00, 3.136E+02,3.442E+00, 1.281E+01, 7.354E-01, 2.782E+00, 1.788E-01/
      DATA (PH2(I,16,16),I=1,7) /1.808E+01, 4.564E+04,1.000E+00, 1.361E+01, 6.385E-01, 9.935E-01, 2.486E-01/
      DATA (PH2(I,18, 1),I=1,7) /1.399E+02, 1.690E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18, 2),I=1,7) /4.468E+02, 3.108E+01,3.039E+01, 2.092E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18, 3),I=1,7) /4.154E+01, 2.135E+01,4.118E+01, 4.945E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18, 4),I=1,7) /1.888E+01, 2.571E+04,1.475E+00, 2.634E+01, 2.909E+01, 2.445E-02, 1.054E-02/
      DATA (PH2(I,18, 5),I=1,7) /1.557E+00, 4.997E-02,5.031E+02, 8.966E+00, 2.938E-01, 4.552E+02, 6.459E+00/
      DATA (PH2(I,18, 6),I=1,7) /3.209E-01, 2.459E-02,2.285E+03, 8.810E+00, 6.692E-01, 2.068E+03, 2.113E+01/
      DATA (PH2(I,18, 7),I=1,7) /5.310E+00, 7.018E-01,1.001E+02, 8.939E+00, 4.987E-01, 1.099E+02, 2.202E-01/
      DATA (PH2(I,18, 8),I=1,7) /1.257E-01, 1.760E+03,1.579E+03, 6.714E+00, 1.975E+00, 3.286E-03, 3.226E-01/
      DATA (PH2(I,18, 9),I=1,7) /1.040E+01, 8.204E+00,1.495E+01, 1.115E+01, 9.203E-04, 3.804E+01, 6.390E-01/
      DATA (PH2(I,18,10),I=1,7) /1.926E-01, 8.279E-01,2.392E+02, 1.121E+01, 1.434E+00, 3.814E+01, 4.649E+00/
      DATA (PH2(I,18,11),I=1,7) /3.884E+00, 3.295E+01,7.082E+02, 4.645E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18,12),I=1,7) /2.966E-02, 3.693E+00,9.951E+03, 7.313E+00, 1.363E-02, 4.383E-04, 2.513E+00/
      DATA (PH2(I,18,13),I=1,7) /5.440E-01, 1.080E+00,9.419E+02, 7.582E+00, 1.107E+01, 1.700E+02, 1.587E+01/
      DATA (PH2(I,18,14),I=1,7) /1.031E+01, 9.946E+00,7.444E+01, 6.261E+00, 4.885E-01, 6.406E+00, 3.659E-03/
      DATA (PH2(I,18,15),I=1,7) /6.953E+00, 2.035E+01,1.400E+01, 9.595E+00, 8.842E-01, 7.501E+00, 1.806E-01/
      DATA (PH2(I,18,16),I=1,7) /1.417E+01, 3.580E+01,3.776E+01, 5.742E+00, 6.316E-01, 2.384E+00, 1.794E+00/
      DATA (PH2(I,18,17),I=1,7) /2.494E+01, 2.503E+01,1.272E+02, 4.288E+00, 5.108E-01, 9.299E-01, 7.195E-01/
      DATA (PH2(I,18,18),I=1,7) /1.709E+01, 2.106E+01,2.645E+02, 4.796E+00, 4.185E-01, 1.688E+00, 8.943E-01/
      DATA (PH2(I,20, 1),I=1,7) /1.729E+02, 1.369E+02,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20, 2),I=1,7) /6.297E+02, 1.936E+01,3.921E+01, 1.862E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20, 3),I=1,7) /9.472E+01, 1.105E+01,3.818E+01, 4.192E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20, 4),I=1,7) /2.618E+01, 2.028E+04,1.456E+00, 2.560E+01, 2.803E+01, 2.402E-02, 9.323E-03/
      DATA (PH2(I,20, 5),I=1,7) /4.293E+00, 1.293E+00,1.691E+01, 1.438E+01, 3.461E-05, 9.363E-01, 4.589E-02/
      DATA (PH2(I,20, 6),I=1,7) /1.309E+02, 5.513E+01,3.828E+02, 2.023E+00, 9.084E-02, 1.833E-02, 9.359E-01/
      DATA (PH2(I,20, 7),I=1,7) /9.980E+00, 1.116E+00,5.918E+01, 9.005E+00, 3.879E+00, 7.104E+01, 5.311E+00/
      DATA (PH2(I,20, 8),I=1,7) /1.008E+01, 1.849E+03,1.792E+04, 2.868E+00, 2.410E+02, 6.138E-03, 6.931E+01/
      DATA (PH2(I,20, 9),I=1,7) /2.345E+01, 1.227E+01,1.312E+01, 9.771E+00, 6.842E-04, 2.417E+01, 5.469E-01/
      DATA (PH2(I,20,10),I=1,7) /2.288E-01, 9.384E-01,2.549E+02, 1.103E+01, 1.390E+00, 2.478E+01, 3.100E+00/
      DATA (PH2(I,20,11),I=1,7) /1.605E+01, 1.437E+01,6.989E+02, 3.857E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20,12),I=1,7) /5.520E-02, 2.076E+02,1.790E+04, 5.893E+00, 1.843E-03, 2.826E-04, 1.657E+00/
      DATA (PH2(I,20,13),I=1,7) /1.366E+00, 6.641E-01,3.188E+02, 8.138E+00, 2.806E-01, 1.039E+02, 3.329E+00/
      DATA (PH2(I,20,14),I=1,7) /8.080E-01, 4.760E-01,3.682E+02, 8.634E+00, 5.720E-01, 1.487E+02, 1.283E+00/
      DATA (PH2(I,20,15),I=1,7) /9.515E+00, 7.642E+01,8.973E+01, 5.141E+00, 2.471E+00, 4.829E+00, 5.824E+00/
      DATA (PH2(I,20,16),I=1,7) /6.882E-01, 1.523E-01,1.502E+02, 1.061E+01, 8.227E+00, 1.210E+02, 3.876E+00/
      DATA (PH2(I,20,17),I=1,7) /4.255E+00, 7.736E+00,1.355E+01, 1.236E+01, 1.369E+00, 1.467E+01, 3.298E-02/
      DATA (PH2(I,20,18),I=1,7) /2.436E+01, 3.815E+01,2.931E+02, 3.944E+00, 3.126E-01, 1.802E+00, 1.233E+00/
      DATA (PH2(I,20,19),I=1,7) /1.553E+01, 1.064E+07,7.790E-01, 2.130E+01, 6.453E-01, 2.161E-03, 6.706E-02/
      DATA (PH2(I,20,20),I=1,7) /1.278E+01, 5.370E+05,3.162E-01, 1.242E+01, 4.477E-01, 1.012E-03, 1.851E-02/
      DATA (PH2(I,26, 1),I=1,7) /2.932E+02, 8.099E+01,3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26, 2),I=1,7) /1.057E+03, 1.195E+01,5.769E+01, 1.718E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26, 3),I=1,7) /7.326E+01, 1.276E+01,4.914E+01, 4.941E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26, 4),I=1,7) /4.575E+01, 2.580E+04,1.358E+00, 2.604E+01, 2.723E+01, 3.582E-02, 8.712E-03/
      DATA (PH2(I,26, 5),I=1,7) /9.713E+00, 7.204E-02,1.853E+02, 8.843E+00, 9.551E-03, 1.702E+02, 4.263E+00/
      DATA (PH2(I,26, 6),I=1,7) /9.243E+00, 1.098E+01,7.637E+01, 7.962E+00, 1.748E+00, 4.446E+01, 3.512E+00/
      DATA (PH2(I,26, 7),I=1,7) /2.011E+01, 4.455E-01,4.236E+01, 9.724E+00, 2.757E+00, 6.847E+01, 3.989E+00/
      DATA (PH2(I,26, 8),I=1,7) /7.519E-04, 6.066E-05,1.606E+06, 8.813E+00, 4.398E+00, 1.915E+06, 3.140E+01/
      DATA (PH2(I,26, 9),I=1,7) /3.190E+01, 2.388E+00,2.186E+01, 9.589E+00, 2.902E-02, 3.805E+01, 4.805E-01/
      DATA (PH2(I,26,10),I=1,7) /3.444E-01, 1.452E+00,3.960E+02, 1.013E+01, 1.264E+00, 2.891E+01, 3.404E+00/
      DATA (PH2(I,26,11),I=1,7) /2.873E+01, 1.207E+01,5.150E+02, 3.846E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26,12),I=1,7) /5.555E-02, 2.108E+02,2.045E+04, 6.033E+00, 1.885E-03, 2.706E-04, 1.628E+00/
      DATA (PH2(I,26,13),I=1,7) /8.509E-01, 1.454E-01,1.239E+03, 8.066E+00, 4.937E-01, 4.505E+02, 2.504E+00/
      DATA (PH2(I,26,14),I=1,7) /1.317E-01, 2.791E-03,2.487E+03, 9.791E+00, 6.938E-01, 2.170E+03, 6.852E-03/
      DATA (PH2(I,26,15),I=1,7) /6.295E+00, 1.738E+00,1.130E+02, 8.037E+00, 3.096E-01, 4.671E+01, 1.425E-01/
      DATA (PH2(I,26,16),I=1,7) /8.284E+00, 3.281E+00,5.360E+01, 8.571E+00, 3.279E-01, 2.971E+01, 5.220E-01/
      DATA (PH2(I,26,17),I=1,7) /6.886E+01, 6.470E+01,2.062E+01, 4.111E+00, 2.778E-04, 1.190E-05, 6.570E-03/
      DATA (PH2(I,26,18),I=1,7) /6.741E+00, 2.687E+01,1.807E+02, 6.290E+00, 2.387E-04, 2.494E+01, 8.251E+00/
      DATA (PH2(I,26,19),I=1,7) /7.098E-02, 1.979E+01,1.745E+04, 6.750E+00, 2.158E+02, 2.542E+03, 4.672E+02/
      DATA (PH2(I,26,20),I=1,7) /5.059E+00, 2.420E+04,4.850E+04, 2.374E+00, 2.516E-03, 4.546E-01, 2.683E+01/
      DATA (PH2(I,26,21),I=1,7) /2.656E+00, 5.259E-01,1.450E+01, 1.632E+01, 1.558E+01, 3.361E+01, 3.743E-03/
      DATA (PH2(I,26,22),I=1,7) /7.256E-01, 1.523E-03,3.736E+01, 1.767E+01, 5.064E+01, 8.871E+01, 5.280E-02/
      DATA (PH2(I,26,23),I=1,7) /2.544E+01, 3.653E+02,8.913E+00, 6.538E+00, 5.602E-01, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26,24),I=1,7) /1.698E-01, 6.107E+00,1.555E+03, 8.055E+00, 8.698E+00, 1.760E+02, 1.847E+01/
      DATA (PH2(I,26,25),I=1,7) /1.761E-01, 4.365E+03,6.298E+03, 5.204E+00, 1.141E+01, 9.272E+01, 1.075E+02/
      DATA (PH2(I,26,26),I=1,7) /5.461E-02, 3.062E-01,2.671E+07, 7.923E+00, 2.069E+01, 1.382E+02, 2.481E-01/
      END

     SUBROUTINE readAugerKM93(W_Auger_tmp)

!     read in the modified table for the fitting parameters for direct
!     collisional ionization from the work of Kaastra + Mewe
!     (1993.A&ASS,97,443)

!     we convert the Xray shells indices of Kaastra to the
!     photoionization shell indices of Verner because we will need to
!     compute the Auger rates from the partial photoionization rates for
!     each shell and these rates are computed using the Verner indices

!     the translation is as follows:
!
!     Kaastra                Verner
!     idx  Xray  Term        s     shell  g
!      1   K     1s_1/2      1     1s     
!      2   L1    2s_1/2      2     2s
!      3   L2    2p^5_1/2    3     2p     1  
!      4   L3    2p_5_3/2    3     2p     2   
!      5   M1    3s_1/2      4     3s
!      6   M2    3p5_1/2     5     3p     1
!      7   M3    3p5_3/2     5     3p     2
!      8   M4    3d9_3/2     6     3d
!      9   M5    3d9_5/2     6     3d
!     10   N1    4s          7     4s
!
!     where g is the statistical weight for the averaging; Kaastra's
!     table goes no higher than idx=7

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      implicit none

      include           'rates.h'
      integer,parameter :: Nrec=1090
      integer :: i,ii,j,k,n,s,nline,Xshnum,z
      REAL (Kind=4) :: W(Imax,Imax,Nshmax,10)
      REAL (Kind=4) :: pr(10),E,EA,eps
      character(LEN=80) :: header
      !character(LEN=255) :: tabfile
!      character (len=240) :: local_dir ='/media/Data/Projects2/ioneqther/v7'
      character (*), parameter :: data_file = '/tab-auger-Wijk.dat'
 
      character(LEN=500) :: tabfile 
!      include           'com-auger.com'
!      include           'com-file.com'
      integer           jAugktot,nAugshtot,Xshell
      COMMON /iauger/   jAugktot(Imax),nAugshtot(Imax,Imax), &
                       Xshell(Imax,Imax,Nshmax)
      double precision :: W_Auger_tmp(Imax,Imax,Nshmax,10) 
      double precision  W_Auger
      COMMON /rauger/   W_Auger(Imax,Imax,Nshmax,10) 
      character (len=240) :: local_dir = '/media/efrain/DATA/softwares/modelosXSPEC/IonEq/thermal/solar/v1.1'
      tabfile = trim(local_dir)//trim(data_file)
 
!     construct file name
!   RN   Z   st   s     I        EA   epsilon  PrEj1  PrEj2  PrEj3  PrEj4  PrEj5  PrEj6  PrEj7  PrEj8  PrEj9 PrEj10
!      tabfile = homedir
!      tabfile = local_dir,'tab-auger-Wijk.dat'
!      tabfile=  tabfile
 !     CALL sappend(tabfile,'tab-auger-Wijk.dat',tabfile)

!     read file

      OPEN(unit=2,file=tabfile,status='old')
      READ(2,*) header
      DO 11 i=1,Nrec
       READ(2,*,END=12) Nline,z,j,Xshnum,E,EA,eps,(pr(ii),ii=1,10)
       k = z
       CALL store_auger(k,j,Xshnum,pr,W)
 11   CONTINUE
 12   CLOSE(unit=2)

!     translate the Xray shell numbers to Verner's shell numbers so we
!     can calculate the Auger rates using the photoionization cross
!     sections from Verner

!     the translation and weighting is given in the header comments of
!     this routine

      DO 21 k=1,Imax
       DO 22 j=1,jAugktot(k)
        DO 23 n=1,nAugshtot(k,j)

         IF (Xshell(k,j,n).eq.1) then
          s = 1
          DO 31 i=1,10
           W_auger(k,j,s,i) = W(k,j,n,i) 
           W_Auger_tmp(k,j,s,i) = W_auger(k,j,s,i)
 31       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.2) then
          s = 2
          DO 32 i=1,10
           W_auger(k,j,s,i) = W(k,j,n,i) 
           W_Auger_tmp(k,j,s,i) = W_auger(k,j,s,i)
 32       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.3) then
          s = 3
          DO 33 i=1,10
           W_auger(k,j,s,i) = (W(k,j,n,i) + 2.0d0*W(k,j,n+1,i))/3.0d0
           W_Auger_tmp(k,j,s,i) = W_auger(k,j,s,i)
 33       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.5) then
          s = 4
          DO 34 i=1,10
           W_auger(k,j,s,i) = W(k,j,n,i) 
           W_Auger_tmp(k,j,s,i) = W_auger(k,j,s,i)
 34       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.6) then
          s = 5
          DO 35 i=1,10
           W_auger(k,j,s,i) = (W(k,j,n,i) + 2.0d0*W(k,j,n+1,i))/3.0d0
           W_Auger_tmp(k,j,s,i) = W_auger(k,j,s,i)
 35       CONTINUE
         END IF

 23     CONTINUE
 22    CONTINUE
 21   CONTINUE


      RETURN

      END


!
!.............................................................................
!

      SUBROUTINE store_auger(k,j,Xshnum,pr,W)

!     this routine simply stores the read in line into the common block
!     array and matrices for the given k,j combination

!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

      implicit none
      include           'rates.h'
      integer :: i,j,k,n,Xshnum
      real :: W(Imax,Imax,Nshmax,10)
      real :: pr(10)
      include           'com-auger.com'


!     store the fitting parameters and subshell data

      n  = nAugshtot(k,j) + 1

      Xshell(k,j,n) = Xshnum

      DO 11 i=1,10
       W(k,j,n,i) = 0.0001*pr(i)
 11   CONTINUE


!     save the shell counters for k,j

      nAugshtot(k,j) = n
      jAugktot(k)    = max(jAugktot(k),j)

      RETURN
      END


!==========================================================================
      block data tabl_auger
!--------------------------------------------------------------------------
!     AUGER COEFFICIENTS DATA Kaastra + Mewe (1993.A&ASS,97,443)
!--------------------------------------------------------------------------
      integer, parameter :: ndr=1090
      real :: tab_auger(17,ndr)
      common/tab_auger/tab_auger
!  Z   st   s     I        EA   epsilon  PrEj1  PrEj2  PrEj3  PrEj4  PrEj5  PrEj6  PrEj7  PrEj8  PrEj9 PrEj10
data (tab_auger(i,1),i=1,16) /4,1,1,115,92.8,46,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,2),i=1,16) /5,1,1,192,163.9,91,6,9994,0,0,0,0,0,0,0,0/
data (tab_auger(i,3),i=1,16) /5,2,1,206,149.9,34,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,4),i=1,16) /6,1,1,288,254.1,39,26,9974,0,0,0,0,0,0,0,0/
data (tab_auger(i,5),i=1,16) /6,2,1,306,240.4,73,19,9981,0,0,0,0,0,0,0,0/
data (tab_auger(i,6),i=1,16) /6,3,1,325,221.3,27,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,7),i=1,16) /7,1,1,403,365.2,48,60,9940,0,0,0,0,0,0,0,0/
data (tab_auger(i,8),i=1,16) /7,2,1,426,348.9,32,71,9929,0,0,0,0,0,0,0,0/
data (tab_auger(i,9),i=1,16) /7,3,1,448,329.7,62,52,9948,0,0,0,0,0,0,0,0/
data (tab_auger(i,10),i=1,16) /7,4,1,471,305.9,22,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,11),i=1,16) /8,1,1,538,491.6,22,94,9906,0,0,0,0,0,0,0,0/
data (tab_auger(i,12),i=1,16) /8,2,1,565,480.1,43,109,9891,0,0,0,0,0,0,0,0/
data (tab_auger(i,13),i=1,16) /8,3,1,592,456.2,27,129,9871,0,0,0,0,0,0,0,0/
data (tab_auger(i,14),i=1,16) /8,4,1,618,432.2,53,96,9904,0,0,0,0,0,0,0,0/
data (tab_auger(i,15),i=1,16) /8,5,1,645,404.9,18,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,16),i=1,16) /9,1,1,694,635,17,133,9867,0,0,0,0,0,0,0,0/
data (tab_auger(i,17),i=1,16) /9,2,1,724,624.7,19,152,9848,0,0,0,0,0,0,0,0/
data (tab_auger(i,18),i=1,16) /9,3,1,754,605.3,38,176,9824,0,0,0,0,0,0,0,0/
data (tab_auger(i,19),i=1,16) /9,4,1,785,573.9,23,208,9792,0,0,0,0,0,0,0,0/
data (tab_auger(i,20),i=1,16) /9,5,1,815,546,47,154,9846,0,0,0,0,0,0,0,0/
data (tab_auger(i,21),i=1,16) /9,6,1,846,517.1,16,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,22),i=1,16) /10,1,1,870.1,787.4,7,182,9818,0,0,0,0,0,0,0,0/
data (tab_auger(i,23),i=1,16) /10,2,1,903,783.7,15,204,9796,0,0,0,0,0,0,0,0/
data (tab_auger(i,24),i=1,16) /10,3,1,937,766.9,17,232,9768,0,0,0,0,0,0,0,0/
data (tab_auger(i,25),i=1,16) /10,4,1,971,740.5,33,268,9732,0,0,0,0,0,0,0,0/
data (tab_auger(i,26),i=1,16) /10,5,1,1005,700.7,21,312,9688,0,0,0,0,0,0,0,0/
data (tab_auger(i,27),i=1,16) /10,6,1,1039,670.1,42,229,9771,0,0,0,0,0,0,0,0/
data (tab_auger(i,28),i=1,16) /10,7,1,1074,642.6,14,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,29),i=1,16) /11,1,1,1075,953.6,10,260,6342,3398,0,0,0,0,0,0,0/
data (tab_auger(i,30),i=1,16) /11,1,2,67,25.3,160,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,31),i=1,16) /11,2,1,1101,948.4,4,263,9737,0,0,0,0,0,0,0,0/
data (tab_auger(i,32),i=1,16) /11,3,1,1139,941.8,14,296,9704,0,0,0,0,0,0,0,0/
data (tab_auger(i,33),i=1,16) /11,4,1,1177,917.8,15,338,9662,0,0,0,0,0,0,0,0/
data (tab_auger(i,34),i=1,16) /11,5,1,1215,883.7,30,393,9607,0,0,0,0,0,0,0,0/
data (tab_auger(i,35),i=1,16) /11,6,1,1253,834.4,18,466,9534,0,0,0,0,0,0,0,0/
data (tab_auger(i,36),i=1,16) /11,7,1,1291,802.6,36,352,9648,0,0,0,0,0,0,0,0/
data (tab_auger(i,37),i=1,16) /11,8,1,1329,781.4,12,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,38),i=1,16) /12,1,1,1308,1149.7,7,1,544,8645,810,0,0,0,0,0,0/
data (tab_auger(i,39),i=1,16) /12,1,2,92,27.7,87,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,40),i=1,16) /12,1,3,54.3,35.9,80,33,9967,0,0,0,0,0,0,0,0/
data (tab_auger(i,41),i=1,16) /12,1,4,54,35.6,81,33,9967,0,0,0,0,0,0,0,0/
data (tab_auger(i,42),i=1,16) /12,2,1,1333,1124.6,7,341,6100,3559,0,0,0,0,0,0,0/
data (tab_auger(i,43),i=1,16) /12,2,2,104,20.1,108,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,44),i=1,16) /12,3,1,1360,1118.2,3,346,9654,0,0,0,0,0,0,0,0/
data (tab_auger(i,45),i=1,16) /12,4,1,1402,1107.5,12,386,9614,0,0,0,0,0,0,0,0/
data (tab_auger(i,46),i=1,16) /12,5,1,1444,1076.8,13,436,9564,0,0,0,0,0,0,0,0/
data (tab_auger(i,47),i=1,16) /12,6,1,1486,1036.1,27,497,9503,0,0,0,0,0,0,0,0/
data (tab_auger(i,48),i=1,16) /12,7,1,1528,980.4,16,569,9431,0,0,0,0,0,0,0,0/
data (tab_auger(i,49),i=1,16) /12,8,1,1570,948.5,33,424,9576,0,0,0,0,0,0,0,0/
data (tab_auger(i,50),i=1,16) /12,9,1,1612,934.5,11,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,51),i=1,16) /13,1,1,1564,1373.3,-7,5,391,6272,3332,0,0,0,0,0,0/
data (tab_auger(i,52),i=1,16) /13,1,2,121,85.1,142,0,272,9728,0,0,0,0,0,0,0/
data (tab_auger(i,53),i=1,16) /13,1,3,77.4,58.2,72,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,54),i=1,16) /13,1,4,77,54.4,31,44,9956,0,0,0,0,0,0,0,0/
data (tab_auger(i,55),i=1,16) /13,2,1,1590,1342.5,4,2,643,9355,0,0,0,0,0,0,0/
data (tab_auger(i,56),i=1,16) /13,2,2,134,21.7,66,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,57),i=1,16) /13,2,3,90.4,48.4,60,64,9936,0,0,0,0,0,0,0,0/
data (tab_auger(i,58),i=1,16) /13,2,4,90,48,61,64,9936,0,0,0,0,0,0,0,0/
data (tab_auger(i,59),i=1,16) /13,3,1,1618,1309.4,2,389,9611,0,0,0,0,0,0,0,0/
data (tab_auger(i,60),i=1,16) /13,3,2,148,11.4,79,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,61),i=1,16) /13,4,1,1646,1302.9,2,397,9603,0,0,0,0,0,0,0,0/
data (tab_auger(i,62),i=1,16) /13,5,1,1692,1288.2,11,442,9558,0,0,0,0,0,0,0,0/
data (tab_auger(i,63),i=1,16) /13,6,1,1738,1250.7,12,499,9501,0,0,0,0,0,0,0,0/
data (tab_auger(i,64),i=1,16) /13,7,1,1784,1202.7,25,568,9432,0,0,0,0,0,0,0,0/
data (tab_auger(i,65),i=1,16) /13,8,1,1830,1138.5,15,650,9350,0,0,0,0,0,0,0,0/
data (tab_auger(i,66),i=1,16) /13,9,1,1876,1106.4,30,484,9516,0,0,0,0,0,0,0,0/
data (tab_auger(i,67),i=1,16) /13,10,1,1922,1100.9,10,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,68),i=1,16) /14,1,1,1844,1652.4,10,12,438,548,8643,359,0,0,0,0,0/
data (tab_auger(i,69),i=1,16) /14,1,2,154,112.7,156,0,304,9696,0,0,0,0,0,0,0/
data (tab_auger(i,70),i=1,16) /14,1,3,104.7,82.8,63,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,71),i=1,16) /14,1,4,104,77.9,70,39,9961,0,0,0,0,0,0,0,0/
data (tab_auger(i,72),i=1,16) /14,2,1,1872,1591,7,7,442,6460,3091,0,0,0,0,0,0/
data (tab_auger(i,73),i=1,16) /14,2,2,169,93,113,0,333,9667,0,0,0,0,0,0,0/
data (tab_auger(i,74),i=1,16) /14,2,3,118.6,75.4,68,7,9993,0,0,0,0,0,0,0,0/
data (tab_auger(i,75),i=1,16) /14,2,4,118,69.8,122,59,9941,0,0,0,0,0,0,0,0/
data (tab_auger(i,76),i=1,16) /14,3,1,1901,1583.1,-1,3,775,9222,0,0,0,0,0,0,0/
data (tab_auger(i,77),i=1,16) /14,3,2,183,12.4,54,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,78),i=1,16) /14,3,3,133.6,61,48,81,9919,0,0,0,0,0,0,0,0/
data (tab_auger(i,79),i=1,16) /14,3,4,133,60.4,48,81,9919,0,0,0,0,0,0,0,0/
data (tab_auger(i,80),i=1,16) /14,4,1,1930,1510.9,1,442,9558,0,0,0,0,0,0,0,0/
data (tab_auger(i,81),i=1,16) /14,5,1,1959,1501.3,1,449,9551,0,0,0,0,0,0,0,0/
data (tab_auger(i,82),i=1,16) /14,6,1,2009,1482.2,11,507,9493,0,0,0,0,0,0,0,0/
data (tab_auger(i,83),i=1,16) /14,7,1,2059,1436.4,12,586,9414,0,0,0,0,0,0,0,0/
data (tab_auger(i,84),i=1,16) /14,8,1,2109,1375.9,23,699,9301,0,0,0,0,0,0,0,0/
data (tab_auger(i,85),i=1,16) /14,9,1,2160,1288.1,13,890,9110,0,0,0,0,0,0,0,0/
data (tab_auger(i,86),i=1,16) /14,10,1,2210,1249.4,23,768,9232,0,0,0,0,0,0,0,0/
data (tab_auger(i,87),i=1,16) /14,11,1,2260,1281.5,9,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,88),i=1,16) /15,1,1,2148,1920.8,14,23,603,692,6022,2660,0,0,0,0,0/
data (tab_auger(i,89),i=1,16) /15,1,2,191,143.9,142,1,307,9692,0,0,0,0,0,0,0/
data (tab_auger(i,90),i=1,16) /15,1,3,134.9,111.5,87,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,91),i=1,16) /15,1,4,134,109.6,78,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,92),i=1,16) /15,2,1,2178,1871.1,8,18,605,643,8395,339,0,0,0,0,0/
data (tab_auger(i,93),i=1,16) /15,2,2,207,123.2,129,1,311,9688,0,0,0,0,0,0,0/
data (tab_auger(i,94),i=1,16) /15,2,3,150.9,103.5,50,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,95),i=1,16) /15,2,4,150,96.5,56,56,9944,0,0,0,0,0,0,0,0/
data (tab_auger(i,96),i=1,16) /15,3,1,2208,1797.8,4,12,614,6960,2414,0,0,0,0,0,0/
data (tab_auger(i,97),i=1,16) /15,3,2,223,96.4,94,0,396,9604,0,0,0,0,0,0,0/
data (tab_auger(i,98),i=1,16) /15,3,3,166.9,92.3,56,11,9989,0,0,0,0,0,0,0,0/
data (tab_auger(i,99),i=1,16) /15,3,4,166,85,105,82,9918,0,0,0,0,0,0,0,0/
data (tab_auger(i,100),i=1,16) /15,4,1,2238,1782.8,-2,7,1019,8974,0,0,0,0,0,0,0/
data (tab_auger(i,101),i=1,16) /15,4,2,239,129.4,31,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,102),i=1,16) /15,4,3,183.9,74,40,108,9892,0,0,0,0,0,0,0,0/
data (tab_auger(i,103),i=1,16) /15,4,4,183,73.1,40,108,9892,0,0,0,0,0,0,0,0/
data (tab_auger(i,104),i=1,16) /15,5,1,2269,1698.1,0,623,9377,0,0,0,0,0,0,0,0/
data (tab_auger(i,105),i=1,16) /15,6,1,2300,1685.7,1,634,9366,0,0,0,0,0,0,0,0/
data (tab_auger(i,106),i=1,16) /15,7,1,2354,1658,10,719,9281,0,0,0,0,0,0,0,0/
data (tab_auger(i,107),i=1,16) /15,8,1,2408,1600,11,833,9167,0,0,0,0,0,0,0,0/
data (tab_auger(i,108),i=1,16) /15,9,1,2462,1523.5,20,998,9002,0,0,0,0,0,0,0,0/
data (tab_auger(i,109),i=1,16) /15,10,1,2517,1412.8,12,1279,8721,0,0,0,0,0,0,0,0/
data (tab_auger(i,110),i=1,16) /15,11,1,2571,1381.6,21,1102,8898,0,0,0,0,0,0,0,0/
data (tab_auger(i,111),i=1,16) /15,12,1,2625,1475.3,9,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,112),i=1,16) /16,1,1,2476,2172,12,46,812,831,5600,2385,326,0,0,0,0/
data (tab_auger(i,113),i=1,16) /16,1,2,232,180.1,122,1,336,9663,0,0,0,0,0,0,0/
data (tab_auger(i,114),i=1,16) /16,1,3,169.2,142.5,56,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,115),i=1,16) /16,1,4,168,140.8,56,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,116),i=1,16) /16,2,1,2508,2141.2,12,36,816,719,5697,2732,0,0,0,0,0/
data (tab_auger(i,117),i=1,16) /16,2,2,249,156.9,121,1,316,9683,0,0,0,0,0,0,0/
data (tab_auger(i,118),i=1,16) /16,2,3,186.2,135.7,76,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,119),i=1,16) /16,2,4,185,133.2,68,11,9989,0,0,0,0,0,0,0,0/
data (tab_auger(i,120),i=1,16) /16,3,1,2540,2083.8,6,28,824,708,8440,0,0,0,0,0,0/
data (tab_auger(i,121),i=1,16) /16,3,2,267,130.7,109,1,339,9660,0,0,0,0,0,0,0/
data (tab_auger(i,122),i=1,16) /16,3,3,204.2,124,41,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,123),i=1,16) /16,3,4,203,114.9,46,74,9926,0,0,0,0,0,0,0,0/
data (tab_auger(i,124),i=1,16) /16,4,1,2571,2014.6,3,19,840,9141,0,0,0,0,0,0,0/
data (tab_auger(i,125),i=1,16) /16,4,2,284,105.3,72,3,1465,8532,0,0,0,0,0,0,0/
data (tab_auger(i,126),i=1,16) /16,4,3,222.2,109.8,47,15,9985,0,0,0,0,0,0,0,0/
data (tab_auger(i,127),i=1,16) /16,4,4,221,100.7,91,107,9893,0,0,0,0,0,0,0,0/
data (tab_auger(i,128),i=1,16) /16,5,1,2603,1977.2,-3,12,1271,8717,0,0,0,0,0,0,0/
data (tab_auger(i,129),i=1,16) /16,5,2,302,149,27,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,130),i=1,16) /16,5,3,241.2,87.5,34,138,9862,0,0,0,0,0,0,0,0/
data (tab_auger(i,131),i=1,16) /16,5,4,240,86.3,34,138,9862,0,0,0,0,0,0,0,0/
data (tab_auger(i,132),i=1,16) /16,6,1,2636,1877.7,0,859,9141,0,0,0,0,0,0,0,0/
data (tab_auger(i,133),i=1,16) /16,7,1,2668,1860.3,1,875,9125,0,0,0,0,0,0,0,0/
data (tab_auger(i,134),i=1,16) /16,8,1,2726,1823.3,8,984,9016,0,0,0,0,0,0,0,0/
data (tab_auger(i,135),i=1,16) /16,9,1,2785,1754.4,10,1127,8873,0,0,0,0,0,0,0,0/
data (tab_auger(i,136),i=1,16) /16,10,1,2843,1664.1,19,1330,8670,0,0,0,0,0,0,0,0/
data (tab_auger(i,137),i=1,16) /16,11,1,2902,1534.9,11,1661,8339,0,0,0,0,0,0,0,0/
data (tab_auger(i,138),i=1,16) /16,12,1,2960,1511.9,19,1446,8554,0,0,0,0,0,0,0,0/
data (tab_auger(i,139),i=1,16) /16,13,1,3018,1683.1,8,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,140),i=1,16) /17,1,1,2829,2450.8,9,69,930,944,5569,2140,348,0,0,0,0/
data (tab_auger(i,141),i=1,16) /17,1,2,277,218.2,86,3,375,9622,0,0,0,0,0,0,0/
data (tab_auger(i,142),i=1,16) /17,1,3,207.6,175.6,39,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,143),i=1,16) /17,1,4,206,173.8,40,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,144),i=1,16) /17,2,1,2862,2422.2,10,57,936,821,5651,2200,335,0,0,0,0/
data (tab_auger(i,145),i=1,16) /17,2,2,295,195.5,104,2,354,9644,0,0,0,0,0,0,0/
data (tab_auger(i,146),i=1,16) /17,2,3,225.6,170.1,50,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,147),i=1,16) /17,2,4,224,167.8,49,7,9993,0,0,0,0,0,0,0,0/
data (tab_auger(i,148),i=1,16) /17,3,1,2896,2385.8,9,44,946,732,5912,2366,0,0,0,0,0/
data (tab_auger(i,149),i=1,16) /17,3,2,314,166.7,104,2,337,9661,0,0,0,0,0,0,0/
data (tab_auger(i,150),i=1,16) /17,3,3,244.6,159.8,65,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,151),i=1,16) /17,3,4,243,156.6,58,16,9984,0,0,0,0,0,0,0,0/
data (tab_auger(i,152),i=1,16) /17,4,1,2929,2367.6,3,36,959,1013,7992,0,0,0,0,0,0/
data (tab_auger(i,153),i=1,16) /17,4,2,333,140.1,84,6,1148,8846,0,0,0,0,0,0,0/
data (tab_auger(i,154),i=1,16) /17,4,3,264.7,145.1,33,6,9994,0,0,0,0,0,0,0,0/
data (tab_auger(i,155),i=1,16) /17,4,4,263,133.7,38,105,9895,0,0,0,0,0,0,0,0/
data (tab_auger(i,156),i=1,16) /17,5,1,2962,2235.6,2,25,979,8996,0,0,0,0,0,0,0/
data (tab_auger(i,157),i=1,16) /17,5,2,352,189.5,52,42,9958,0,0,0,0,0,0,0,0/
data (tab_auger(i,158),i=1,16) /17,5,3,284.7,128.1,40,21,9979,0,0,0,0,0,0,0,0/
data (tab_auger(i,159),i=1,16) /17,5,4,283,116.7,79,152,9848,0,0,0,0,0,0,0,0/
data (tab_auger(i,160),i=1,16) /17,6,1,2996,2195.4,-3,19,1417,8564,0,0,0,0,0,0,0/
data (tab_auger(i,161),i=1,16) /17,6,2,371,168.2,24,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,162),i=1,16) /17,6,3,304.7,100.5,29,196,9804,0,0,0,0,0,0,0,0/
data (tab_auger(i,163),i=1,16) /17,6,4,303,98.8,29,196,9804,0,0,0,0,0,0,0,0/
data (tab_auger(i,164),i=1,16) /17,7,1,3030,2080.3,0,1002,8998,0,0,0,0,0,0,0,0/
data (tab_auger(i,165),i=1,16) /17,8,1,3064,2059.6,0,1019,8981,0,0,0,0,0,0,0,0/
data (tab_auger(i,166),i=1,16) /17,9,1,3126,2013.3,8,1147,8853,0,0,0,0,0,0,0,0/
data (tab_auger(i,167),i=1,16) /17,10,1,3189,1931.3,9,1315,8685,0,0,0,0,0,0,0,0/
data (tab_auger(i,168),i=1,16) /17,11,1,3252,1825.6,17,1550,8450,0,0,0,0,0,0,0,0/
data (tab_auger(i,169),i=1,16) /17,12,1,3315,1674.9,10,1926,8074,0,0,0,0,0,0,0,0/
data (tab_auger(i,170),i=1,16) /17,13,1,3377,1664.4,18,1656,8344,0,0,0,0,0,0,0,0/
data (tab_auger(i,171),i=1,16) /17,14,1,3439,1905.2,7,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,172),i=1,16) /18,1,1,3206.3,2701.1,6,102,1172,1001,5172,2009,544,0,0,0,0/
data (tab_auger(i,173),i=1,16) /18,1,2,326.5,256.6,69,5,361,9634,0,0,0,0,0,0,0/
data (tab_auger(i,174),i=1,16) /18,1,3,250.7,212.7,35,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,175),i=1,16) /18,1,4,248.6,210.8,36,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,176),i=1,16) /18,2,1,3241,2673.1,8,87,1181,885,5243,2072,532,0,0,0,0/
data (tab_auger(i,177),i=1,16) /18,2,2,345,236.1,74,4,348,9648,0,0,0,0,0,0,0/
data (tab_auger(i,178),i=1,16) /18,2,3,269.1,205.7,34,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,179),i=1,16) /18,2,4,267,203.5,35,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,180),i=1,16) /18,3,1,3276,2641.7,9,72,1193,777,5423,2535,0,0,0,0,0/
data (tab_auger(i,181),i=1,16) /18,3,2,365,208.3,92,4,334,9662,0,0,0,0,0,0,0/
data (tab_auger(i,182),i=1,16) /18,3,3,289.1,197.5,43,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,183),i=1,16) /18,3,4,287,194.6,43,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,184),i=1,16) /18,4,1,3311,2624.3,7,56,1211,1005,7728,0,0,0,0,0,0/
data (tab_auger(i,185),i=1,16) /18,4,2,386,175.5,84,4,413,9583,0,0,0,0,0,0,0/
data (tab_auger(i,186),i=1,16) /18,4,3,310.1,184.6,57,6,9994,0,0,0,0,0,0,0,0/
data (tab_auger(i,187),i=1,16) /18,4,4,308,180.4,51,20,9980,0,0,0,0,0,0,0,0/
data (tab_auger(i,188),i=1,16) /18,5,1,3346,2573.4,2,47,1225,1092,7636,0,0,0,0,0,0/
data (tab_auger(i,189),i=1,16) /18,5,2,406,229.1,26,102,9898,0,0,0,0,0,0,0,0/
data (tab_auger(i,190),i=1,16) /18,5,3,332.2,166.8,28,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,191),i=1,16) /18,5,4,330,152.9,32,125,9875,0,0,0,0,0,0,0,0/
data (tab_auger(i,192),i=1,16) /18,6,1,3381,2428.6,2,35,1249,8716,0,0,0,0,0,0,0/
data (tab_auger(i,193),i=1,16) /18,6,2,426,211.5,47,77,9923,0,0,0,0,0,0,0,0/
data (tab_auger(i,194),i=1,16) /18,6,3,353.2,145.9,36,28,9972,0,0,0,0,0,0,0,0/
data (tab_auger(i,195),i=1,16) /18,6,4,351,132.5,70,178,9822,0,0,0,0,0,0,0,0/
data (tab_auger(i,196),i=1,16) /18,7,1,3417,2382.4,-3,29,1689,8282,0,0,0,0,0,0,0/
data (tab_auger(i,197),i=1,16) /18,7,2,447,188,22,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,198),i=1,16) /18,7,3,375.2,114.3,25,230,9770,0,0,0,0,0,0,0,0/
data (tab_auger(i,199),i=1,16) /18,7,4,373,112.2,25,230,9770,0,0,0,0,0,0,0,0/
data (tab_auger(i,200),i=1,16) /18,8,1,3452,2248.2,0,1284,8716,0,0,0,0,0,0,0,0/
data (tab_auger(i,201),i=1,16) /18,9,1,3488,2223.5,0,1305,8695,0,0,0,0,0,0,0,0/
data (tab_auger(i,202),i=1,16) /18,10,1,3554,2167.6,7,1447,8553,0,0,0,0,0,0,0,0/
data (tab_auger(i,203),i=1,16) /18,11,1,3621,2076.1,7,1623,8377,0,0,0,0,0,0,0,0/
data (tab_auger(i,204),i=1,16) /18,12,1,3688,1965.4,15,1843,8157,0,0,0,0,0,0,0,0/
data (tab_auger(i,205),i=1,16) /18,13,1,3755,1827.9,9,2116,7884,0,0,0,0,0,0,0,0/
data (tab_auger(i,206),i=1,16) /18,14,1,3821,1859.4,19,1671,8329,0,0,0,0,0,0,0,0/
data (tab_auger(i,207),i=1,16) /18,15,1,3887,2140.4,7,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,208),i=1,16) /19,1,1,3610,3057.1,14,110,861,969,3083,3005,1496,476,0,0,0/
data (tab_auger(i,209),i=1,16) /19,1,2,381,290.7,116,7,26,4364,5603,0,0,0,0,0,0/
data (tab_auger(i,210),i=1,16) /19,1,3,298.7,248.5,58,0,7843,2157,0,0,0,0,0,0,0/
data (tab_auger(i,211),i=1,16) /19,1,4,296,246,59,0,7813,2187,0,0,0,0,0,0,0/
data (tab_auger(i,212),i=1,16) /19,1,5,37,12.1,301,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,213),i=1,16) /19,2,1,3644,3040,5,111,1120,1170,5182,1930,487,0,0,0,0/
data (tab_auger(i,214),i=1,16) /19,2,2,399,273.8,57,7,383,9610,0,0,0,0,0,0,0/
data (tab_auger(i,215),i=1,16) /19,2,3,316.7,242.8,24,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,216),i=1,16) /19,2,4,314,240.3,25,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,217),i=1,16) /19,3,1,3681,3012,7,95,1126,1055,5450,2274,0,0,0,0,0/
data (tab_auger(i,218),i=1,16) /19,3,2,420,250.5,64,6,368,9626,0,0,0,0,0,0,0/
data (tab_auger(i,219),i=1,16) /19,3,3,337.7,236,29,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,220),i=1,16) /19,3,4,335,233.1,30,6,9994,0,0,0,0,0,0,0,0/
data (tab_auger(i,221),i=1,16) /19,4,1,3718,2991.5,6,78,1143,1306,7473,0,0,0,0,0,0/
data (tab_auger(i,222),i=1,16) /19,4,2,442,228.4,78,12,745,9243,0,0,0,0,0,0,0/
data (tab_auger(i,223),i=1,16) /19,4,3,359.8,225.7,39,6,9994,0,0,0,0,0,0,0,0/
data (tab_auger(i,224),i=1,16) /19,4,4,357,222,38,11,9989,0,0,0,0,0,0,0,0/
data (tab_auger(i,225),i=1,16) /19,5,1,3754,2943.7,6,61,1157,1222,7560,0,0,0,0,0,0/
data (tab_auger(i,226),i=1,16) /19,5,2,464,275.1,49,156,9844,0,0,0,0,0,0,0,0/
data (tab_auger(i,227),i=1,16) /19,5,3,382.8,210.1,51,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,228),i=1,16) /19,5,4,380,205,45,27,9973,0,0,0,0,0,0,0,0/
data (tab_auger(i,229),i=1,16) /19,6,1,3791,2885.2,2,54,1171,1342,7433,0,0,0,0,0,0/
data (tab_auger(i,230),i=1,16) /19,6,2,485,253.6,23,138,9862,0,0,0,0,0,0,0,0/
data (tab_auger(i,231),i=1,16) /19,6,3,405.8,188.5,25,11,9989,0,0,0,0,0,0,0,0/
data (tab_auger(i,232),i=1,16) /19,6,4,403,171.8,28,173,9827,0,0,0,0,0,0,0,0/
data (tab_auger(i,233),i=1,16) /19,7,1,3828,2723.4,1,42,1199,8759,0,0,0,0,0,0,0/
data (tab_auger(i,234),i=1,16) /19,7,2,507,234.2,43,106,9894,0,0,0,0,0,0,0,0/
data (tab_auger(i,235),i=1,16) /19,7,3,428.8,164.8,31,37,9963,0,0,0,0,0,0,0,0/
data (tab_auger(i,236),i=1,16) /19,7,4,426,148.4,63,250,9750,0,0,0,0,0,0,0,0/
data (tab_auger(i,237),i=1,16) /19,8,1,3866,2667.4,-3,40,1754,8206,0,0,0,0,0,0,0/
data (tab_auger(i,238),i=1,16) /19,8,2,530,208.6,20,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,239),i=1,16) /19,8,3,452.9,128.1,22,328,9672,0,0,0,0,0,0,0,0/
data (tab_auger(i,240),i=1,16) /19,8,4,450,125.4,23,328,9672,0,0,0,0,0,0,0,0/
data (tab_auger(i,241),i=1,16) /19,9,1,3902,2514.4,0,1228,8772,0,0,0,0,0,0,0,0/
data (tab_auger(i,242),i=1,16) /19,10,1,3940,2483.9,0,1253,8747,0,0,0,0,0,0,0,0/
data (tab_auger(i,243),i=1,16) /19,11,1,4010,2422.3,6,1393,8607,0,0,0,0,0,0,0,0/
data (tab_auger(i,244),i=1,16) /19,12,1,4081,2321.4,7,1568,8432,0,0,0,0,0,0,0,0/
data (tab_auger(i,245),i=1,16) /19,13,1,4152,2199.5,15,1787,8213,0,0,0,0,0,0,0,0/
data (tab_auger(i,246),i=1,16) /19,14,1,4223,2047.5,9,2064,7936,0,0,0,0,0,0,0,0/
data (tab_auger(i,247),i=1,16) /19,15,1,4293,2083.1,18,1626,8374,0,0,0,0,0,0,0,0/
data (tab_auger(i,248),i=1,16) /19,16,1,4363,2389,6,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,249),i=1,16) /20,1,1,4041,3348.2,21,1,167,1317,1294,4435,1713,785,288,0,0/
data (tab_auger(i,250),i=1,16) /20,1,2,441,339.9,155,0,40,523,8318,1119,0,0,0,0,0/
data (tab_auger(i,251),i=1,16) /20,1,3,352.6,299.2,69,0,227,9660,113,0,0,0,0,0,0/
data (tab_auger(i,252),i=1,16) /20,1,4,349,295.8,71,0,229,9658,113,0,0,0,0,0,0/
data (tab_auger(i,253),i=1,16) /20,1,5,46,11.7,206,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,254),i=1,16) /20,1,6,28.4,15.7,187,50,9950,0,0,0,0,0,0,0,0/
data (tab_auger(i,255),i=1,16) /20,1,7,28,15.3,189,50,9950,0,0,0,0,0,0,0,0/
data (tab_auger(i,256),i=1,16) /20,2,1,4075,3312.9,15,140,1011,1030,2966,3005,1455,393,0,0,0/
data (tab_auger(i,257),i=1,16) /20,2,2,459,310.9,118,10,32,4185,5773,0,0,0,0,0,0/
data (tab_auger(i,258),i=1,16) /20,2,3,369.6,281.9,57,0,7749,2251,0,0,0,0,0,0,0/
data (tab_auger(i,259),i=1,16) /20,2,4,366,278.6,58,0,7721,2279,0,0,0,0,0,0,0/
data (tab_auger(i,260),i=1,16) /20,2,5,58,7.1,209,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,261),i=1,16) /20,3,1,4110,3291.1,4,141,1327,1196,5147,1810,379,0,0,0,0/
data (tab_auger(i,262),i=1,16) /20,3,2,478,288.1,49,9,450,9541,0,0,0,0,0,0,0/
data (tab_auger(i,263),i=1,16) /20,3,3,389.6,274.2,20,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,264),i=1,16) /20,3,4,386,270.9,21,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,265),i=1,16) /20,4,1,4149,3257.3,6,121,1339,1089,5531,1920,0,0,0,0,0/
data (tab_auger(i,266),i=1,16) /20,4,2,501,272.3,60,17,835,9148,0,0,0,0,0,0,0/
data (tab_auger(i,267),i=1,16) /20,4,3,412.6,265.6,25,6,9994,0,0,0,0,0,0,0,0/
data (tab_auger(i,268),i=1,16) /20,4,4,409,261.9,26,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,269),i=1,16) /20,5,1,4187,3231.6,5,100,1362,1327,7211,0,0,0,0,0,0/
data (tab_auger(i,270),i=1,16) /20,5,2,525,241.7,65,27,1368,8605,0,0,0,0,0,0,0/
data (tab_auger(i,271),i=1,16) /20,5,3,436.6,253.5,35,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,272),i=1,16) /20,5,4,433,248.8,34,15,9985,0,0,0,0,0,0,0,0/
data (tab_auger(i,273),i=1,16) /20,6,1,4225,3177.6,6,79,1381,1255,7285,0,0,0,0,0,0/
data (tab_auger(i,274),i=1,16) /20,6,2,548,302.6,45,185,9815,0,0,0,0,0,0,0,0/
data (tab_auger(i,275),i=1,16) /20,6,3,461.6,235.3,47,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,276),i=1,16) /20,6,4,458,229,41,35,9965,0,0,0,0,0,0,0,0/
data (tab_auger(i,277),i=1,16) /20,7,1,4264,3110.7,2,73,1398,1429,7100,0,0,0,0,0,0/
data (tab_auger(i,278),i=1,16) /20,7,2,571,278.8,20,164,9836,0,0,0,0,0,0,0,0/
data (tab_auger(i,279),i=1,16) /20,7,3,486.7,210.9,22,14,9986,0,0,0,0,0,0,0,0/
data (tab_auger(i,280),i=1,16) /20,7,4,483,190.9,24,223,9777,0,0,0,0,0,0,0,0/
data (tab_auger(i,281),i=1,16) /20,8,1,4303,2935.2,1,59,1433,8508,0,0,0,0,0,0,0/
data (tab_auger(i,282),i=1,16) /20,8,2,595,257.7,40,124,9876,0,0,0,0,0,0,0,0/
data (tab_auger(i,283),i=1,16) /20,8,3,511.7,184.6,28,47,9953,0,0,0,0,0,0,0,0/
data (tab_auger(i,284),i=1,16) /20,8,4,508,164.9,56,326,9674,0,0,0,0,0,0,0,0/
data (tab_auger(i,285),i=1,16) /20,9,1,4342,2870.6,-3,63,2024,7913,0,0,0,0,0,0,0/
data (tab_auger(i,286),i=1,16) /20,9,2,620,230.2,18,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,287),i=1,16) /20,9,3,537.7,142.4,20,433,9567,0,0,0,0,0,0,0,0/
data (tab_auger(i,288),i=1,16) /20,9,4,534,138.9,20,433,9567,0,0,0,0,0,0,0,0/
data (tab_auger(i,289),i=1,16) /20,10,1,4380,2699.5,0,1475,8525,0,0,0,0,0,0,0,0/
data (tab_auger(i,290),i=1,16) /20,11,1,4420,2664.4,-1,1505,8495,0,0,0,0,0,0,0,0/
data (tab_auger(i,291),i=1,16) /20,12,1,4494,2589.3,6,1672,8328,0,0,0,0,0,0,0,0/
data (tab_auger(i,292),i=1,16) /20,13,1,4569,2470.8,6,1878,8122,0,0,0,0,0,0,0,0/
data (tab_auger(i,293),i=1,16) /20,14,1,4644,2327.4,13,2140,7860,0,0,0,0,0,0,0,0/
data (tab_auger(i,294),i=1,16) /20,15,1,4719,2148.3,8,2478,7522,0,0,0,0,0,0,0,0/
data (tab_auger(i,295),i=1,16) /20,16,1,4793,2210.2,16,1984,8016,0,0,0,0,0,0,0,0/
data (tab_auger(i,296),i=1,16) /20,17,1,4867,2652,6,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,297),i=1,16) /21,1,1,4494,3523.8,20,0,205,1426,1218,3457,2509,958,227,0,0/
data (tab_auger(i,298),i=1,16) /21,1,2,503,406.1,181,0,17,155,5263,4527,38,0,0,0,0/
data (tab_auger(i,299),i=1,16) /21,1,3,407.7,356.7,99,3,14,8078,1905,0,0,0,0,0,0/
data (tab_auger(i,300),i=1,16) /21,1,4,403,351.3,102,1,9,8001,1989,0,0,0,0,0,0/
data (tab_auger(i,301),i=1,16) /21,1,5,55,33.7,402,0,358,9642,0,0,0,0,0,0,0/
data (tab_auger(i,302),i=1,16) /21,1,6,33.5,18.7,135,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,303),i=1,16) /21,1,7,33,18.2,138,3,9997,0,0,0,0,0,0,0,0/
data (tab_auger(i,304),i=1,16) /21,2,1,4530,3503.9,17,0,350,1748,1754,4493,1364,291,0,0,0/
data (tab_auger(i,305),i=1,16) /21,2,2,523,369.8,141,0,47,1982,7971,0,0,0,0,0,0/
data (tab_auger(i,306),i=1,16) /21,2,3,427.7,344.2,81,2,1046,8952,0,0,0,0,0,0,0/
data (tab_auger(i,307),i=1,16) /21,2,4,423,339,83,1,750,9249,0,0,0,0,0,0,0/
data (tab_auger(i,308),i=1,16) /21,2,5,68,7.7,300,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,309),i=1,16) /21,2,6,46.5,18.9,214,4,9996,0,0,0,0,0,0,0,0/
data (tab_auger(i,310),i=1,16) /21,2,7,46,18.4,217,5,9995,0,0,0,0,0,0,0,0/
data (tab_auger(i,311),i=1,16) /21,3,1,4567,3470.3,9,205,1790,1225,4991,1587,202,0,0,0,0/
data (tab_auger(i,312),i=1,16) /21,3,2,543,336.2,84,13,510,9477,0,0,0,0,0,0,0/
data (tab_auger(i,313),i=1,16) /21,3,3,447.7,318.8,43,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,314),i=1,16) /21,3,4,443,313.1,45,7,9993,0,0,0,0,0,0,0,0/
data (tab_auger(i,315),i=1,16) /21,4,1,4604,3453.4,1,205,1802,1596,6397,0,0,0,0,0,0/
data (tab_auger(i,316),i=1,16) /21,4,2,564,314,44,30,1046,8924,0,0,0,0,0,0,0/
data (tab_auger(i,317),i=1,16) /21,4,3,469.7,306.7,17,7,9993,0,0,0,0,0,0,0,0/
data (tab_auger(i,318),i=1,16) /21,4,4,465,302.5,19,7,9993,0,0,0,0,0,0,0,0/
data (tab_auger(i,319),i=1,16) /21,5,1,4644,3411.4,4,176,1825,1458,6541,0,0,0,0,0,0/
data (tab_auger(i,320),i=1,16) /21,5,2,589,361.1,17,274,9726,0,0,0,0,0,0,0,0/
data (tab_auger(i,321),i=1,16) /21,5,3,494.7,296.3,22,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,322),i=1,16) /21,5,4,490,291.5,22,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,323),i=1,16) /21,6,1,4684,3356.2,5,146,1852,1321,6681,0,0,0,0,0,0/
data (tab_auger(i,324),i=1,16) /21,6,2,614,347.3,17,263,9737,0,0,0,0,0,0,0,0/
data (tab_auger(i,325),i=1,16) /21,6,3,520.7,282.4,31,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,326),i=1,16) /21,6,4,516,276.5,31,19,9981,0,0,0,0,0,0,0,0/
data (tab_auger(i,327),i=1,16) /21,7,1,4724,3295.7,5,116,1884,1248,6752,0,0,0,0,0,0/
data (tab_auger(i,328),i=1,16) /21,7,2,639,329.7,41,247,9753,0,0,0,0,0,0,0,0/
data (tab_auger(i,329),i=1,16) /21,7,3,547.7,261.6,42,13,9987,0,0,0,0,0,0,0,0/
data (tab_auger(i,330),i=1,16) /21,7,4,543,253.8,37,44,9956,0,0,0,0,0,0,0,0/
data (tab_auger(i,331),i=1,16) /21,8,1,4764,3218.3,2,112,1908,1478,6502,0,0,0,0,0,0/
data (tab_auger(i,332),i=1,16) /21,8,2,664,303.8,18,219,9781,0,0,0,0,0,0,0,0/
data (tab_auger(i,333),i=1,16) /21,8,3,574.7,234.4,19,18,9982,0,0,0,0,0,0,0,0/
data (tab_auger(i,334),i=1,16) /21,8,4,570,210.3,21,279,9721,0,0,0,0,0,0,0,0/
data (tab_auger(i,335),i=1,16) /21,9,1,4805,3035,1,95,1955,7950,0,0,0,0,0,0,0/
data (tab_auger(i,336),i=1,16) /21,9,2,690,281.4,37,167,9833,0,0,0,0,0,0,0,0/
data (tab_auger(i,337),i=1,16) /21,9,3,601.7,205.1,26,60,9940,0,0,0,0,0,0,0,0/
data (tab_auger(i,338),i=1,16) /21,9,4,597,181.5,51,404,9596,0,0,0,0,0,0,0,0/
data (tab_auger(i,339),i=1,16) /21,10,1,4846,2962.3,-3,107,2544,7349,0,0,0,0,0,0,0/
data (tab_auger(i,340),i=1,16) /21,10,2,717,252.4,17,7,9993,0,0,0,0,0,0,0,0/
data (tab_auger(i,341),i=1,16) /21,10,3,629.7,157.4,18,533,9467,0,0,0,0,0,0,0,0/
data (tab_auger(i,342),i=1,16) /21,10,4,625,153,18,533,9467,0,0,0,0,0,0,0,0/
data (tab_auger(i,343),i=1,16) /21,11,1,4886,2774.3,0,2036,7964,0,0,0,0,0,0,0,0/
data (tab_auger(i,344),i=1,16) /21,12,1,4928,2735.3,-1,2073,7927,0,0,0,0,0,0,0,0/
data (tab_auger(i,345),i=1,16) /21,13,1,5006,2633,5,2304,7696,0,0,0,0,0,0,0,0/
data (tab_auger(i,346),i=1,16) /21,14,1,5085,2479.8,6,2596,7404,0,0,0,0,0,0,0,0/
data (tab_auger(i,347),i=1,16) /21,15,1,5164,2289.6,11,2979,7021,0,0,0,0,0,0,0,0/
data (tab_auger(i,348),i=1,16) /21,16,1,5243,2039.2,6,3521,6479,0,0,0,0,0,0,0,0/
data (tab_auger(i,349),i=1,16) /21,17,1,5321,2140.1,13,2963,7037,0,0,0,0,0,0,0,0/
data (tab_auger(i,350),i=1,16) /21,18,1,5399,2930,6,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,351),i=1,16) /22,1,1,4970,3769,19,2,261,322,1895,2094,3959,1168,278,21,0/
data (tab_auger(i,352),i=1,16) /22,1,2,567,470.6,200,0,25,60,3288,6045,582,0,0,0,0/
data (tab_auger(i,353),i=1,16) /22,1,3,465,432.6,132,18,228,1853,7838,63,0,0,0,0,0/
data (tab_auger(i,354),i=1,16) /22,1,4,459,426,140,4,48,1485,8395,68,0,0,0,0,0/
data (tab_auger(i,355),i=1,16) /22,1,5,64,42.2,391,0,385,9615,0,0,0,0,0,0,0/
data (tab_auger(i,356),i=1,16) /22,1,6,38.7,23.2,127,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,357),i=1,16) /22,1,7,38,22.6,132,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,358),i=1,16) /22,2,1,5008,3740.2,16,2,260,2061,1303,4021,1820,498,35,0,0/
data (tab_auger(i,359),i=1,16) /22,2,2,589,426.3,131,0,22,859,8474,645,0,0,0,0,0/
data (tab_auger(i,360),i=1,16) /22,2,3,488,395.8,91,17,225,9651,107,0,0,0,0,0,0/
data (tab_auger(i,361),i=1,16) /22,2,4,482,387.1,92,4,45,9922,29,0,0,0,0,0,0/
data (tab_auger(i,362),i=1,16) /22,2,5,79,36,511,0,205,9795,0,0,0,0,0,0,0/
data (tab_auger(i,363),i=1,16) /22,2,6,53.7,24.1,213,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,364),i=1,16) /22,2,7,53,23.4,216,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,365),i=1,16) /22,3,1,5047,3714.8,9,2,2263,938,2030,4127,640,0,0,0,0/
data (tab_auger(i,366),i=1,16) /22,3,2,611,395.6,69,0,556,6673,2771,0,0,0,0,0,0/
data (tab_auger(i,367),i=1,16) /22,3,3,510,369,14,17,9870,113,0,0,0,0,0,0,0/
data (tab_auger(i,368),i=1,16) /22,3,4,504,360.3,13,4,9807,189,0,0,0,0,0,0,0/
data (tab_auger(i,369),i=1,16) /22,3,5,93,35.1,150,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,370),i=1,16) /22,3,6,68.7,11.9,203,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,371),i=1,16) /22,3,7,68,11.2,205,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,372),i=1,16) /22,4,1,5086,3680.3,7,243,2066,1264,5038,1389,0,0,0,0,0/
data (tab_auger(i,373),i=1,16) /22,4,2,633,368.8,68,36,1038,8926,0,0,0,0,0,0,0/
data (tab_auger(i,374),i=1,16) /22,4,3,533,353.6,40,17,9983,0,0,0,0,0,0,0,0/
data (tab_auger(i,375),i=1,16) /22,4,4,527,347.1,41,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,376),i=1,16) /22,5,1,5125,3653.9,1,242,2081,1593,6084,0,0,0,0,0,0/
data (tab_auger(i,377),i=1,16) /22,5,2,656,401.5,12,371,9629,0,0,0,0,0,0,0,0/
data (tab_auger(i,378),i=1,16) /22,5,3,556.9,339.1,14,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,379),i=1,16) /22,5,4,551,333.8,15,8,9992,0,0,0,0,0,0,0,0/
data (tab_auger(i,380),i=1,16) /22,6,1,5167,3606.9,3,208,2111,1451,6230,0,0,0,0,0,0/
data (tab_auger(i,381),i=1,16) /22,6,2,683,390.6,15,363,9637,0,0,0,0,0,0,0,0/
data (tab_auger(i,382),i=1,16) /22,6,3,583.9,327.8,19,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,383),i=1,16) /22,6,4,578,321.8,20,13,9987,0,0,0,0,0,0,0,0/
data (tab_auger(i,384),i=1,16) /22,7,1,5209,3545.4,4,173,2146,1321,6360,0,0,0,0,0,0/
data (tab_auger(i,385),i=1,16) /22,7,2,710,375.7,16,350,9650,0,0,0,0,0,0,0,0/
data (tab_auger(i,386),i=1,16) /22,7,3,611.9,312.2,29,12,9988,0,0,0,0,0,0,0,0/
data (tab_auger(i,387),i=1,16) /22,7,4,606,304.9,29,23,9977,0,0,0,0,0,0,0,0/
data (tab_auger(i,388),i=1,16) /22,8,1,5251,3477.9,5,139,2186,1266,6409,0,0,0,0,0,0/
data (tab_auger(i,389),i=1,16) /22,8,2,737,356.5,38,329,9671,0,0,0,0,0,0,0,0/
data (tab_auger(i,390),i=1,16) /22,8,3,640.9,288.9,39,16,9984,0,0,0,0,0,0,0,0/
data (tab_auger(i,391),i=1,16) /22,8,4,635,279.3,34,54,9946,0,0,0,0,0,0,0,0/
data (tab_auger(i,392),i=1,16) /22,9,1,5292,3390,1,141,2210,1559,6090,0,0,0,0,0,0/
data (tab_auger(i,393),i=1,16) /22,9,2,764,328.8,16,293,9707,0,0,0,0,0,0,0,0/
data (tab_auger(i,394),i=1,16) /22,9,3,669.9,258.8,17,22,9978,0,0,0,0,0,0,0,0/
data (tab_auger(i,395),i=1,16) /22,9,4,664,230.1,19,340,9660,0,0,0,0,0,0,0,0/
data (tab_auger(i,396),i=1,16) /22,10,1,5335,3196.5,1,125,2264,7611,0,0,0,0,0,0,0/
data (tab_auger(i,397),i=1,16) /22,10,2,792,305.2,34,224,9776,0,0,0,0,0,0,0,0/
data (tab_auger(i,398),i=1,16) /22,10,3,698.8,226.5,24,74,9926,0,0,0,0,0,0,0,0/
data (tab_auger(i,399),i=1,16) /22,10,4,693,198.4,47,488,9512,0,0,0,0,0,0,0,0/
data (tab_auger(i,400),i=1,16) /22,11,1,5378,3115.7,-3,149,2864,6987,0,0,0,0,0,0,0/
data (tab_auger(i,401),i=1,16) /22,11,2,821,275.2,16,9,9991,0,0,0,0,0,0,0,0/
data (tab_auger(i,402),i=1,16) /22,11,3,728.8,172.7,16,639,9361,0,0,0,0,0,0,0,0/
data (tab_auger(i,403),i=1,16) /22,11,4,723,167.3,17,639,9361,0,0,0,0,0,0,0,0/
data (tab_auger(i,404),i=1,16) /22,12,1,5420,2908.7,-1,2369,7631,0,0,0,0,0,0,0,0/
data (tab_auger(i,405),i=1,16) /22,13,1,5464,2865.1,-1,2411,7589,0,0,0,0,0,0,0,0/
data (tab_auger(i,406),i=1,16) /22,14,1,5546,2744.5,5,2670,7330,0,0,0,0,0,0,0,0/
data (tab_auger(i,407),i=1,16) /22,15,1,5629,2569,5,2994,7006,0,0,0,0,0,0,0,0/
data (tab_auger(i,408),i=1,16) /22,16,1,5712,2351.7,10,3418,6582,0,0,0,0,0,0,0,0/
data (tab_auger(i,409),i=1,16) /22,17,1,5795,2065.7,6,4016,5984,0,0,0,0,0,0,0,0/
data (tab_auger(i,410),i=1,16) /22,18,1,5877,2191.1,12,3438,6562,0,0,0,0,0,0,0,0/
data (tab_auger(i,411),i=1,16) /22,19,1,5959,3221,6,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,412),i=1,16) /23,1,1,5470,4002.9,19,6,340,488,1479,1675,4056,1595,331,30,0/
data (tab_auger(i,413),i=1,16) /23,1,2,633,549.2,215,0,32,93,834,6643,2394,4,0,0,0/
data (tab_auger(i,414),i=1,16) /23,1,3,525.5,492.5,116,50,576,2342,5697,1335,0,0,0,0,0/
data (tab_auger(i,415),i=1,16) /23,1,4,518,485.7,128,12,137,2031,6295,1525,0,0,0,0,0/
data (tab_auger(i,416),i=1,16) /23,1,5,72,50.8,392,0,454,9546,0,0,0,0,0,0,0/
data (tab_auger(i,417),i=1,16) /23,1,6,43.9,28,127,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,418),i=1,16) /23,1,7,43,27.2,132,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,419),i=1,16) /23,2,1,5511,3969.7,5,6,389,2573,3664,2485,760,117,6,0,0/
data (tab_auger(i,420),i=1,16) /23,2,2,658,488.3,63,0,29,4761,4499,711,0,0,0,0,0/
data (tab_auger(i,421),i=1,16) /23,2,3,550.5,458.6,50,59,1012,8618,311,0,0,0,0,0,0/
data (tab_auger(i,422),i=1,16) /23,2,4,543,446,52,14,239,9593,154,0,0,0,0,0,0/
data (tab_auger(i,423),i=1,16) /23,2,5,89,29.3,346,0,21,9979,0,0,0,0,0,0,0/
data (tab_auger(i,424),i=1,16) /23,2,6,60.9,29.7,210,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,425),i=1,16) /23,2,7,60,28.8,214,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,426),i=1,16) /23,3,1,5551,3933.5,7,5,2212,1222,3261,2542,683,75,0,0,0/
data (tab_auger(i,427),i=1,16) /23,3,2,682,461,43,0,128,8192,1680,0,0,0,0,0,0/
data (tab_auger(i,428),i=1,16) /23,3,3,575.4,427,21,50,8399,1551,0,0,0,0,0,0,0/
data (tab_auger(i,429),i=1,16) /23,3,4,568,415.3,21,11,8248,1741,0,0,0,0,0,0,0/
data (tab_auger(i,430),i=1,16) /23,3,5,106,45.2,148,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,431),i=1,16) /23,3,6,77.9,17.6,201,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,432),i=1,16) /23,3,7,77,16.7,203,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,433),i=1,16) /23,4,1,5592,3902.1,7,281,1922,1266,3650,2666,215,0,0,0,0/
data (tab_auger(i,434),i=1,16) /23,4,2,706,421.9,65,33,190,8758,1019,0,0,0,0,0,0/
data (tab_auger(i,435),i=1,16) /23,4,3,600.4,406.8,15,38,8321,1641,0,0,0,0,0,0,0/
data (tab_auger(i,436),i=1,16) /23,4,4,593,397.5,15,9,8195,1796,0,0,0,0,0,0,0/
data (tab_auger(i,437),i=1,16) /23,4,5,122,24.9,125,521,9479,0,0,0,0,0,0,0,0/
data (tab_auger(i,438),i=1,16) /23,5,1,5633,3865.3,6,282,2348,1226,4815,1329,0,0,0,0,0/
data (tab_auger(i,439),i=1,16) /23,5,2,730,419.3,37,118,2611,7271,0,0,0,0,0,0,0/
data (tab_auger(i,440),i=1,16) /23,5,3,625.4,388.5,36,30,9970,0,0,0,0,0,0,0,0/
data (tab_auger(i,441),i=1,16) /23,5,4,618,381.2,37,14,9986,0,0,0,0,0,0,0,0/
data (tab_auger(i,442),i=1,16) /23,6,1,5674,3837.7,1,280,2366,1555,5799,0,0,0,0,0,0/
data (tab_auger(i,443),i=1,16) /23,6,2,755,431.5,10,469,9531,0,0,0,0,0,0,0,0/
data (tab_auger(i,444),i=1,16) /23,6,3,651.4,372.5,11,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,445),i=1,16) /23,6,4,644,366,13,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,446),i=1,16) /23,7,1,5717,3784.9,3,241,2404,1421,5934,0,0,0,0,0,0/
data (tab_auger(i,447),i=1,16) /23,7,2,784,420,13,459,9541,0,0,0,0,0,0,0,0/
data (tab_auger(i,448),i=1,16) /23,7,3,680.3,360.3,17,12,9988,0,0,0,0,0,0,0,0/
data (tab_auger(i,449),i=1,16) /23,7,4,673,352.9,18,16,9984,0,0,0,0,0,0,0,0/
data (tab_auger(i,450),i=1,16) /23,8,1,5761,3717.7,4,201,2446,1301,6052,0,0,0,0,0,0/
data (tab_auger(i,451),i=1,16) /23,8,2,813,403.8,14,443,9557,0,0,0,0,0,0,0,0/
data (tab_auger(i,452),i=1,16) /23,8,3,710.3,342.9,26,15,9985,0,0,0,0,0,0,0,0/
data (tab_auger(i,453),i=1,16) /23,8,4,703,334,26,28,9972,0,0,0,0,0,0,0,0/
data (tab_auger(i,454),i=1,16) /23,9,1,5805,3643.8,5,163,2492,1271,6074,0,0,0,0,0,0/
data (tab_auger(i,455),i=1,16) /23,9,2,842,383.2,36,417,9583,0,0,0,0,0,0,0,0/
data (tab_auger(i,456),i=1,16) /23,9,3,741.3,317.1,36,19,9981,0,0,0,0,0,0,0,0/
data (tab_auger(i,457),i=1,16) /23,9,4,734,305.5,32,66,9934,0,0,0,0,0,0,0,0/
data (tab_auger(i,458),i=1,16) /23,10,1,5848,3546.7,1,174,2513,1628,5685,0,0,0,0,0,0/
data (tab_auger(i,459),i=1,16) /23,10,2,871,353.7,15,371,9629,0,0,0,0,0,0,0,0/
data (tab_auger(i,460),i=1,16) /23,10,3,772.3,284.1,16,27,9973,0,0,0,0,0,0,0,0/
data (tab_auger(i,461),i=1,16) /23,10,4,765,250,17,413,9587,0,0,0,0,0,0,0,0/
data (tab_auger(i,462),i=1,16) /23,11,1,5893,3344.6,1,163,2571,7266,0,0,0,0,0,0,0/
data (tab_auger(i,463),i=1,16) /23,11,2,901,329.2,32,285,9715,0,0,0,0,0,0,0,0/
data (tab_auger(i,464),i=1,16) /23,11,3,803.2,248.7,22,92,9908,0,0,0,0,0,0,0,0/
data (tab_auger(i,465),i=1,16) /23,11,4,796,215.4,44,592,9408,0,0,0,0,0,0,0,0/
data (tab_auger(i,466),i=1,16) /23,12,1,5938,3256.3,-2,206,3178,6616,0,0,0,0,0,0,0/
data (tab_auger(i,467),i=1,16) /23,12,2,932,298.8,14,11,9989,0,0,0,0,0,0,0,0/
data (tab_auger(i,468),i=1,16) /23,12,3,835.2,187.9,15,773,9227,0,0,0,0,0,0,0,0/
data (tab_auger(i,469),i=1,16) /23,12,4,828,181.3,15,773,9227,0,0,0,0,0,0,0,0/
data (tab_auger(i,470),i=1,16) /23,13,1,5982,3029.9,-1,2705,7295,0,0,0,0,0,0,0,0/
data (tab_auger(i,471),i=1,16) /23,14,1,6028,2981.9,-1,2751,7249,0,0,0,0,0,0,0,0/
data (tab_auger(i,472),i=1,16) /23,15,1,6114,2843.7,4,3031,6969,0,0,0,0,0,0,0,0/
data (tab_auger(i,473),i=1,16) /23,16,1,6201,2646.8,5,3378,6622,0,0,0,0,0,0,0,0/
data (tab_auger(i,474),i=1,16) /23,17,1,6288,2406.6,9,3826,6174,0,0,0,0,0,0,0,0/
data (tab_auger(i,475),i=1,16) /23,18,1,6375,2095.2,5,4442,5558,0,0,0,0,0,0,0,0/
data (tab_auger(i,476),i=1,16) /23,19,1,6461,2248.9,10,3838,6162,0,0,0,0,0,0,0,0/
data (tab_auger(i,477),i=1,16) /23,20,1,6547,3526,5,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,478),i=1,16) /24,1,1,5995,4201.6,9,17,607,723,1364,1118,3777,2070,304,23,0/
data (tab_auger(i,479),i=1,16) /24,1,2,702,603.2,164,0,66,400,1439,4481,3564,50,0,0,0/
data (tab_auger(i,480),i=1,16) /24,1,3,589.2,545.9,79,94,1557,2791,4486,1064,8,0,0,0,0/
data (tab_auger(i,481),i=1,16) /24,1,4,580,536.8,86,52,979,2809,4926,1225,9,0,0,0,0/
data (tab_auger(i,482),i=1,16) /24,1,5,80,56.9,376,0,352,9648,0,0,0,0,0,0,0/
data (tab_auger(i,483),i=1,16) /24,1,6,49.2,32,124,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,484),i=1,16) /24,1,7,48,30.8,127,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,485),i=1,16) /24,2,1,6038,4194.3,7,17,608,743,2738,3759,1860,253,24,0,0/
data (tab_auger(i,486),i=1,16) /24,2,2,730,561.8,98,0,64,502,6549,2885,0,0,0,0,0/
data (tab_auger(i,487),i=1,16) /24,2,3,616.2,526.6,69,94,1561,2860,5485,0,0,0,0,0,0/
data (tab_auger(i,488),i=1,16) /24,2,4,607,514.1,73,52,978,2894,6076,0,0,0,0,0,0/
data (tab_auger(i,489),i=1,16) /24,2,5,99,36.7,340,0,135,9865,0,0,0,0,0,0,0/
data (tab_auger(i,490),i=1,16) /24,2,6,68.2,34.7,205,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,491),i=1,16) /24,2,7,67,33.5,208,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,492),i=1,16) /24,3,1,6080,4155.1,2,13,441,2540,5392,1335,225,53,0,0,0/
data (tab_auger(i,493),i=1,16) /24,3,2,756,516.8,52,0,53,6270,3677,0,0,0,0,0,0/
data (tab_auger(i,494),i=1,16) /24,3,3,643.1,486.3,46,106,994,7653,1247,0,0,0,0,0,0/
data (tab_auger(i,495),i=1,16) /24,3,4,634,471.8,48,25,269,8316,1390,0,0,0,0,0,0/
data (tab_auger(i,496),i=1,16) /24,3,5,117,52,132,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,497),i=1,16) /24,3,6,87.2,22.5,177,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,498),i=1,16) /24,3,7,86,21.3,179,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,499),i=1,16) /24,4,1,6122,4124.7,5,11,2474,1621,4398,1100,396,0,0,0,0/
data (tab_auger(i,500),i=1,16) /24,4,2,782,482.9,34,0,158,9842,0,0,0,0,0,0,0/
data (tab_auger(i,501),i=1,16) /24,4,3,671.1,464.5,18,90,8397,1513,0,0,0,0,0,0,0/
data (tab_auger(i,502),i=1,16) /24,4,4,662,452.6,18,20,8286,1694,0,0,0,0,0,0,0/
data (tab_auger(i,503),i=1,16) /24,4,5,136,35.1,133,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,504),i=1,16) /24,4,6,105.2,4.9,172,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,505),i=1,16) /24,4,7,104,3.7,174,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,506),i=1,16) /24,5,1,6165,4087.5,7,321,2609,1200,4025,1845,0,0,0,0,0/
data (tab_auger(i,507),i=1,16) /24,5,2,808,468.7,50,86,1625,8289,0,0,0,0,0,0,0/
data (tab_auger(i,508),i=1,16) /24,5,3,698.1,443,10,69,9931,0,0,0,0,0,0,0,0/
data (tab_auger(i,509),i=1,16) /24,5,4,689,433.5,9,15,9985,0,0,0,0,0,0,0,0/
data (tab_auger(i,510),i=1,16) /24,5,5,154,11.6,112,565,9435,0,0,0,0,0,0,0,0/
data (tab_auger(i,511),i=1,16) /24,6,1,6208,4060.7,6,320,2629,1526,5525,0,0,0,0,0,0/
data (tab_auger(i,512),i=1,16) /24,6,2,834,437.7,35,149,2594,7257,0,0,0,0,0,0,0/
data (tab_auger(i,513),i=1,16) /24,6,3,725,424.2,32,49,9951,0,0,0,0,0,0,0,0/
data (tab_auger(i,514),i=1,16) /24,6,4,716,416.3,34,19,9981,0,0,0,0,0,0,0,0/
data (tab_auger(i,515),i=1,16) /24,7,1,6250,4017.4,1,316,2634,1528,5522,0,0,0,0,0,0/
data (tab_auger(i,516),i=1,16) /24,7,2,861,460.2,8,592,9408,0,0,0,0,0,0,0,0/
data (tab_auger(i,517),i=1,16) /24,7,3,753,407,9,12,9988,0,0,0,0,0,0,0,0/
data (tab_auger(i,518),i=1,16) /24,7,4,744,399,11,12,9988,0,0,0,0,0,0,0,0/
data (tab_auger(i,519),i=1,16) /24,8,1,6295,3959.8,3,272,2679,1400,5649,0,0,0,0,0,0/
data (tab_auger(i,520),i=1,16) /24,8,2,892,447.9,12,578,9422,0,0,0,0,0,0,0,0/
data (tab_auger(i,521),i=1,16) /24,8,3,784,393.5,15,14,9986,0,0,0,0,0,0,0,0/
data (tab_auger(i,522),i=1,16) /24,8,4,775,384.6,16,19,9981,0,0,0,0,0,0,0,0/
data (tab_auger(i,523),i=1,16) /24,9,1,6341,3887.2,4,227,2727,1291,5755,0,0,0,0,0,0/
data (tab_auger(i,524),i=1,16) /24,9,2,923,430.9,13,558,9442,0,0,0,0,0,0,0,0/
data (tab_auger(i,525),i=1,16) /24,9,3,815.9,374.6,25,18,9982,0,0,0,0,0,0,0,0/
data (tab_auger(i,526),i=1,16) /24,9,4,807,363.8,24,34,9966,0,0,0,0,0,0,0,0/
data (tab_auger(i,527),i=1,16) /24,10,1,6387,3806.8,4,186,2779,1283,5752,0,0,0,0,0,0/
data (tab_auger(i,528),i=1,16) /24,10,2,954,409,33,525,9475,0,0,0,0,0,0,0,0/
data (tab_auger(i,529),i=1,16) /24,10,3,848.9,346.2,33,23,9977,0,0,0,0,0,0,0,0/
data (tab_auger(i,530),i=1,16) /24,10,4,840,332.2,29,79,9921,0,0,0,0,0,0,0,0/
data (tab_auger(i,531),i=1,16) /24,11,1,6432,3700,1,211,2793,1696,5300,0,0,0,0,0,0/
data (tab_auger(i,532),i=1,16) /24,11,2,985,378.3,14,468,9532,0,0,0,0,0,0,0,0/
data (tab_auger(i,533),i=1,16) /24,11,3,881.8,310.6,15,33,9967,0,0,0,0,0,0,0,0/
data (tab_auger(i,534),i=1,16) /24,11,4,873,270.5,16,492,9508,0,0,0,0,0,0,0,0/
data (tab_auger(i,535),i=1,16) /24,12,1,6479,3489.9,1,206,2853,6941,0,0,0,0,0,0,0/
data (tab_auger(i,536),i=1,16) /24,12,2,1017,352.9,30,358,9642,0,0,0,0,0,0,0,0/
data (tab_auger(i,537),i=1,16) /24,12,3,914.8,271.6,21,112,9888,0,0,0,0,0,0,0,0/
data (tab_auger(i,538),i=1,16) /24,12,4,906,232.4,40,703,9297,0,0,0,0,0,0,0,0/
data (tab_auger(i,539),i=1,16) /24,13,1,6526,3393.9,-2,272,3463,6265,0,0,0,0,0,0,0/
data (tab_auger(i,540),i=1,16) /24,13,2,1050,323.1,14,14,9986,0,0,0,0,0,0,0,0/
data (tab_auger(i,541),i=1,16) /24,13,3,948.8,203.2,14,917,9083,0,0,0,0,0,0,0,0/
data (tab_auger(i,542),i=1,16) /24,13,4,940,195.3,14,917,9083,0,0,0,0,0,0,0,0/
data (tab_auger(i,543),i=1,16) /24,14,1,6572,3147.8,-1,3019,6981,0,0,0,0,0,0,0,0/
data (tab_auger(i,544),i=1,16) /24,15,1,6621,3096.1,-1,3068,6932,0,0,0,0,0,0,0,0/
data (tab_auger(i,545),i=1,16) /24,16,1,6711,2939.5,4,3365,6635,0,0,0,0,0,0,0,0/
data (tab_auger(i,546),i=1,16) /24,17,1,6802,2723,4,3731,6269,0,0,0,0,0,0,0,0/
data (tab_auger(i,547),i=1,16) /24,18,1,6893,2459.8,8,4195,5805,0,0,0,0,0,0,0,0/
data (tab_auger(i,548),i=1,16) /24,19,1,6984,2122.5,4,4826,5174,0,0,0,0,0,0,0,0/
data (tab_auger(i,549),i=1,16) /24,20,1,7074,2300.6,9,4214,5786,0,0,0,0,0,0,0,0/
data (tab_auger(i,550),i=1,16) /24,21,1,7164,3845,5,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,551),i=1,16) /25,1,1,6544,4425.9,13,27,706,838,1437,809,1538,3121,1252,262,10/
data (tab_auger(i,552),i=1,16) /25,1,2,775,687,189,0,94,473,1555,3400,3019,1458,1,0,0/
data (tab_auger(i,553),i=1,16) /25,1,3,656,612.8,76,133,1805,2968,4061,983,50,0,0,0,0/
data (tab_auger(i,554),i=1,16) /25,1,4,645,603.1,88,78,1026,2906,4775,1157,58,0,0,0,0/
data (tab_auger(i,555),i=1,16) /25,1,5,89,67.1,370,0,590,9410,0,0,0,0,0,0,0/
data (tab_auger(i,556),i=1,16) /25,1,6,54.5,36.6,96,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,557),i=1,16) /25,1,7,53,35.1,98,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,558),i=1,16) /25,2,1,6589,4416.2,13,27,707,829,1545,1738,3528,1419,206,1,0/
data (tab_auger(i,559),i=1,16) /25,2,2,804,634,167,0,90,462,1951,5257,2205,35,0,0,0/
data (tab_auger(i,560),i=1,16) /25,2,3,685,592,77,134,1809,2943,4278,835,1,0,0,0,0/
data (tab_auger(i,561),i=1,16) /25,2,4,674,578,87,79,1026,2876,5027,991,1,0,0,0,0/
data (tab_auger(i,562),i=1,16) /25,2,5,109,49.4,358,0,371,9629,0,0,0,0,0,0,0/
data (tab_auger(i,563),i=1,16) /25,2,6,75.5,41.1,198,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,564),i=1,16) /25,2,7,74,39.6,202,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,565),i=1,16) /25,3,1,6633,4376.1,6,27,708,1952,2204,3940,1023,146,0,0,0/
data (tab_auger(i,566),i=1,16) /25,3,2,833,593.1,85,0,114,2355,6756,775,0,0,0,0,0/
data (tab_auger(i,567),i=1,16) /25,3,3,713.9,563.6,58,134,1814,5669,2383,0,0,0,0,0,0/
data (tab_auger(i,568),i=1,16) /25,3,4,703,545.3,58,79,1025,6992,1904,0,0,0,0,0,0/
data (tab_auger(i,569),i=1,16) /25,3,5,129,14,321,0,128,9872,0,0,0,0,0,0,0/
data (tab_auger(i,570),i=1,16) /25,3,6,95.5,31.2,215,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,571),i=1,16) /25,3,7,94,29.6,219,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,572),i=1,16) /25,4,1,6677,4342.7,1,21,2755,1612,4159,1237,216,0,0,0,0/
data (tab_auger(i,573),i=1,16) /25,4,2,861,547.1,37,0,172,9828,0,0,0,0,0,0,0/
data (tab_auger(i,574),i=1,16) /25,4,3,743.9,528.7,19,152,8542,1306,0,0,0,0,0,0,0/
data (tab_auger(i,575),i=1,16) /25,4,4,733,510.2,19,37,8336,1627,0,0,0,0,0,0,0/
data (tab_auger(i,576),i=1,16) /25,4,5,149,42.6,118,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,577),i=1,16) /25,4,6,116.5,10.5,151,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,578),i=1,16) /25,4,7,115,9,153,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,579),i=1,16) /25,5,1,6721,4300.2,4,366,2867,1172,4390,1205,0,0,0,0,0/
data (tab_auger(i,580),i=1,16) /25,5,2,889,526.2,33,72,1234,8694,0,0,0,0,0,0,0/
data (tab_auger(i,581),i=1,16) /25,5,3,773.9,505.3,12,132,9868,0,0,0,0,0,0,0,0/
data (tab_auger(i,582),i=1,16) /25,5,4,763,489.9,11,31,9969,0,0,0,0,0,0,0,0/
data (tab_auger(i,583),i=1,16) /25,5,5,169,21.1,115,242,9758,0,0,0,0,0,0,0,0/
data (tab_auger(i,584),i=1,16) /25,6,1,6766,4268.9,5,365,2889,1509,5236,0,0,0,0,0,0/
data (tab_auger(i,585),i=1,16) /25,6,2,917,486.6,47,103,1590,8307,0,0,0,0,0,0,0/
data (tab_auger(i,586),i=1,16) /25,6,3,802.8,482.7,8,113,9887,0,0,0,0,0,0,0,0/
data (tab_auger(i,587),i=1,16) /25,6,4,792,470.9,7,34,9966,0,0,0,0,0,0,0,0/
data (tab_auger(i,588),i=1,16) /25,7,1,6810,4227.8,6,360,2896,1509,5235,0,0,0,0,0,0/
data (tab_auger(i,589),i=1,16) /25,7,2,945,511,25,673,9327,0,0,0,0,0,0,0,0/
data (tab_auger(i,590),i=1,16) /25,7,3,831.8,461.2,30,71,9929,0,0,0,0,0,0,0,0/
data (tab_auger(i,591),i=1,16) /25,7,4,821,451.8,31,26,9974,0,0,0,0,0,0,0,0/
data (tab_auger(i,592),i=1,16) /25,8,1,6854,4181.2,0,353,2904,1500,5243,0,0,0,0,0,0/
data (tab_auger(i,593),i=1,16) /25,8,2,974,487.7,7,737,9263,0,0,0,0,0,0,0,0/
data (tab_auger(i,594),i=1,16) /25,8,3,861.8,442.1,8,15,9985,0,0,0,0,0,0,0,0/
data (tab_auger(i,595),i=1,16) /25,8,4,851,432.6,10,15,9985,0,0,0,0,0,0,0,0/
data (tab_auger(i,596),i=1,16) /25,9,1,6901,4119.2,2,305,2955,1380,5360,0,0,0,0,0,0/
data (tab_auger(i,597),i=1,16) /25,9,2,1007,474.7,10,720,9280,0,0,0,0,0,0,0,0/
data (tab_auger(i,598),i=1,16) /25,9,3,894.7,427.6,14,18,9982,0,0,0,0,0,0,0,0/
data (tab_auger(i,599),i=1,16) /25,9,4,884,417,15,23,9977,0,0,0,0,0,0,0,0/
data (tab_auger(i,600),i=1,16) /25,10,1,6949,4041.3,3,255,3010,1281,5454,0,0,0,0,0,0/
data (tab_auger(i,601),i=1,16) /25,10,2,1040,457,12,695,9305,0,0,0,0,0,0,0,0/
data (tab_auger(i,602),i=1,16) /25,10,3,928.7,407.1,23,22,9978,0,0,0,0,0,0,0,0/
data (tab_auger(i,603),i=1,16) /25,10,4,918,394.2,23,42,9958,0,0,0,0,0,0,0,0/
data (tab_auger(i,604),i=1,16) /25,11,1,6997,3954.5,4,211,3068,1298,5423,0,0,0,0,0,0/
data (tab_auger(i,605),i=1,16) /25,11,2,1073,434.2,31,655,9345,0,0,0,0,0,0,0,0/
data (tab_auger(i,606),i=1,16) /25,11,3,963.7,376.3,31,28,9972,0,0,0,0,0,0,0,0/
data (tab_auger(i,607),i=1,16) /25,11,4,953,359.6,28,98,9902,0,0,0,0,0,0,0,0/
data (tab_auger(i,608),i=1,16) /25,12,1,7045,3838.4,1,256,3070,1768,4906,0,0,0,0,0,0/
data (tab_auger(i,609),i=1,16) /25,12,2,1106,402,13,585,9415,0,0,0,0,0,0,0,0/
data (tab_auger(i,610),i=1,16) /25,12,3,998.6,337.5,14,40,9960,0,0,0,0,0,0,0,0/
data (tab_auger(i,611),i=1,16) /25,12,4,988,290.7,14,591,9409,0,0,0,0,0,0,0,0/
data (tab_auger(i,612),i=1,16) /25,13,1,7094,3622,1,261,3129,6610,0,0,0,0,0,0,0/
data (tab_auger(i,613),i=1,16) /25,13,2,1140,376.5,28,448,9552,0,0,0,0,0,0,0,0/
data (tab_auger(i,614),i=1,16) /25,13,3,1033.6,295.3,19,137,9863,0,0,0,0,0,0,0,0/
data (tab_auger(i,615),i=1,16) /25,13,4,1023,248.9,37,847,9153,0,0,0,0,0,0,0,0/
data (tab_auger(i,616),i=1,16) /25,14,1,7143,3516.8,-2,365,3741,5894,0,0,0,0,0,0,0/
data (tab_auger(i,617),i=1,16) /25,14,2,1175,348.1,13,17,9983,0,0,0,0,0,0,0,0/
data (tab_auger(i,618),i=1,16) /25,14,3,1069.6,217.5,13,1113,8887,0,0,0,0,0,0,0,0/
data (tab_auger(i,619),i=1,16) /25,14,4,1059,208.1,13,1113,8887,0,0,0,0,0,0,0,0/
data (tab_auger(i,620),i=1,16) /25,15,1,7191,3252.4,-1,3335,6665,0,0,0,0,0,0,0,0/
data (tab_auger(i,621),i=1,16) /25,16,1,7242,3196.4,-1,3386,6614,0,0,0,0,0,0,0,0/
data (tab_auger(i,622),i=1,16) /25,17,1,7336,3022.6,3,3698,6302,0,0,0,0,0,0,0,0/
data (tab_auger(i,623),i=1,16) /25,18,1,7431,2785.8,4,4077,5923,0,0,0,0,0,0,0,0/
data (tab_auger(i,624),i=1,16) /25,19,1,7526,2501.4,7,4552,5448,0,0,0,0,0,0,0,0/
data (tab_auger(i,625),i=1,16) /25,20,1,7621,2142.6,4,5184,4816,0,0,0,0,0,0,0,0/
data (tab_auger(i,626),i=1,16) /25,21,1,7715,2346.7,8,4562,5438,0,0,0,0,0,0,0,0/
data (tab_auger(i,627),i=1,16) /25,22,1,7809,4179,5,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,628),i=1,16) /26,1,1,7117,4661.6,11,47,1005,984,1438,970,1598,2007,1729,205,22/
data (tab_auger(i,629),i=1,16) /26,1,2,851,769.3,166,0,162,893,2000,3317,2495,1122,11,0,0/
data (tab_auger(i,630),i=1,16) /26,1,3,726,683.4,61,158,2601,3096,3304,802,39,0,0,0,0/
data (tab_auger(i,631),i=1,16) /26,1,4,713,668.7,71,152,1768,3086,3981,966,47,0,0,0,0/
data (tab_auger(i,632),i=1,16) /26,1,5,98,76.8,322,0,641,9359,0,0,0,0,0,0,0/
data (tab_auger(i,633),i=1,16) /26,1,6,60.9,43,101,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,634),i=1,16) /26,1,7,59,41.1,104,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,635),i=1,16) /26,2,1,7164,4657.2,13,47,1007,976,1564,1859,3242,1154,143,7,0/
data (tab_auger(i,636),i=1,16) /26,2,2,882,730.4,154,0,170,953,2581,4829,1409,58,0,0,0/
data (tab_auger(i,637),i=1,16) /26,2,3,757,670.1,66,158,2606,3076,3542,617,1,0,0,0,0/
data (tab_auger(i,638),i=1,16) /26,2,4,744,651.1,76,153,1770,3061,4256,759,1,0,0,0,0/
data (tab_auger(i,639),i=1,16) /26,2,5,119,61.6,352,0,437,9563,0,0,0,0,0,0,0/
data (tab_auger(i,640),i=1,16) /26,2,6,82.9,48.3,148,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,641),i=1,16) /26,2,7,81,46.4,151,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,642),i=1,16) /26,3,1,7210,4614.9,7,47,1009,2220,2139,3601,858,123,3,0,0/
data (tab_auger(i,643),i=1,16) /26,3,2,913,673.1,84,0,194,3687,4976,1139,4,0,0,0,0/
data (tab_auger(i,644),i=1,16) /26,3,3,787.9,641.6,50,159,2612,5912,1317,0,0,0,0,0,0/
data (tab_auger(i,645),i=1,16) /26,3,4,775,618.1,55,153,1773,7056,1018,0,0,0,0,0,0/
data (tab_auger(i,646),i=1,16) /26,3,5,141,24.6,301,0,217,9783,0,0,0,0,0,0,0/
data (tab_auger(i,647),i=1,16) /26,3,6,104.9,40.8,203,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,648),i=1,16) /26,3,7,103,38.8,207,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,649),i=1,16) /26,4,1,7256,4569.7,6,41,3058,1243,3313,1850,398,97,0,0,0/
data (tab_auger(i,650),i=1,16) /26,4,2,943,629.8,70,0,275,8295,1426,4,0,0,0,0,0/
data (tab_auger(i,651),i=1,16) /26,4,3,819.9,604.7,33,180,8721,1045,54,0,0,0,0,0,0/
data (tab_auger(i,652),i=1,16) /26,4,4,807,579.9,40,113,8487,1333,67,0,0,0,0,0,0/
data (tab_auger(i,653),i=1,16) /26,4,5,162,55.8,146,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,654),i=1,16) /26,4,6,126.9,20.8,187,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,655),i=1,16) /26,4,7,125,18.9,190,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,656),i=1,16) /26,5,1,7301,4509.2,1,415,2624,1625,3961,1180,195,0,0,0,0/
data (tab_auger(i,657),i=1,16) /26,5,2,973,589.4,35,65,194,9741,0,0,0,0,0,0,0/
data (tab_auger(i,658),i=1,16) /26,5,3,851.9,573,17,208,8562,1230,0,0,0,0,0,0,0/
data (tab_auger(i,659),i=1,16) /26,5,4,839,549.2,16,54,8313,1633,0,0,0,0,0,0,0/
data (tab_auger(i,660),i=1,16) /26,5,5,184,29,104,177,9823,0,0,0,0,0,0,0,0/
data (tab_auger(i,661),i=1,16) /26,6,1,7348,4465,3,411,3117,1167,4220,1085,0,0,0,0,0/
data (tab_auger(i,662),i=1,16) /26,6,2,1003,545.8,31,85,1225,8690,0,0,0,0,0,0,0/
data (tab_auger(i,663),i=1,16) /26,6,3,883.8,548.9,10,183,9817,0,0,0,0,0,0,0,0/
data (tab_auger(i,664),i=1,16) /26,6,4,871,529.9,8,44,9956,0,0,0,0,0,0,0,0/
data (tab_auger(i,665),i=1,16) /26,6,5,205,3.9,108,347,9653,0,0,0,0,0,0,0,0/
data (tab_auger(i,666),i=1,16) /26,7,1,7394,4430.7,5,409,3143,1492,4956,0,0,0,0,0,0/
data (tab_auger(i,667),i=1,16) /26,7,2,1033,506.1,41,162,2048,7790,0,0,0,0,0,0,0/
data (tab_auger(i,668),i=1,16) /26,7,3,914.8,523.1,8,156,9844,0,0,0,0,0,0,0,0/
data (tab_auger(i,669),i=1,16) /26,7,4,902,509,7,46,9954,0,0,0,0,0,0,0,0/
data (tab_auger(i,670),i=1,16) /26,8,1,7440,4386.1,5,401,3153,1487,4959,0,0,0,0,0,0/
data (tab_auger(i,671),i=1,16) /26,8,2,1063,540.3,22,813,9187,0,0,0,0,0,0,0,0/
data (tab_auger(i,672),i=1,16) /26,8,3,945.8,498.3,27,99,9901,0,0,0,0,0,0,0,0/
data (tab_auger(i,673),i=1,16) /26,8,4,933,487.4,28,33,9967,0,0,0,0,0,0,0,0/
data (tab_auger(i,674),i=1,16) /26,9,1,7486,4336.8,0,389,3165,1471,4975,0,0,0,0,0,0/
data (tab_auger(i,675),i=1,16) /26,9,2,1094,514,6,895,9105,0,0,0,0,0,0,0,0/
data (tab_auger(i,676),i=1,16) /26,9,3,977.8,478.3,7,17,9983,0,0,0,0,0,0,0,0/
data (tab_auger(i,677),i=1,16) /26,9,4,965,467.2,8,17,9983,0,0,0,0,0,0,0,0/
data (tab_auger(i,678),i=1,16) /26,10,1,7535,4270.1,2,336,3222,1357,5085,0,0,0,0,0,0/
data (tab_auger(i,679),i=1,16) /26,10,2,1129,500.4,9,876,9124,0,0,0,0,0,0,0,0/
data (tab_auger(i,680),i=1,16) /26,10,3,1012.7,462.7,12,20,9980,0,0,0,0,0,0,0,0/
data (tab_auger(i,681),i=1,16) /26,10,4,1000,450.1,13,27,9973,0,0,0,0,0,0,0,0/
data (tab_auger(i,682),i=1,16) /26,11,1,7585,4187,3,282,3283,1267,5168,0,0,0,0,0,0/
data (tab_auger(i,683),i=1,16) /26,11,2,1164,482.4,11,846,9154,0,0,0,0,0,0,0,0/
data (tab_auger(i,684),i=1,16) /26,11,3,1048.7,440.9,22,25,9975,0,0,0,0,0,0,0,0/
data (tab_auger(i,685),i=1,16) /26,11,4,1036,425.7,21,48,9952,0,0,0,0,0,0,0,0/
data (tab_auger(i,686),i=1,16) /26,12,1,7636,4094.7,4,236,3345,1299,5120,0,0,0,0,0,0/
data (tab_auger(i,687),i=1,16) /26,12,2,1199,458.1,28,798,9202,0,0,0,0,0,0,0,0/
data (tab_auger(i,688),i=1,16) /26,12,3,1085.7,406.8,29,33,9967,0,0,0,0,0,0,0,0/
data (tab_auger(i,689),i=1,16) /26,12,4,1073,387.2,25,112,9888,0,0,0,0,0,0,0,0/
data (tab_auger(i,690),i=1,16) /26,13,1,7686,3969.6,1,297,3336,1792,4575,0,0,0,0,0,0/
data (tab_auger(i,691),i=1,16) /26,13,2,1234,425,12,714,9286,0,0,0,0,0,0,0,0/
data (tab_auger(i,692),i=1,16) /26,13,3,1122.6,365.3,13,46,9954,0,0,0,0,0,0,0,0/
data (tab_auger(i,693),i=1,16) /26,13,4,1110,311.8,13,662,9338,0,0,0,0,0,0,0,0/
data (tab_auger(i,694),i=1,16) /26,14,1,7737,3746.4,1,313,3390,6297,0,0,0,0,0,0,0/
data (tab_auger(i,695),i=1,16) /26,14,2,1270,399.8,27,548,9452,0,0,0,0,0,0,0,0/
data (tab_auger(i,696),i=1,16) /26,14,3,1159.6,319.7,18,157,9843,0,0,0,0,0,0,0,0/
data (tab_auger(i,697),i=1,16) /26,14,4,1147,266.8,35,954,9046,0,0,0,0,0,0,0,0/
data (tab_auger(i,698),i=1,16) /26,15,1,7788,3633.6,-2,453,3985,5562,0,0,0,0,0,0,0/
data (tab_auger(i,699),i=1,16) /26,15,2,1307,373.7,12,19,9981,0,0,0,0,0,0,0,0/
data (tab_auger(i,700),i=1,16) /26,15,3,1197.6,232.9,12,1263,8737,0,0,0,0,0,0,0,0/
data (tab_auger(i,701),i=1,16) /26,15,4,1185,221.9,12,1263,8737,0,0,0,0,0,0,0,0/
data (tab_auger(i,702),i=1,16) /26,16,1,7838,3349.1,-1,3638,6362,0,0,0,0,0,0,0,0/
data (tab_auger(i,703),i=1,16) /26,17,1,7891,3288.9,-1,3692,6308,0,0,0,0,0,0,0,0/
data (tab_auger(i,704),i=1,16) /26,18,1,7989,3097.2,3,4015,5985,0,0,0,0,0,0,0,0/
data (tab_auger(i,705),i=1,16) /26,19,1,8088,2840.5,3,4404,5596,0,0,0,0,0,0,0,0/
data (tab_auger(i,706),i=1,16) /26,20,1,8187,2535.3,7,4887,5113,0,0,0,0,0,0,0,0/
data (tab_auger(i,707),i=1,16) /26,21,1,8286,2153.6,4,5521,4479,0,0,0,0,0,0,0,0/
data (tab_auger(i,708),i=1,16) /26,22,1,8384,2380.1,8,4903,5097,0,0,0,0,0,0,0,0/
data (tab_auger(i,709),i=1,16) /26,23,1,8482,4525,4,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,710),i=1,16) /27,1,1,7715,4903.6,10,73,1249,1082,1493,1083,1462,1513,1379,633,37/
data (tab_auger(i,711),i=1,16) /27,1,2,931,832.3,130,0,224,1242,2248,3175,2005,1021,85,0,0/
data (tab_auger(i,712),i=1,16) /27,1,3,800.1,752.7,53,197,2916,3075,3060,717,35,0,0,0,0/
data (tab_auger(i,713),i=1,16) /27,1,4,785,733.3,58,240,2404,3115,3401,801,39,0,0,0,0/
data (tab_auger(i,714),i=1,16) /27,1,5,107,85.4,321,0,739,9261,0,0,0,0,0,0,0/
data (tab_auger(i,715),i=1,16) /27,1,6,68.3,49.5,90,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,716),i=1,16) /27,1,7,66,47.1,93,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,717),i=1,16) /27,2,1,7763,4905.3,9,77,1418,1082,1595,3298,1196,1210,123,1,0/
data (tab_auger(i,718),i=1,16) /27,2,2,963,786.9,91,0,263,1703,2706,4987,332,9,0,0,0/
data (tab_auger(i,719),i=1,16) /27,2,3,833.1,743.6,50,178,3293,3047,3034,448,0,0,0,0,0/
data (tab_auger(i,720),i=1,16) /27,2,4,818,720.1,52,263,2988,3070,3173,506,0,0,0,0,0/
data (tab_auger(i,721),i=1,16) /27,2,5,130,67.2,287,0,366,9634,0,0,0,0,0,0,0/
data (tab_auger(i,722),i=1,16) /27,2,6,91.3,56.1,168,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,723),i=1,16) /27,2,7,89,53.7,172,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,724),i=1,16) /27,3,1,7811,4861.5,9,74,1254,1104,3435,1893,1981,257,2,0,0/
data (tab_auger(i,725),i=1,16) /27,3,2,996,744.1,77,0,267,1691,7095,931,16,0,0,0,0/
data (tab_auger(i,726),i=1,16) /27,3,3,866.1,715.8,50,197,2927,3131,3745,0,0,0,0,0,0/
data (tab_auger(i,727),i=1,16) /27,3,4,851,691.8,54,241,2412,3184,4163,0,0,0,0,0,0/
data (tab_auger(i,728),i=1,16) /27,3,5,153,36.1,294,0,304,9696,0,0,0,0,0,0,0/
data (tab_auger(i,729),i=1,16) /27,3,6,114.3,45.6,141,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,730),i=1,16) /27,3,7,112,43.3,144,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,731),i=1,16) /27,4,1,7859,4804.8,7,69,1070,2601,3023,2444,737,56,0,0,0/
data (tab_auger(i,732),i=1,16) /27,4,2,1028,717.3,62,0,322,8114,1525,39,0,0,0,0,0/
data (tab_auger(i,733),i=1,16) /27,4,3,899.1,675.7,45,220,2496,6418,866,0,0,0,0,0,0/
data (tab_auger(i,734),i=1,16) /27,4,4,884,649.6,49,206,1777,7029,988,0,0,0,0,0,0/
data (tab_auger(i,735),i=1,16) /27,4,5,176,69.8,140,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,736),i=1,16) /27,4,6,137.3,31.1,179,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,737),i=1,16) /27,4,7,135,28.8,182,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,738),i=1,16) /27,5,1,7906,4746.6,5,60,3260,1355,3109,1996,206,14,0,0,0/
data (tab_auger(i,739),i=1,16) /27,5,2,1060,666.4,66,0,357,8612,1031,0,0,0,0,0,0/
data (tab_auger(i,740),i=1,16) /27,5,3,933.1,645.8,30,249,8644,1100,7,0,0,0,0,0,0/
data (tab_auger(i,741),i=1,16) /27,5,4,918,620.6,36,152,8471,1368,9,0,0,0,0,0,0/
data (tab_auger(i,742),i=1,16) /27,5,5,199,43.9,133,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,743),i=1,16) /27,5,6,161.3,6.3,164,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,744),i=1,16) /27,5,7,159,4,167,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,745),i=1,16) /27,6,1,7954,4693.4,0,457,3315,1160,3997,1071,0,0,0,0,0/
data (tab_auger(i,746),i=1,16) /27,6,2,1092,606.8,30,76,1095,8829,0,0,0,0,0,0,0/
data (tab_auger(i,747),i=1,16) /27,6,3,967,612.2,12,285,9715,0,0,0,0,0,0,0,0/
data (tab_auger(i,748),i=1,16) /27,6,4,952,590.2,11,73,9927,0,0,0,0,0,0,0,0/
data (tab_auger(i,749),i=1,16) /27,6,5,222,11.8,93,158,9842,0,0,0,0,0,0,0,0/
data (tab_auger(i,750),i=1,16) /27,7,1,8002,4656,2,456,3342,1485,4717,0,0,0,0,0,0/
data (tab_auger(i,751),i=1,16) /27,7,2,1124,564.4,26,135,1639,8226,0,0,0,0,0,0,0/
data (tab_auger(i,752),i=1,16) /27,7,3,1001,586.4,9,264,9736,0,0,0,0,0,0,0,0/
data (tab_auger(i,753),i=1,16) /27,7,4,986,569.1,7,76,9924,0,0,0,0,0,0,0,0/
data (tab_auger(i,754),i=1,16) /27,8,1,8050,4606.6,4,447,3351,1479,4723,0,0,0,0,0,0/
data (tab_auger(i,755),i=1,16) /27,8,2,1156,598.5,3,848,9152,0,0,0,0,0,0,0,0/
data (tab_auger(i,756),i=1,16) /27,8,3,1034,559.9,6,212,9788,0,0,0,0,0,0,0,0/
data (tab_auger(i,757),i=1,16) /27,8,4,1019,546.8,5,61,9939,0,0,0,0,0,0,0,0/
data (tab_auger(i,758),i=1,16) /27,9,1,8098,4560.7,5,435,3364,1462,4739,0,0,0,0,0,0/
data (tab_auger(i,759),i=1,16) /27,9,2,1188,570,20,925,9075,0,0,0,0,0,0,0,0/
data (tab_auger(i,760),i=1,16) /27,9,3,1067,534.7,24,133,9867,0,0,0,0,0,0,0,0/
data (tab_auger(i,761),i=1,16) /27,9,4,1052,523.6,26,44,9956,0,0,0,0,0,0,0,0/
data (tab_auger(i,762),i=1,16) /27,10,1,8146,4509.6,0,419,3380,1432,4769,0,0,0,0,0,0/
data (tab_auger(i,763),i=1,16) /27,10,2,1221,542.3,5,1014,8986,0,0,0,0,0,0,0,0/
data (tab_auger(i,764),i=1,16) /27,10,3,1101,515.6,5,22,9978,0,0,0,0,0,0,0,0/
data (tab_auger(i,765),i=1,16) /27,10,4,1086,502.7,7,22,9978,0,0,0,0,0,0,0,0/
data (tab_auger(i,766),i=1,16) /27,11,1,8198,4438.8,2,362,3441,1326,4871,0,0,0,0,0,0/
data (tab_auger(i,767),i=1,16) /27,11,2,1258,528.2,8,993,9007,0,0,0,0,0,0,0,0/
data (tab_auger(i,768),i=1,16) /27,11,3,1138,499.2,12,26,9974,0,0,0,0,0,0,0,0/
data (tab_auger(i,769),i=1,16) /27,11,4,1123,484.4,13,34,9966,0,0,0,0,0,0,0,0/
data (tab_auger(i,770),i=1,16) /27,12,1,8250,4350.6,3,306,3506,1248,4940,0,0,0,0,0,0/
data (tab_auger(i,771),i=1,16) /27,12,2,1295,508.8,11,960,9040,0,0,0,0,0,0,0,0/
data (tab_auger(i,772),i=1,16) /27,12,3,1176,474.6,20,32,9968,0,0,0,0,0,0,0,0/
data (tab_auger(i,773),i=1,16) /27,12,4,1161,456.7,20,61,9939,0,0,0,0,0,0,0,0/
data (tab_auger(i,774),i=1,16) /27,13,1,8303,4251.4,4,260,3571,1320,4849,0,0,0,0,0,0/
data (tab_auger(i,775),i=1,16) /27,13,2,1332,483.7,26,907,9093,0,0,0,0,0,0,0,0/
data (tab_auger(i,776),i=1,16) /27,13,3,1215,438.6,27,41,9959,0,0,0,0,0,0,0,0/
data (tab_auger(i,777),i=1,16) /27,13,4,1200,415.4,24,142,9858,0,0,0,0,0,0,0,0/
data (tab_auger(i,778),i=1,16) /27,14,1,8355,4113.6,1,363,3534,1909,4194,0,0,0,0,0,0/
data (tab_auger(i,779),i=1,16) /27,14,2,1369,450.3,11,814,9186,0,0,0,0,0,0,0,0/
data (tab_auger(i,780),i=1,16) /27,14,3,1253.9,394.4,12,59,9941,0,0,0,0,0,0,0,0/
data (tab_auger(i,781),i=1,16) /27,14,4,1239,329.9,12,840,9160,0,0,0,0,0,0,0,0/
data (tab_auger(i,782),i=1,16) /27,15,1,8408,3885.4,1,400,3588,6012,0,0,0,0,0,0,0/
data (tab_auger(i,783),i=1,16) /27,15,2,1407,424.2,25,630,9370,0,0,0,0,0,0,0,0/
data (tab_auger(i,784),i=1,16) /27,15,3,1292.9,344.6,17,200,9800,0,0,0,0,0,0,0,0/
data (tab_auger(i,785),i=1,16) /27,15,4,1278,279.9,32,1196,8804,0,0,0,0,0,0,0,0/
data (tab_auger(i,786),i=1,16) /27,16,1,8461,3762.6,-2,598,4184,5218,0,0,0,0,0,0,0/
data (tab_auger(i,787),i=1,16) /27,16,2,1446,400.1,11,27,9973,0,0,0,0,0,0,0,0/
data (tab_auger(i,788),i=1,16) /27,16,3,1332.9,244.7,10,1560,8440,0,0,0,0,0,0,0,0/
data (tab_auger(i,789),i=1,16) /27,16,4,1318,232.1,10,1560,8440,0,0,0,0,0,0,0,0/
data (tab_auger(i,790),i=1,16) /27,17,1,8513,3462,-1,3887,6113,0,0,0,0,0,0,0,0/
data (tab_auger(i,791),i=1,16) /27,18,1,8569,3398.5,-1,3942,6058,0,0,0,0,0,0,0,0/
data (tab_auger(i,792),i=1,16) /27,19,1,8671,3187.7,3,4276,5724,0,0,0,0,0,0,0,0/
data (tab_auger(i,793),i=1,16) /27,20,1,8774,2909.2,3,4678,5322,0,0,0,0,0,0,0,0/
data (tab_auger(i,794),i=1,16) /27,21,1,8877,2576.9,6,5177,4823,0,0,0,0,0,0,0,0/
data (tab_auger(i,795),i=1,16) /27,22,1,8980,2157,3,5839,4161,0,0,0,0,0,0,0,0/
data (tab_auger(i,796),i=1,16) /27,23,1,9082,2384.3,7,5267,4733,0,0,0,0,0,0,0,0/
data (tab_auger(i,797),i=1,16) /27,24,1,9184,4887,4,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,798),i=1,16) /28,1,1,8338,4999.7,11,114,1571,1234,1600,1166,1411,1366,1084,412,44/
data (tab_auger(i,799),i=1,16) /28,1,2,1015,909.9,113,0,292,1558,2420,3049,1778,843,60,0,0/
data (tab_auger(i,800),i=1,16) /28,1,3,877.4,823.1,47,240,3137,3068,2870,655,30,0,0,0,0/
data (tab_auger(i,801),i=1,16) /28,1,4,860,797.1,48,344,2929,3081,2939,676,31,0,0,0,0/
data (tab_auger(i,802),i=1,16) /28,1,5,117,93,292,0,795,9205,0,0,0,0,0,0,0/
data (tab_auger(i,803),i=1,16) /28,1,6,75.8,55.5,73,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,804),i=1,16) /28,1,7,73,52.7,76,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,805),i=1,16) /28,2,1,8387,4968.9,7,116,1737,1227,1711,3057,1820,287,43,2,0/
data (tab_auger(i,806),i=1,16) /28,2,2,1048,882.3,77,0,376,2267,3113,4065,179,0,0,0,0/
data (tab_auger(i,807),i=1,16) /28,2,3,912.4,816.1,44,220,3452,3049,2889,390,0,0,0,0,0/
data (tab_auger(i,808),i=1,16) /28,2,4,895,787.7,44,363,3441,3010,2792,394,0,0,0,0,0/
data (tab_auger(i,809),i=1,16) /28,2,5,141,74.9,283,0,407,9593,0,0,0,0,0,0,0/
data (tab_auger(i,810),i=1,16) /28,2,6,99.8,61.4,149,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,811),i=1,16) /28,2,7,97,58.5,153,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,812),i=1,16) /28,3,1,8436,4920.5,6,114,1577,1262,3403,2771,726,141,6,0,0/
data (tab_auger(i,813),i=1,16) /28,3,2,1082,826.6,66,0,362,2134,6688,808,8,0,0,0,0/
data (tab_auger(i,814),i=1,16) /28,3,3,947.4,791,45,241,3147,3136,3476,0,0,0,0,0,0/
data (tab_auger(i,815),i=1,16) /28,3,4,930,763.8,46,345,2937,3153,3565,0,0,0,0,0,0/
data (tab_auger(i,816),i=1,16) /28,3,5,166,43.1,241,0,359,9641,0,0,0,0,0,0,0/
data (tab_auger(i,817),i=1,16) /28,3,6,124.8,52.8,145,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,818),i=1,16) /28,3,7,122,50,148,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,819),i=1,16) /28,4,1,8486,4863.1,5,110,1398,2738,3725,1386,543,100,0,0,0/
data (tab_auger(i,820),i=1,16) /28,4,2,1116,792.2,45,0,431,8200,1343,26,0,0,0,0,0/
data (tab_auger(i,821),i=1,16) /28,4,3,982.4,753.7,38,265,2798,6155,782,0,0,0,0,0,0/
data (tab_auger(i,822),i=1,16) /28,4,4,965,726.7,40,316,2381,6465,838,0,0,0,0,0,0/
data (tab_auger(i,823),i=1,16) /28,4,5,190,78,97,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,824),i=1,16) /28,4,6,148.8,36.8,124,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,825),i=1,16) /28,4,7,146,34,126,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,826),i=1,16) /28,5,1,8535,4805.5,5,103,3712,1387,2626,1668,473,31,0,0,0/
data (tab_auger(i,827),i=1,16) /28,5,2,1150,740.8,57,0,468,8819,713,0,0,0,0,0,0/
data (tab_auger(i,828),i=1,16) /28,5,3,1017.4,717.6,18,295,8696,1006,3,0,0,0,0,0,0/
data (tab_auger(i,829),i=1,16) /28,5,4,1000,691.1,17,270,8560,1166,4,0,0,0,0,0,0/
data (tab_auger(i,830),i=1,16) /28,5,5,215,59.5,131,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,831),i=1,16) /28,5,6,173.8,18.4,162,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,832),i=1,16) /28,5,7,171,15.6,164,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,833),i=1,16) /28,6,1,8585,4749.4,4,555,3705,1132,3529,991,88,0,0,0,0/
data (tab_auger(i,834),i=1,16) /28,6,2,1184,684.3,61,72,1102,8826,0,0,0,0,0,0,0/
data (tab_auger(i,835),i=1,16) /28,6,3,1053.4,686.3,27,332,9668,0,0,0,0,0,0,0,0/
data (tab_auger(i,836),i=1,16) /28,6,4,1036,661.7,32,200,9800,0,0,0,0,0,0,0,0/
data (tab_auger(i,837),i=1,16) /28,6,5,239,27.9,125,26,9974,0,0,0,0,0,0,0,0/
data (tab_auger(i,838),i=1,16) /28,7,1,8635,4703.5,-1,542,3739,1433,4286,0,0,0,0,0,0/
data (tab_auger(i,839),i=1,16) /28,7,2,1218,625.7,27,122,1415,8463,0,0,0,0,0,0,0/
data (tab_auger(i,840),i=1,16) /28,7,3,1089.4,650.3,10,392,9608,0,0,0,0,0,0,0,0/
data (tab_auger(i,841),i=1,16) /28,7,4,1072,631.4,9,113,9887,0,0,0,0,0,0,0,0/
data (tab_auger(i,842),i=1,16) /28,8,1,8685,4651.6,1,533,3750,1427,4290,0,0,0,0,0,0/
data (tab_auger(i,843),i=1,16) /28,8,2,1252,656.9,4,906,9094,0,0,0,0,0,0,0,0/
data (tab_auger(i,844),i=1,16) /28,8,3,1125.4,623.8,7,347,9653,0,0,0,0,0,0,0,0/
data (tab_auger(i,845),i=1,16) /28,8,4,1108,608.9,6,97,9903,0,0,0,0,0,0,0,0/
data (tab_auger(i,846),i=1,16) /28,9,1,8735,4601.7,4,519,3764,1413,4304,0,0,0,0,0,0/
data (tab_auger(i,847),i=1,16) /28,9,2,1286,627.9,2,978,9022,0,0,0,0,0,0,0,0/
data (tab_auger(i,848),i=1,16) /28,9,3,1160.4,596.4,4,276,9724,0,0,0,0,0,0,0,0/
data (tab_auger(i,849),i=1,16) /28,9,4,1143,584.8,3,78,9922,0,0,0,0,0,0,0,0/
data (tab_auger(i,850),i=1,16) /28,10,1,8785,4555,5,501,3782,1389,4328,0,0,0,0,0,0/
data (tab_auger(i,851),i=1,16) /28,10,2,1320,597.9,19,1063,8937,0,0,0,0,0,0,0,0/
data (tab_auger(i,852),i=1,16) /28,10,3,1195.4,571.7,22,172,9828,0,0,0,0,0,0,0,0/
data (tab_auger(i,853),i=1,16) /28,10,4,1178,560.6,24,54,9946,0,0,0,0,0,0,0,0/
data (tab_auger(i,854),i=1,16) /28,11,1,8835,4504.3,0,476,3804,1348,4372,0,0,0,0,0,0/
data (tab_auger(i,855),i=1,16) /28,11,2,1355,568.5,5,1161,8839,0,0,0,0,0,0,0,0/
data (tab_auger(i,856),i=1,16) /28,11,3,1231.5,554,5,25,9975,0,0,0,0,0,0,0,0/
data (tab_auger(i,857),i=1,16) /28,11,4,1214,539,7,25,9975,0,0,0,0,0,0,0,0/
data (tab_auger(i,858),i=1,16) /28,12,1,8889,4431.3,2,413,3875,1250,4462,0,0,0,0,0,0/
data (tab_auger(i,859),i=1,16) /28,12,2,1394,553,7,1137,8863,0,0,0,0,0,0,0,0/
data (tab_auger(i,860),i=1,16) /28,12,3,1270.5,535.6,10,30,9970,0,0,0,0,0,0,0,0/
data (tab_auger(i,861),i=1,16) /28,12,4,1253,518.5,12,40,9960,0,0,0,0,0,0,0,0/
data (tab_auger(i,862),i=1,16) /28,13,1,8943,4340.3,3,351,3949,1185,4515,0,0,0,0,0,0/
data (tab_auger(i,863),i=1,16) /28,13,2,1433,534.1,9,1100,8900,0,0,0,0,0,0,0,0/
data (tab_auger(i,864),i=1,16) /28,13,3,1310.5,510.9,19,37,9963,0,0,0,0,0,0,0,0/
data (tab_auger(i,865),i=1,16) /28,13,4,1293,490,19,71,9929,0,0,0,0,0,0,0,0/
data (tab_auger(i,866),i=1,16) /28,14,1,8998,4235.5,4,303,4019,1282,4396,0,0,0,0,0,0/
data (tab_auger(i,867),i=1,16) /28,14,2,1472,509.1,26,1039,8961,0,0,0,0,0,0,0,0/
data (tab_auger(i,868),i=1,16) /28,14,3,1351.5,472.3,27,49,9951,0,0,0,0,0,0,0,0/
data (tab_auger(i,869),i=1,16) /28,14,4,1334,445.3,23,167,9833,0,0,0,0,0,0,0,0/
data (tab_auger(i,870),i=1,16) /28,15,1,9053,4087.6,1,454,3945,1902,3699,0,0,0,0,0,0/
data (tab_auger(i,871),i=1,16) /28,15,2,1511,473.6,11,934,9066,0,0,0,0,0,0,0,0/
data (tab_auger(i,872),i=1,16) /28,15,3,1392.5,424.3,12,69,9931,0,0,0,0,0,0,0,0/
data (tab_auger(i,873),i=1,16) /28,15,4,1375,347.9,12,994,9006,0,0,0,0,0,0,0,0/
data (tab_auger(i,874),i=1,16) /28,16,1,9108,3861.1,1,515,3988,5497,0,0,0,0,0,0,0/
data (tab_auger(i,875),i=1,16) /28,16,2,1551,447.5,24,725,9275,0,0,0,0,0,0,0,0/
data (tab_auger(i,876),i=1,16) /28,16,3,1433.5,370,16,236,9764,0,0,0,0,0,0,0,0/
data (tab_auger(i,877),i=1,16) /28,16,4,1416,293.5,30,1400,8600,0,0,0,0,0,0,0,0/
data (tab_auger(i,878),i=1,16) /28,17,1,9163,3730.6,-2,779,4517,4704,0,0,0,0,0,0,0/
data (tab_auger(i,879),i=1,16) /28,17,2,1592,426.3,11,37,9963,0,0,0,0,0,0,0,0/
data (tab_auger(i,880),i=1,16) /28,17,3,1475.5,257,9,1802,8198,0,0,0,0,0,0,0,0/
data (tab_auger(i,881),i=1,16) /28,17,4,1458,242.7,10,1802,8198,0,0,0,0,0,0,0,0/
data (tab_auger(i,882),i=1,16) /28,18,1,9217,3414.3,-1,4382,5618,0,0,0,0,0,0,0,0/
data (tab_auger(i,883),i=1,16) /28,19,1,9275,3347.4,-1,4438,5562,0,0,0,0,0,0,0,0/
data (tab_auger(i,884),i=1,16) /28,20,1,9381,3117.9,2,4783,5217,0,0,0,0,0,0,0,0/
data (tab_auger(i,885),i=1,16) /28,21,1,9488,2821,3,5193,4807,0,0,0,0,0,0,0,0/
data (tab_auger(i,886),i=1,16) /28,22,1,9595,2470.9,5,5695,4305,0,0,0,0,0,0,0,0/
data (tab_auger(i,887),i=1,16) /28,23,1,9702,2032,3,6355,3645,0,0,0,0,0,0,0,0/
data (tab_auger(i,888),i=1,16) /28,24,1,9808,2257.6,6,5836,4164,0,0,0,0,0,0,0,0/
data (tab_auger(i,889),i=1,16) /28,25,1,9914,5262,4,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,890),i=1,16) /29,1,1,8986,5177.4,8,170,2012,1352,1651,1344,1495,1357,585,35,0/
data (tab_auger(i,891),i=1,16) /29,1,2,1103,985.2,88,0,423,2202,2695,2761,1424,494,1,0,0/
data (tab_auger(i,892),i=1,16) /29,1,3,958,896,36,274,3805,3091,2296,514,20,0,0,0,0/
data (tab_auger(i,893),i=1,16) /29,1,4,938,855.4,35,514,3839,2986,2155,487,19,0,0,0,0/
data (tab_auger(i,894),i=1,16) /29,1,5,127,98.5,260,0,905,9095,0,0,0,0,0,0,0/
data (tab_auger(i,895),i=1,16) /29,1,6,83.3,61.7,78,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,896),i=1,16) /29,1,7,80,58.4,81,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,897),i=1,16) /29,2,1,9035,5135.8,6,171,2015,1345,1771,2220,2059,418,2,0,0/
data (tab_auger(i,898),i=1,16) /29,2,2,1136,957.9,70,0,500,2621,3238,3291,350,0,0,0,0/
data (tab_auger(i,899),i=1,16) /29,2,3,995,889.8,38,274,3810,3079,2549,288,0,0,0,0,0/
data (tab_auger(i,900),i=1,16) /29,2,4,975,849.8,37,515,3843,2974,2395,273,0,0,0,0,0/
data (tab_auger(i,901),i=1,16) /29,2,5,153,80.4,248,0,644,9356,0,0,0,0,0,0,0/
data (tab_auger(i,902),i=1,16) /29,2,6,109.3,65,117,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,903),i=1,16) /29,2,7,106,61.7,121,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,904),i=1,16) /29,3,1,9086,5085.6,6,171,1867,1393,2680,2720,1155,13,1,0,0/
data (tab_auger(i,905),i=1,16) /29,3,2,1171,907.2,63,0,520,2673,5225,1569,13,0,0,0,0/
data (tab_auger(i,906),i=1,16) /29,3,3,1032.1,869.1,38,298,3545,3180,2977,0,0,0,0,0,0/
data (tab_auger(i,907),i=1,16) /29,3,4,1012,831.1,39,502,3406,3139,2953,0,0,0,0,0,0/
data (tab_auger(i,908),i=1,16) /29,3,5,179,51.9,235,0,577,9423,0,0,0,0,0,0,0/
data (tab_auger(i,909),i=1,16) /29,3,6,135.3,58.8,128,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,910),i=1,16) /29,3,7,132,55.5,132,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,911),i=1,16) /29,4,1,9137,5024.9,4,169,1702,2305,3290,2310,193,30,1,0,0/
data (tab_auger(i,912),i=1,16) /29,4,2,1207,868,52,0,595,5800,3356,249,0,0,0,0,0/
data (tab_auger(i,913),i=1,16) /29,4,3,1069.1,833.2,38,326,3247,4314,2113,0,0,0,0,0,0/
data (tab_auger(i,914),i=1,16) /29,4,4,1049,796.6,36,479,2923,5428,1170,0,0,0,0,0,0/
data (tab_auger(i,915),i=1,16) /29,4,5,205,88.5,100,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,916),i=1,16) /29,4,6,161.4,44.8,127,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,917),i=1,16) /29,4,7,158,41.5,130,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,918),i=1,16) /29,5,1,9188,4968.6,4,163,2914,2081,3637,930,208,67,0,0,0/
data (tab_auger(i,919),i=1,16) /29,5,2,1243,815.4,42,0,612,8267,1121,0,0,0,0,0,0/
data (tab_auger(i,920),i=1,16) /29,5,3,1106.1,793.7,32,360,4199,4767,674,0,0,0,0,0,0/
data (tab_auger(i,921),i=1,16) /29,5,4,1086,760.7,24,439,6915,2568,78,0,0,0,0,0,0/
data (tab_auger(i,922),i=1,16) /29,5,5,231,68.5,88,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,923),i=1,16) /29,5,6,187.4,24.9,109,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,924),i=1,16) /29,5,7,184,21.5,111,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,925),i=1,16) /29,6,1,9240,4917.2,4,486,3622,1599,3093,911,289,0,0,0,0/
data (tab_auger(i,926),i=1,16) /29,6,2,1279,760.1,53,48,532,9420,0,0,0,0,0,0,0/
data (tab_auger(i,927),i=1,16) /29,6,3,1143.1,764.2,15,401,8645,954,0,0,0,0,0,0,0/
data (tab_auger(i,928),i=1,16) /29,6,4,1123,731.8,16,377,8491,1132,0,0,0,0,0,0,0/
data (tab_auger(i,929),i=1,16) /29,6,5,257,45.9,121,10,9990,0,0,0,0,0,0,0,0/
data (tab_auger(i,930),i=1,16) /29,6,6,213.4,2.4,145,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,931),i=1,16) /29,7,1,9294,4871.6,4,634,3951,1456,3959,0,0,0,0,0,0/
data (tab_auger(i,932),i=1,16) /29,7,2,1315,702.8,58,119,1485,8396,0,0,0,0,0,0,0/
data (tab_auger(i,933),i=1,16) /29,7,3,1181.2,730.8,25,453,9547,0,0,0,0,0,0,0,0/
data (tab_auger(i,934),i=1,16) /29,7,4,1161,703.8,30,281,9719,0,0,0,0,0,0,0,0/
data (tab_auger(i,935),i=1,16) /29,7,5,283,10,113,29,9971,0,0,0,0,0,0,0,0/
data (tab_auger(i,936),i=1,16) /29,8,1,9344,4806.8,-2,613,3969,1430,3988,0,0,0,0,0,0/
data (tab_auger(i,937),i=1,16) /29,8,2,1351,716.9,4,975,9025,0,0,0,0,0,0,0,0/
data (tab_auger(i,938),i=1,16) /29,8,3,1219.2,689.2,8,534,9466,0,0,0,0,0,0,0,0/
data (tab_auger(i,939),i=1,16) /29,8,4,1199,672.5,6,157,9843,0,0,0,0,0,0,0,0/
data (tab_auger(i,940),i=1,16) /29,9,1,9396,4753.2,1,600,3984,1419,3997,0,0,0,0,0,0/
data (tab_auger(i,941),i=1,16) /29,9,2,1387,687.6,4,1050,8950,0,0,0,0,0,0,0,0/
data (tab_auger(i,942),i=1,16) /29,9,3,1257.2,661.7,7,477,9523,0,0,0,0,0,0,0,0/
data (tab_auger(i,943),i=1,16) /29,9,4,1237,649.3,6,135,9865,0,0,0,0,0,0,0,0/
data (tab_auger(i,944),i=1,16) /29,10,1,9448,4701.3,3,581,4003,1399,4017,0,0,0,0,0,0/
data (tab_auger(i,945),i=1,16) /29,10,2,1423,655.8,1,1137,8863,0,0,0,0,0,0,0,0/
data (tab_auger(i,946),i=1,16) /29,10,3,1294.2,632.8,3,383,9617,0,0,0,0,0,0,0,0/
data (tab_auger(i,947),i=1,16) /29,10,4,1274,623.1,1,107,9893,0,0,0,0,0,0,0,0/
data (tab_auger(i,948),i=1,16) /29,11,1,9500,4653.6,4,553,4028,1365,4054,0,0,0,0,0,0/
data (tab_auger(i,949),i=1,16) /29,11,2,1459,623.4,18,1238,8762,0,0,0,0,0,0,0,0/
data (tab_auger(i,950),i=1,16) /29,11,3,1331.3,609.1,21,240,9760,0,0,0,0,0,0,0,0/
data (tab_auger(i,951),i=1,16) /29,11,4,1311,598.4,23,74,9926,0,0,0,0,0,0,0,0/
data (tab_auger(i,952),i=1,16) /29,12,1,9552,4602.5,0,515,4061,1304,4120,0,0,0,0,0,0/
data (tab_auger(i,953),i=1,16) /29,12,2,1496,590.8,4,1357,8643,0,0,0,0,0,0,0,0/
data (tab_auger(i,954),i=1,16) /29,12,3,1369.3,593.1,4,32,9968,0,0,0,0,0,0,0,0/
data (tab_auger(i,955),i=1,16) /29,12,4,1349,575.8,6,32,9968,0,0,0,0,0,0,0,0/
data (tab_auger(i,956),i=1,16) /29,13,1,9609,4526.6,2,448,4137,1214,4201,0,0,0,0,0,0/
data (tab_auger(i,957),i=1,16) /29,13,2,1537,574.2,6,1329,8671,0,0,0,0,0,0,0,0/
data (tab_auger(i,958),i=1,16) /29,13,3,1410.3,573.7,9,38,9962,0,0,0,0,0,0,0,0/
data (tab_auger(i,959),i=1,16) /29,13,4,1390,553.8,10,49,9951,0,0,0,0,0,0,0,0/
data (tab_auger(i,960),i=1,16) /29,14,1,9665,4431.1,3,383,4214,1164,4239,0,0,0,0,0,0/
data (tab_auger(i,961),i=1,16) /29,14,2,1578,555.4,8,1287,8713,0,0,0,0,0,0,0,0/
data (tab_auger(i,962),i=1,16) /29,14,3,1452.3,547.2,18,47,9953,0,0,0,0,0,0,0,0/
data (tab_auger(i,963),i=1,16) /29,14,4,1432,522.8,18,88,9912,0,0,0,0,0,0,0,0/
data (tab_auger(i,964),i=1,16) /29,15,1,9722,4319.9,3,338,4284,1292,4086,0,0,0,0,0,0/
data (tab_auger(i,965),i=1,16) /29,15,2,1619,529.9,24,1217,8783,0,0,0,0,0,0,0,0/
data (tab_auger(i,966),i=1,16) /29,15,3,1495.4,505.6,25,60,9940,0,0,0,0,0,0,0,0/
data (tab_auger(i,967),i=1,16) /29,15,4,1475,473.7,22,206,9794,0,0,0,0,0,0,0,0/
data (tab_auger(i,968),i=1,16) /29,16,1,9779,4160.3,1,545,4168,1953,3334,0,0,0,0,0,0/
data (tab_auger(i,969),i=1,16) /29,16,2,1660,493.8,9,1097,8903,0,0,0,0,0,0,0,0/
data (tab_auger(i,970),i=1,16) /29,16,3,1538.4,454.1,10,86,9914,0,0,0,0,0,0,0,0/
data (tab_auger(i,971),i=1,16) /29,16,4,1518,363.6,10,1188,8812,0,0,0,0,0,0,0,0/
data (tab_auger(i,972),i=1,16) /29,17,1,9836,3932.3,1,637,4198,5165,0,0,0,0,0,0,0/
data (tab_auger(i,973),i=1,16) /29,17,2,1702,469.3,22,857,9143,0,0,0,0,0,0,0,0/
data (tab_auger(i,974),i=1,16) /29,17,3,1581.4,395.7,15,293,9707,0,0,0,0,0,0,0,0/
data (tab_auger(i,975),i=1,16) /29,17,4,1561,305.2,27,1659,8341,0,0,0,0,0,0,0,0/
data (tab_auger(i,976),i=1,16) /29,18,1,9893,3792.3,-2,979,4685,4336,0,0,0,0,0,0,0/
data (tab_auger(i,977),i=1,16) /29,18,2,1745,453.4,10,49,9951,0,0,0,0,0,0,0,0/
data (tab_auger(i,978),i=1,16) /29,18,3,1625.4,266.7,9,2120,7880,0,0,0,0,0,0,0,0/
data (tab_auger(i,979),i=1,16) /29,18,4,1605,250.6,9,2120,7880,0,0,0,0,0,0,0,0/
data (tab_auger(i,980),i=1,16) /29,19,1,9949,3462.6,-1,4677,5323,0,0,0,0,0,0,0,0/
data (tab_auger(i,981),i=1,16) /29,20,1,10009,3392.4,-1,4734,5266,0,0,0,0,0,0,0,0/
data (tab_auger(i,982),i=1,16) /29,21,1,10119,3145.4,2,5085,4915,0,0,0,0,0,0,0,0/
data (tab_auger(i,983),i=1,16) /29,22,1,10230,2829.4,3,5499,4501,0,0,0,0,0,0,0,0/
data (tab_auger(i,984),i=1,16) /29,23,1,10341,2456.1,5,6007,3993,0,0,0,0,0,0,0,0/
data (tab_auger(i,985),i=1,16) /29,24,1,10452,1986.2,2,6678,3322,0,0,0,0,0,0,0,0/
data (tab_auger(i,986),i=1,16) /29,25,1,10562,2201.7,5,6215,3785,0,0,0,0,0,0,0,0/
data (tab_auger(i,987),i=1,16) /29,26,1,10672,5651,4,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,988),i=1,16) /30,1,1,9663,5582.7,12,232,2013,1417,1659,1338,1341,994,684,258,65/
data (tab_auger(i,989),i=1,16) /30,1,2,1198,1064.8,87,0,526,2185,2752,2683,1303,505,46,0,0/
data (tab_auger(i,990),i=1,16) /30,1,3,1047.2,973.4,32,362,3884,3118,2174,444,18,0,0,0,0/
data (tab_auger(i,991),i=1,16) /30,1,4,1024,917.7,32,699,3739,3012,2100,432,18,0,0,0,0/
data (tab_auger(i,992),i=1,16) /30,1,5,141,108.4,202,0,1524,8476,0,0,0,0,0,0,0/
data (tab_auger(i,993),i=1,16) /30,1,6,94.9,72.6,53,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,994),i=1,16) /30,1,7,91,68.7,55,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,995),i=1,16) /30,2,1,9708,5556.6,14,232,2016,1411,1719,1333,1604,1416,257,13,0/
data (tab_auger(i,996),i=1,16) /30,2,2,1228,1027.5,90,0,556,2315,2955,2700,1212,262,0,0,0/
data (tab_auger(i,997),i=1,16) /30,2,3,1081.2,966.8,38,362,3888,3106,2282,351,11,0,0,0,0/
data (tab_auger(i,998),i=1,16) /30,2,4,1058,911.2,37,701,3742,3001,2204,342,10,0,0,0,0/
data (tab_auger(i,999),i=1,16) /30,2,5,165,93.1,232,0,1184,8816,0,0,0,0,0,0,0/
data (tab_auger(i,1000),i=1,16) /30,2,6,118.9,71.1,83,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,1001),i=1,16) /30,2,7,115,67.1,85,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,1002),i=1,16) /30,3,1,9760,5504,11,233,2018,1453,1883,2805,1583,22,3,0,0/
data (tab_auger(i,1003),i=1,16) /30,3,2,1264,992.4,74,0,650,2891,3168,3168,123,0,0,0,0/
data (tab_auger(i,1004),i=1,16) /30,3,3,1120.3,954.8,42,363,3892,3194,2551,0,0,0,0,0,0/
data (tab_auger(i,1005),i=1,16) /30,3,4,1097,899.5,41,701,3746,3086,2467,0,0,0,0,0,0/
data (tab_auger(i,1006),i=1,16) /30,3,5,192,65.4,254,0,767,9233,0,0,0,0,0,0,0/
data (tab_auger(i,1007),i=1,16) /30,3,6,145.9,66.5,135,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,1008),i=1,16) /30,3,7,142,62.6,139,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,1009),i=1,16) /30,4,1,9813,5444.1,10,233,1877,1449,3325,2876,203,36,1,0,0/
data (tab_auger(i,1010),i=1,16) /30,4,2,1301,954.5,63,0,791,3086,5542,581,0,0,0,0,0/
data (tab_auger(i,1011),i=1,16) /30,4,3,1159.3,923.1,42,395,3636,3190,2779,0,0,0,0,0,0/
data (tab_auger(i,1012),i=1,16) /30,4,4,1136,870.6,42,685,3333,3148,2834,0,0,0,0,0,0/
data (tab_auger(i,1013),i=1,16) /30,4,5,220,101.2,105,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,1014),i=1,16) /30,4,6,174,55.2,133,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,1015),i=1,16) /30,4,7,170,51.2,136,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,1016),i=1,16) /30,5,1,9866,5386.5,9,230,1712,2856,4030,871,275,24,2,0,0/
data (tab_auger(i,1017),i=1,16) /30,5,2,1339,900.2,50,0,757,8118,1125,0,0,0,0,0,0/
data (tab_auger(i,1018),i=1,16) /30,5,3,1198.3,885.1,38,433,3326,5630,611,0,0,0,0,0,0/
data (tab_auger(i,1019),i=1,16) /30,5,4,1175,835.4,39,654,2859,5833,654,0,0,0,0,0,0/
data (tab_auger(i,1020),i=1,16) /30,5,5,247,81.8,103,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,1021),i=1,16) /30,5,6,201,35.8,126,1,9999,0,0,0,0,0,0,0,0/
data (tab_auger(i,1022),i=1,16) /30,5,7,197,31.8,129,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,1023),i=1,16) /30,6,1,9919,5334.7,9,222,3988,1621,2717,1231,221,0,0,0,0/
data (tab_auger(i,1024),i=1,16) /30,6,2,1377,846.9,45,0,789,9211,0,0,0,0,0,0,0/
data (tab_auger(i,1025),i=1,16) /30,6,3,1237.3,851.3,23,479,8702,819,0,0,0,0,0,0,0/
data (tab_auger(i,1026),i=1,16) /30,6,4,1214,806.4,24,602,8457,941,0,0,0,0,0,0,0/
data (tab_auger(i,1027),i=1,16) /30,6,5,275,58.9,95,0,10000,0,0,0,0,0,0,0,0/
data (tab_auger(i,1028),i=1,16) /30,6,6,229,13,114,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,1029),i=1,16) /30,6,7,225,9,116,2,9998,0,0,0,0,0,0,0,0/
data (tab_auger(i,1030),i=1,16) /30,7,1,9973,5295.6,9,713,3919,1312,3195,861,0,0,0,0,0/
data (tab_auger(i,1031),i=1,16) /30,7,2,1415,786.6,57,85,1345,8140,430,0,0,0,0,0,0/
data (tab_auger(i,1032),i=1,16) /30,7,3,1276.4,819.5,20,535,9465,0,0,0,0,0,0,0,0/
data (tab_auger(i,1033),i=1,16) /30,7,4,1253,780.3,18,518,9482,0,0,0,0,0,0,0,0/
data (tab_auger(i,1034),i=1,16) /30,7,5,303,34.9,132,17,9983,0,0,0,0,0,0,0,0/
data (tab_auger(i,1035),i=1,16) /30,8,1,10028,5235.5,10,696,3934,1544,3826,0,0,0,0,0,0/
data (tab_auger(i,1036),i=1,16) /30,8,2,1453,807.2,26,993,9007,0,0,0,0,0,0,0,0/
data (tab_auger(i,1037),i=1,16) /30,8,3,1316.4,776.8,27,621,9379,0,0,0,0,0,0,0,0/
data (tab_auger(i,1038),i=1,16) /30,8,4,1293,749.9,33,410,9590,0,0,0,0,0,0,0,0/
data (tab_auger(i,1039),i=1,16) /30,9,1,10082,5170.9,5,658,3966,1505,3871,0,0,0,0,0,0/
data (tab_auger(i,1040),i=1,16) /30,9,2,1491,763.5,14,1122,8878,0,0,0,0,0,0,0,0/
data (tab_auger(i,1041),i=1,16) /30,9,3,1356.4,739,17,712,9288,0,0,0,0,0,0,0,0/
data (tab_auger(i,1042),i=1,16) /30,9,4,1333,728.6,17,215,9785,0,0,0,0,0,0,0,0/
data (tab_auger(i,1043),i=1,16) /30,10,1,10136,5118.3,8,640,3984,1491,3885,0,0,0,0,0,0/
data (tab_auger(i,1044),i=1,16) /30,10,2,1529,731.3,12,1211,8789,0,0,0,0,0,0,0,0/
data (tab_auger(i,1045),i=1,16) /30,10,3,1396.5,710.5,16,641,9359,0,0,0,0,0,0,0,0/
data (tab_auger(i,1046),i=1,16) /30,10,4,1373,704.6,16,184,9816,0,0,0,0,0,0,0,0/
data (tab_auger(i,1047),i=1,16) /30,11,1,10190,5072.7,11,615,4007,1464,3914,0,0,0,0,0,0/
data (tab_auger(i,1048),i=1,16) /30,11,2,1567,699,11,1314,8686,0,0,0,0,0,0,0,0/
data (tab_auger(i,1049),i=1,16) /30,11,3,1435.5,684.9,13,520,9480,0,0,0,0,0,0,0,0/
data (tab_auger(i,1050),i=1,16) /30,11,4,1412,680.3,13,146,9854,0,0,0,0,0,0,0,0/
data (tab_auger(i,1051),i=1,16) /30,12,1,10244,5027,12,577,4040,1415,3968,0,0,0,0,0,0/
data (tab_auger(i,1052),i=1,16) /30,12,2,1605,666.6,29,1434,8566,0,0,0,0,0,0,0,0/
data (tab_auger(i,1053),i=1,16) /30,12,3,1474.5,667.2,34,329,9671,0,0,0,0,0,0,0,0/
data (tab_auger(i,1054),i=1,16) /30,12,4,1451,658.9,37,99,9901,0,0,0,0,0,0,0,0/
data (tab_auger(i,1055),i=1,16) /30,13,1,10298,4971.3,8,523,4085,1327,4065,0,0,0,0,0,0/
data (tab_auger(i,1056),i=1,16) /30,13,2,1644,630.6,16,1575,8425,0,0,0,0,0,0,0,0/
data (tab_auger(i,1057),i=1,16) /30,13,3,1514.6,656.9,19,39,9961,0,0,0,0,0,0,0,0/
data (tab_auger(i,1058),i=1,16) /30,13,4,1491,636.6,22,39,9961,0,0,0,0,0,0,0,0/
data (tab_auger(i,1059),i=1,16) /30,14,1,10357,4883,8,457,4160,1242,4141,0,0,0,0,0,0/
data (tab_auger(i,1060),i=1,16) /30,14,2,1687,611.9,16,1545,8455,0,0,0,0,0,0,0,0/
data (tab_auger(i,1061),i=1,16) /30,14,3,1557.6,633.9,22,47,9953,0,0,0,0,0,0,0,0/
data (tab_auger(i,1062),i=1,16) /30,14,4,1534,608.8,23,61,9939,0,0,0,0,0,0,0,0/
data (tab_auger(i,1063),i=1,16) /30,15,1,10415,4792.8,11,394,4235,1204,4167,0,0,0,0,0,0/
data (tab_auger(i,1064),i=1,16) /30,15,2,1730,592,18,1498,8502,0,0,0,0,0,0,0,0/
data (tab_auger(i,1065),i=1,16) /30,15,3,1601.6,603.7,30,58,9942,0,0,0,0,0,0,0,0/
data (tab_auger(i,1066),i=1,16) /30,15,4,1578,575.2,30,109,9891,0,0,0,0,0,0,0,0/
data (tab_auger(i,1067),i=1,16) /30,16,1,10474,4681.9,12,357,4301,1372,3970,0,0,0,0,0,0/
data (tab_auger(i,1068),i=1,16) /30,16,2,1773,556.5,26,1421,8579,0,0,0,0,0,0,0,0/
data (tab_auger(i,1069),i=1,16) /30,16,3,1646.7,545.4,27,75,9925,0,0,0,0,0,0,0,0/
data (tab_auger(i,1070),i=1,16) /30,16,4,1623,508.5,25,254,9746,0,0,0,0,0,0,0,0/
data (tab_auger(i,1071),i=1,16) /30,17,1,10533,4505.8,10,617,4152,2096,3135,0,0,0,0,0,0/
data (tab_auger(i,1072),i=1,16) /30,17,2,1816,533.4,20,1287,8713,0,0,0,0,0,0,0,0/
data (tab_auger(i,1073),i=1,16) /30,17,3,1691.7,508.8,24,106,9894,0,0,0,0,0,0,0,0/
data (tab_auger(i,1074),i=1,16) /30,17,4,1668,398.5,21,1409,8591,0,0,0,0,0,0,0,0/
data (tab_auger(i,1075),i=1,16) /30,18,1,10592,4266.9,10,738,4178,5084,0,0,0,0,0,0,0/
data (tab_auger(i,1076),i=1,16) /30,18,2,1860,510.3,32,1013,8987,0,0,0,0,0,0,0,0/
data (tab_auger(i,1077),i=1,16) /30,18,3,1736.7,443.1,27,362,9638,0,0,0,0,0,0,0,0/
data (tab_auger(i,1078),i=1,16) /30,18,4,1713,333.5,36,1948,8052,0,0,0,0,0,0,0,0/
data (tab_auger(i,1079),i=1,16) /30,19,1,10651,4101.8,6,1147,4673,4180,0,0,0,0,0,0,0/
data (tab_auger(i,1080),i=1,16) /30,19,2,1905,503.1,21,64,9936,0,0,0,0,0,0,0,0/
data (tab_auger(i,1081),i=1,16) /30,19,3,1782.7,291.2,17,2470,7530,0,0,0,0,0,0,0,0/
data (tab_auger(i,1082),i=1,16) /30,19,4,1759,273.4,18,2470,7530,0,0,0,0,0,0,0,0/
data (tab_auger(i,1083),i=1,16) /30,20,1,10709,3741.1,5,4702,5298,0,0,0,0,0,0,0,0/
data (tab_auger(i,1084),i=1,16) /30,21,1,10772,3666.2,4,4758,5242,0,0,0,0,0,0,0,0/
data (tab_auger(i,1085),i=1,16) /30,22,1,10886,3370.4,5,5112,4888,0,0,0,0,0,0,0,0/
data (tab_auger(i,1086),i=1,16) /30,23,1,11001,3026.7,5,5532,4468,0,0,0,0,0,0,0,0/
data (tab_auger(i,1087),i=1,16) /30,24,1,11116,2612.9,6,6049,3951,0,0,0,0,0,0,0,0/
data (tab_auger(i,1088),i=1,16) /30,25,1,11231,2101.2,4,6742,3258,0,0,0,0,0,0,0,0/
data (tab_auger(i,1089),i=1,16) /30,26,1,11345,2310.7,6,6322,3678,0,0,0,0,0,0,0,0/
data (tab_auger(i,1090),i=1,16) /30,27,1,11459,6123,10,0,10000,0,0,0,0,0,0,0,0/


end


!     funcs-auger.f

!     compute quantities related to Auger ionization and store

!     routines:
!     SUBROUTINE getQs(k,j) - obtain Auger ionization yields
!     FUNCTION AUG_Rate(k,j) - obtain the Auger ionization rates

!
!.........................................................................
!

      SUBROUTINE getQs(k,j,auger_rate2)


!     compute the Auger Q coefficients

!     see Eq 10 from Verner + Yakovlev (1990,A&SS,165,27)

!     k   = the species number
!     j   = the ionization stage of the pre ionized ion
!     m   = the ionization stage of the post ionized ion
!     s   = the shell number
!     ne  = the number of electrons ejected going from j to m

!     m is always greater than j+2 because j+1 yields ne=1, and ne=1 is
!     the non-Auger single electron photoionization, which is stored
!     globally as R_ph and was computed in function PH_Rate

!     the partial photoionization rate for shell "s", R_phs, was
!     also computed in function PH_Rate

!     the quantity W_auger is the probability that ne electrons are
!     ejected from shell "s" from initial ion k,j (which subsequently
!     becomes ion k,m)
       
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

      implicit none
      include           'rates.h'
      integer           k,j,m,s,ne,nek,nout
      double precision  asum,auger_rate2
      include           'rates.com'
      include           'com-auger.com'

!     we need the total number of populated shells, so include the ntot
!     common block

      integer           ntot
      COMMON/ntot/      ntot(30)
      auger_rate2=auger_rate2

!     for the initial ion k,j, determine the number of bound electrons,
!     nek, and the index of the most outer shell populated by one or
!     more electron(s) for this number of electrons, nout; nout is
!     determined from Verner's ntot array; it is assumed the ion is in
!     the ground state

      nek  = k - j + 1
      nout = ntot(nek)

!     loop over the final ionization stages (m>=j+2) following the Auger
!     process from initial ionization stage j; compute the number of
!     ejected electrons; null the partial rate running sum (asum), loop
!     over the populated shells of the initial ion (j) and sum the
!     partial rates to obtain the total rate coefficient, Q

      DO 11 m=j+2,k+1
       ne    = m - j
       asum  = 0.0d0
       IF (ne.le.10) then  ! tables stop at 10 ejected electrons
        DO 15 s=1,nout
!         asum = asum + W_auger(k,j,s,ne)*R_phs(s) !ORIGINAL
 15     CONTINUE
       END IF
!       Q(k,j,m) = asum !ORIGINAL
 11   CONTINUE


      END


!
!.........................................................................
!

      DOUBLE PRECISION FUNCTION AUG_Rate(k,j)

!     compute the Auger rate for destrcution of k,j

!     k   = the species number
!     j   = the ionization stage of the pre ionized ion
!     m   = the ionization stage of the post ionized ion

!     m is always greater than j+2 because j+1 yields ne=1, and ne=1 is
!     the non-Auger single electron photoionization, which is stored
!     globally as R_ph and was computed in function PH_Rate

!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

      implicit none
      include           'rates.h'
      integer           k,j,m
      double precision  asum
      include           'rates.com'


      asum =  0.0d0

!     sum over the final ionization stages m out of which initial
!     ionization stage j can go to following Auger ionization of 2 or
!     more electrons (where the photoionized electron is the first
!     electron)

      DO 11 m=j+2,k+1
       asum = asum + Q(k,j,m)
 11   CONTINUE

      AUG_Rate = asum

      RETURN
      END




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine optical_depth_convolved_ioneqther( nemax, bxs2 ,energy,vturb,cross_section_convolved)

implicit none  
integer,parameter :: nion=168, out_unit=20
integer :: i,nemax,aa
double precision ::   bxs2(0:nemax), vturb 

!Constants for fftw3
integer, parameter :: fftw_r2hc=0, fftw_hc2r=1, fftw_dht=2, fftw_redft00=3, fftw_redft01=4
integer, parameter :: fftw_redft10=5, fftw_redft11=6, fftw_rodft00=7, fftw_rodft01=8, fftw_rodft10=9
integer, parameter :: fftw_rodft11=10, fftw_forward=-1, fftw_backward=+1, fftw_measure=0, fftw_destroy_input=1
integer, parameter :: fftw_unaligned=2, fftw_conserve_memory=4, fftw_exhaustive=8, fftw_preserve_input=16
integer, parameter :: fftw_patient=32, fftw_estimate=64, fftw_wisdom_only=2097152, fftw_estimate_patient=128
integer, parameter :: fftw_believe_pcost=256, fftw_no_dft_r2hc=512, fftw_no_nonthreaded=1024, fftw_no_buffering=2048
integer, parameter :: fftw_no_indirect_op=4096, fftw_allow_large_generic=8192, fftw_no_rank_splits=16384
integer, parameter :: fftw_no_vrecurse=65536, fftw_no_simd=131072, fftw_no_slow=262144  
integer, parameter :: fftw_no_fixed_radix_large_n=524288, fftw_allow_pruning=1048576, fftw_no_vrank_splits=32768

!For fourier transform and convolution f*g
integer :: plan
double complex,allocatable :: f_o(:),g_o(:),f_1(:),g_1(:),h_o(:),h_1(:)
double precision :: gauss_norm,delta_v,A,cc,x,B, tmp3,y_o 
double precision,parameter ::   PI=3.141592654, c=299792.458 
double precision :: energy(0:nemax) , gauss(-nemax:nemax),cross_section(nemax),cross_section_convolved(0:nemax)

!Preparing arrays for fourier transform only once (I will use complex numbers with z=a+0*i) 
logical :: startup=.true.

if(startup)then
 allocate(f_o(nemax))
 allocate(g_o(nemax))
 allocate(f_1(nemax))
 allocate(g_1(nemax))
 allocate(h_o(nemax))
 allocate(h_1(nemax))
 startup=.false.  
endif 

!Gaussian to do the convolution
A=energy(1)
y_o=energy(nemax/2)
B=10**(((log10(energy(nemax)))-(log10(energy(1))))/nemax) 
delta_v=(vturb/c)*y_o 
gauss_norm=1/(delta_v*sqrt(2*PI)) 
!Gaussian a*exp(-((x-b)**2)/(2*cc*cc))
do i=1,nemax
 cc=1/(2*delta_v*delta_v) 
 if(energy(i).lt.energy(nemax/2))then
  x=((energy(nemax/2)-energy(i))) 
 else if (energy(i).ge.energy(nemax/2))then
  x=((energy(i)-energy(nemax/2)))
 endif
  gauss(i)=gauss_norm*exp(-((x)**2)*cc)
enddo 
!!!!!!!!!FOURIER TRANSFORM!!!!!!!!!!!!!!!!!!! 
!Fourier transform of cross_section 
do i=1,nemax
 cross_section(i) =bxs2(i)
enddo 
do i=1,nemax
 f_o(i) = dcmplx(cross_section(i),0.0)
enddo 
call dfftw_plan_dft_1d(plan,nemax,f_o,f_1,fftw_forward,fftw_estimate)
call dfftw_execute_dft(plan,f_o,f_1)
call dfftw_destroy_plan(plan) 
!Fourier transform of gaussian
tmp3=0.0 
do i=1,nemax
 tmp3 =tmp3+ gauss(i)
enddo 
do i=1,nemax
 g_o(i) = dcmplx(gauss(i)/tmp3,0.0)
enddo 
call dfftw_plan_dft_1d(plan,nemax,g_o,g_1,fftw_forward,fftw_estimate)
call dfftw_execute_dft(plan,g_o,g_1)
call dfftw_destroy_plan(plan)   
!Multiplication
do i=1,nemax
 h_o(i)=f_1(i)*(conjg(g_1(i)))
enddo   
! Inverse Fourier transform
call dfftw_plan_dft_1d(plan,nemax,h_o,h_1,fftw_backward,fftw_estimate)
call dfftw_execute_dft(plan,h_o,h_1)
call dfftw_destroy_plan(plan) 
aa=nemax/2   
do i=1,nemax/2
 cross_section_convolved(i) = real(h_1(aa)/((real(nemax))))
 aa=aa+1
enddo 
aa=1
do i=nemax/2,nemax
 cross_section_convolved(i) = real(h_1(aa)/((real(nemax))))
 aa=aa+1
enddo 
call dfftw_cleanup()

end subroutine optical_depth_convolved_ioneqther 


! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine interpolate_cross_section_ioneqther(xs,j,ener,maxi,bener2,bxs_restored,zn,ii)
! xs -> old cross section
! j --> ion
! ener --> old energy
! maxi --> size of old energy
! bener2 --> new energy
! bxs_restored -->new cross-section
! zn --> atomic number
! ii --> atomic state (1=neutral, 2=singly ionized, 3=double ionized, etc.)

implicit none

integer,parameter :: nene=650000, out_unit=20,nion=168
integer :: maxi,   i, j, k, r, zn, ii
double precision :: stemp, etemp,s ,etemp2 
double precision :: ener(0:nion,maxi), xs(0:nion,maxi),bxs_restored(0:nion,nene) 
double precision :: bener2(nene) 
 
real SS,thresh
thresh=0.0

k=2  
do i=1,nene-2 
 r=k
	!! Linear interpolation
 etemp2=(bener2(i)+bener2(i-1))/2
 s=xs(j,k-1)+(xs(j,k)-xs(j,k-1))*(etemp2-ener(j,k-1))/(ener(j,k)-ener(j,k-1))
 stemp=0.d0
 etemp=0.d0  
 do while (bener2(i).gt.ener(j,k).and.ener(j,k).lt.bener2(i+1).and.bener2(i).lt.ener(j,maxi))
           ! Integral under area!!!    
  stemp=stemp+(xs(j,k))*(ener(j,k)-ener(j,k-1))
  etemp=etemp+(ener(j,k)-ener(j,k-1))
  k=k+1
 enddo   
 if(r.lt.k)then
  s=real(stemp/etemp)
 endif 
 if(bener2(i).eq.ener(j,k))then
  s=xs(j,k)
 endif
 if(bener2(i).lt.ener(j,1))then
  call phfit2_b(zn,zn-ii+1,1,REAL(bener2(i)),SS,thresh) 
  s=DBLE(SS*1.d-18)
 endif
 if(bener2(i).gt.ener(j,maxi))then
  call phfit2_b(zn,zn-ii+1,1,REAL(bener2(i)),SS,thresh) 
  s=DBLE(SS*1.d-18)
 endif
 bxs_restored(j,i)=s 
enddo 

end subroutine interpolate_cross_section_ioneqther

