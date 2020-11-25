!     com-file.com

!     global variables (strings) for important path locations
!     (these variables are defined in the 'rates.inp' file)

!     homedir    = directory where tab*.dat files live (fitting coefficients)
!     seduvbdir  = directory where UVB SED library lives
!     sedsb99dir = directory where Starburt99 SED library lives 

      character(LEN=255) :: homedir,seduvbdir,sedsb99dir
      COMMON/tabdir/  homedir,seduvbdir,sedsb99dir

! eof
