module def_parameters
  !
  ! Parameters that do not change
  !
  integer,     parameter  :: ip    = 4               ! 4-byte integer
  integer,     parameter  :: rp    = 8               ! Double precision 
  integer(ip), parameter  :: ndime = 3               ! Spatial dimension  
  real(rp),    parameter  :: zeror = epsilon(1.0_rp) ! Almost zero
  integer(ip), parameter  :: TET04 = 30              ! 3D 
  integer(ip), parameter  :: TET10 = 31              ! 3D 
  integer(ip), parameter  :: PYR05 = 32              ! 3D 
  integer(ip), parameter  :: PYR14 = 33              ! 3D 
  integer(ip), parameter  :: PEN06 = 34              ! 3D  
  integer(ip), parameter  :: PEN15 = 35              ! 3D 
  integer(ip), parameter  :: PEN18 = 36              ! 3D 
  integer(ip), parameter  :: HEX08 = 37              ! 3D 
  integer(ip), parameter  :: HEX20 = 38              ! 3D 
  integer(ip), parameter  :: HEX27 = 39              ! 3D 
  integer(ip), parameter  :: HEX64 = 40              ! 3D 
  integer(ip), parameter  :: SHELL = 51              ! 3D shell element
  integer(ip), parameter  :: BAR3D = 52              ! 3D bar element
  !
  ! Size of vector: "VECTOR_SIZE" elements are assembled at the same time
  !
  integer(ip), parameter  :: VECTOR_SIZE = 4         ! Size for vectorization
  
end module def_parameters
