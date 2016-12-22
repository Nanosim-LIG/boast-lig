program grenoble

  use def_parameters
  use mod_nsi_assembly, only : nsi_assembly
  implicit none  
  integer(ip) :: pnode,pgaus,itype,mnode,ntens,ncomp_nsi
  integer(ip) :: ielem
  real(rp)    :: amatr(1000),rhsid(1000)
  !
  ! Choose element type
  !
  itype     = HEX08

  mnode     = 8
  ntens     = 3*ndime-3
  ncomp_nsi = 3

  if(      itype == HEX08 ) then
     ! HEX08
     pnode = 8
     pgaus = 8
  else if( itype == TET04 ) then
     ! TET04
     pnode = 4
     pgaus = 1
  else if( itype == PYR05 ) then
     ! PYR05
     pnode = 5
     pgaus = 5
  else if( itype == PEN06 ) then
     ! PEN06
     pnode = 6
     pgaus = 6
  end if
  !
  ! Assembly
  !
  do ielem = 1,20000
     call nsi_assembly(mnode,pnode,pgaus,ntens,ncomp_nsi,amatr,rhsid)
  end do

end program grenoble

