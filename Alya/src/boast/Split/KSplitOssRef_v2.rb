
require_relative './KSplitOss.rb'
 
def generate_mocks
  mocks = <<EOF
  function pabdf_nsi(x) result(y)
    integer,intent(in) :: x 
    real :: y
    y = 1.0
  end function pabdf_nsi
EOF
return mocks  
end
def generate_ref_declaration
     decl_ref =<<EOF
     subroutine nsi_element_assembly_split_oss_ref(ndime,&
     kfl_lumped,mnode,ntens,kfl_duatss,fact_duatss,&
     kfl_stabi_nsi,fvins_nsi,fcons_nsi,bemol_nsi,&
     kfl_regim_nsi,fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi,&
     kfl_p1ve2_nsi,kfl_linea_nsi,kfl_confi_nsi,nbdfp_nsi,&
     kfl_sgsti_nsi,kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi,&
     penal_nsi,kfl_bubbl_nsi,pnode,pgaus,gpden,gpvis,gppor,&
     gpsp1,gpsp2,gpvol,gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,&
     gprhc,gpvel,gpsgs,elvel,elpre,elbub,wgrgr,agrau,elauu,&
     elaup,elapp,elapu,elrbu,elrbp,dtinv_loc,dtsgs,pbubl,&
     gpsha_bub,gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq,elrbq)

     #{gen_def_parameters}

     integer(ip), intent(in)    :: ndime 
     integer(ip), intent(in)    :: kfl_lumped
     integer(ip), intent(in)    :: mnode
     integer(ip), intent(in)    :: ntens
     integer(ip), intent(in)    :: kfl_duatss
     integer(ip), intent(in)    :: fact_duatss
     integer(ip), intent(in)    :: kfl_stabi_nsi
     real(rp),    intent(in)    :: fvins_nsi
     real(rp),    intent(in)    :: fcons_nsi
     real(rp),    intent(in)    :: bemol_nsi
     integer(ip), intent(in)    :: kfl_regim_nsi
     real(rp),    intent(in)    :: fvela_nsi(3)
     integer(ip), intent(in)    :: kfl_rmom2_nsi
     integer(ip), intent(in)    :: kfl_press_nsi
     integer(ip), intent(in)    :: kfl_p1ve2_nsi
     integer(ip), intent(in)    :: kfl_linea_nsi
     integer(ip), intent(in)    :: kfl_confi_nsi
     integer(ip), intent(in)    :: nbdfp_nsi
     integer(ip), intent(in)    :: kfl_sgsti_nsi
     integer(ip), intent(in)    :: kfl_nota1_nsi
     integer(ip), intent(in)    :: kfl_limit_nsi
     integer(ip), intent(in)    :: kfl_penal_nsi
     real(rp),    intent(in)    :: penal_nsi
     integer(ip), intent(in)    :: kfl_bubbl_nsi

     integer(ip), intent(in)    :: pnode,pgaus
     real(rp),    intent(in)    :: gpden(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpvis(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gppor(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpsp1(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpsp2(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpvol(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpsha(#{$p_vector_length},pnode,pgaus)
     real(rp),    intent(in)    :: gpcar(#{$p_vector_length},ndime,mnode,pgaus)
     real(rp),    intent(in)    :: gpadv(#{$p_vector_length},ndime,pgaus)
     real(rp),    intent(inout) :: gpvep(#{$p_vector_length},ndime,pgaus)
     real(rp),    intent(inout) :: gpgrp(#{$p_vector_length},ndime,pgaus)
     real(rp),    intent(inout) :: gprhs(#{$p_vector_length},ndime,pgaus)
     real(rp),    intent(inout) :: gprhc(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpvel(#{$p_vector_length},ndime,pgaus,*)
     real(rp),    intent(in)    :: gpsgs(#{$p_vector_length},ndime,pgaus,*)
     real(rp),    intent(in)    :: elvel(#{$p_vector_length},ndime,pnode,*)
     real(rp),    intent(in)    :: elpre(#{$p_vector_length},pnode,*)
     real(rp),    intent(in)    :: elbub(#{$p_vector_length})

     real(rp),    intent(out)   :: wgrgr(#{$p_vector_length},pnode,pnode,pgaus)
     real(rp),    intent(out)   :: agrau(#{$p_vector_length},pnode,pgaus)

     ! Matrices
     real(rp),    intent(out)   :: elauu(#{$p_vector_length},pnode*ndime,pnode*ndime)
     real(rp),    intent(out)   :: elaup(#{$p_vector_length},pnode*ndime,pnode)
     real(rp),    intent(out)   :: elapp(#{$p_vector_length},pnode,pnode)
     real(rp),    intent(out)   :: elapu(#{$p_vector_length},pnode,pnode*ndime)
     real(rp),    intent(out)   :: elrbu(#{$p_vector_length},ndime,pnode)
     real(rp),    intent(out)   :: elrbp(#{$p_vector_length},pnode)
     ! Others
     real(rp),    intent(in)    :: dtinv_loc(#{$p_vector_length})
     real(rp),    intent(in)    :: dtsgs(#{$p_vector_length})
     integer(ip), intent(in)    :: pbubl(#{$p_vector_length})
     real(rp),    intent(in)    :: gpsha_bub(#{$p_vector_length},pgaus)
     real(rp),    intent(in)    :: gpcar_bub(#{$p_vector_length},ndime,pgaus)
     ! Enrichement Element matrices
     real(rp),    intent(out)   :: elauq(#{$p_vector_length},pnode*ndime,1)
     real(rp),    intent(out)   :: elapq(#{$p_vector_length},pnode,1)
     real(rp),    intent(out)   :: elaqu(#{$p_vector_length},1,pnode*ndime)
     real(rp),    intent(out)   :: elaqp(#{$p_vector_length},1,pnode)
     real(rp),    intent(out)   :: elaqq(#{$p_vector_length},1,1)
     real(rp),    intent(out)   :: elrbq(#{$p_vector_length},1)
     ! Local arrays
     real(rp)                   :: gpsp1_p(#{$p_vector_length},pgaus)
     real(rp)                   :: gpsp1_v(#{$p_vector_length},pgaus)
     real(rp)                   :: gpsp2_v(#{$p_vector_length},pgaus)
     real(rp)                   :: c1(#{$p_vector_length})
     real(rp)                   :: c2(#{$p_vector_length})
     real(rp)                   :: c3(#{$p_vector_length})
     real(rp)                   :: c4(#{$p_vector_length})
     real(rp)                   :: alpha(#{$p_vector_length})
     real(rp)                   :: beta(#{$p_vector_length})
     real(rp)                   :: fact0(#{$p_vector_length})
     real(rp)                   :: fact1(#{$p_vector_length})
     real(rp)                   :: fact2(#{$p_vector_length})
     real(rp)                   :: fact3(#{$p_vector_length})
     real(rp)                   :: fact4(#{$p_vector_length})
     real(rp)                   :: fact5(#{$p_vector_length})
     real(rp)                   :: fact6(#{$p_vector_length})
     real(rp)                   :: fact7(#{$p_vector_length})
     real(rp)                   :: fact8(#{$p_vector_length})
     real(rp)                   :: gpveo(#{$p_vector_length},3)
     real(rp)                   :: fact1_p(#{$p_vector_length})
     real(rp)                   :: dtinv_mod(#{$p_vector_length})
     integer(ip)                :: inode,jnode,jdime
     integer(ip)                :: idofv,jdof2,jdof3,ivect
     integer(ip)                :: idof1,idof3,idof2,igaus
     integer(ip)                :: idime,jdof1,jdofv,itime
EOF
  if @opts[:preprocessor] then
    decl_ref = decl_ref + <<EOF
    #ifdef OPENACC
    #define DEF_VECT ivect
    #else
    #define DEF_VECT 1:#{$p_vector_length}
    #endif
EOF
  end
  return decl_ref 
end
def gen_def_parameters
  parameters =<<EOF
  integer,     parameter  :: ip    = 4               ! 4-byte integer
  integer,     parameter  :: rp    = 8               ! Double precision 

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
EOF
  return parameters
end
def generate_ref_initialization
init = <<EOF
   dtinv_mod = dtinv_loc

   gpsp1_p = gpsp1
   gpsp1_v = gpsp1
   gpsp2_v = gpsp2

   if( kfl_nota1_nsi == 1 ) gpsp1_v = 0.0_rp 

   elrbp = 0.0_rp
   elrbu = 0.0_rp
   elapp = 0.0_rp
   elauu = 0.0_rp
   elaup = 0.0_rp
   elapu = 0.0_rp
EOF
  return init
end 

class KSplitOssRef_v2 < KSplitOss
 def generate
   # opts = {:vector_length => 1, :preprocessor => false, :nests => (1..10).to_a}
   # opts.update(options)
   macros = ""
   
   if @opts[:preprocessor] then
     $p_vector_length = "VECTOR_SIZE"
     $p_def_vect = "DEF_VECT"
     macros = <<EOF
     #define VECTOR_SIZE #{@opts[:vector_length]}
EOF
   else
     $p_vector_length = "#{@opts[:vector_length]}"
     $p_def_vect = "1:#{@opts[:vector_length]}"
   end
 
   nests = []
   nests.push <<EOF
   agrau(#{$p_def_vect},:,:)   = 0.0_rp
       wgrgr(#{$p_def_vect},:,:,:) = 0.0_rp 
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                agrau(#{$p_def_vect},inode,igaus) =  agrau(#{$p_def_vect},inode,igaus) + &
                                              gpadv(#{$p_def_vect},idime,igaus) * gpcar(#{$p_def_vect},idime,inode,igaus)
             end do
             agrau(#{$p_def_vect},inode,igaus) =  gpden(#{$p_def_vect},igaus) * agrau(#{$p_def_vect},inode,igaus) 
             do jnode = 1,pnode
                do idime = 1,ndime
                   wgrgr(#{$p_def_vect},inode,jnode,igaus) = wgrgr(#{$p_def_vect},inode,jnode,igaus) + &
                                                      gpcar(#{$p_def_vect},idime,inode,igaus)*gpcar(#{$p_def_vect},idime,jnode,igaus)
                end do
             end do
          end do
       end do
EOF
   nests.push <<EOF
   do igaus = 1,pgaus
   
      fact0(#{$p_def_vect}) = gpsp2_v(#{$p_def_vect},igaus) * gpvol(#{$p_def_vect},igaus)
      fact6(#{$p_def_vect}) = gpvis(#{$p_def_vect},igaus)   * gpvol(#{$p_def_vect},igaus)
      fact7(#{$p_def_vect}) = gpsp1_v(#{$p_def_vect},igaus) * gpvol(#{$p_def_vect},igaus)
      fact8(#{$p_def_vect}) = pabdf_nsi(1) * gpden(#{$p_def_vect},igaus) * dtinv_mod(#{$p_def_vect}) + gppor(#{$p_def_vect},igaus)
   
      do inode = 1,pnode
         do idime = 1,ndime
   
            idofv           = (inode-1)*ndime+idime
            fact1(#{$p_def_vect}) = fact0(#{$p_def_vect}) * gpcar(#{$p_def_vect},idime,inode,igaus)      
   
            do jnode = 1,pnode    
               do jdime = 1,ndime                   
                  jdofv                       = (jnode-1)*ndime+jdime
                  elauu(#{$p_def_vect},idofv,jdofv) = elauu(#{$p_def_vect},idofv,jdofv) + fact1(#{$p_def_vect}) * gpcar(#{$p_def_vect},jdime,jnode,igaus)                   
               end do
   
               jdofv           = (jnode-1)*ndime+idime
               fact4(#{$p_def_vect}) = gpsha(#{$p_def_vect},inode,igaus) * gpvol(#{$p_def_vect},igaus)
               fact5(#{$p_def_vect}) = fact4(#{$p_def_vect}) * ( agrau(#{$p_def_vect},jnode,igaus) + fact8(#{$p_def_vect}) * gpsha(#{$p_def_vect},jnode,igaus) ) + fact6(#{$p_def_vect}) *   wgrgr(#{$p_def_vect},inode,jnode,igaus) + fact7(#{$p_def_vect}) *   agrau(#{$p_def_vect},jnode,igaus) * agrau(#{$p_def_vect},inode,igaus)   
               elauu(#{$p_def_vect},idofv,jdofv) = elauu(#{$p_def_vect},idofv,jdofv) + fact5(#{$p_def_vect})
   
            end do
         end do
      end do
   end do
EOF
   nests.push <<EOF
   if( fvins_nsi > 0.9_rp ) then
      do igaus = 1,pgaus
         do inode = 1,pnode
            do idime = 1,ndime
               idofv = (inode-1)*ndime + idime
               do jnode = 1,pnode
                  fact1(#{$p_def_vect}) = gpvis(#{$p_def_vect},igaus) * gpvol(#{$p_def_vect},igaus) * gpcar(#{$p_def_vect},idime,jnode,igaus)     
                  do jdime = 1,ndime
                     jdofv                       = (jnode-1)*ndime + jdime
                     elauu(#{$p_def_vect},idofv,jdofv) = elauu(#{$p_def_vect},idofv,jdofv) + fact1(#{$p_def_vect}) * gpcar(#{$p_def_vect},jdime,inode,igaus)
                  end do
               end do
               if( fvins_nsi == 2.0_rp ) then
                  fact1(#{$p_def_vect}) = -2.0_rp / 3.0_rp * gpvis(#{$p_def_vect},igaus) * gpvol(#{$p_def_vect},igaus) * gpcar(#{$p_def_vect},idime,inode,igaus)
                  do jnode = 1,pnode
                     do jdime = 1,ndime
                        jdofv                       = (jnode-1)*ndime + jdime
                        elauu(#{$p_def_vect},idofv,jdofv) = elauu(#{$p_def_vect},idofv,jdofv) + fact1(#{$p_def_vect}) * gpcar(#{$p_def_vect},jdime,jnode,igaus)
                     end do
                  end do
               end if
            end do
         end do
      end do
   end if
EOF
   nests.push <<EOF
   if( kfl_lumped == 1 ) then 
      if( ndime == 2 ) then
         stop
      else
         do igaus = 1,pgaus
            gpveo(#{$p_def_vect},1:3) = 0.0_rp
            do inode = 1,pnode
               do idime = 1,ndime
                  gpveo(#{$p_def_vect},idime) = gpveo(#{$p_def_vect},idime) + elvel(#{$p_def_vect},idime,inode,2) * gpsha(#{$p_def_vect},inode,igaus)
               end do
            end do
            do inode = 1,pnode
               idof1                       = 3*inode-2
               idof2                       = 3*inode-1
               idof3                       = 3*inode
               fact0(#{$p_def_vect})             = gpvol(#{$p_def_vect},igaus) * gpden(#{$p_def_vect},igaus) * gpsha(#{$p_def_vect},inode,igaus) * dtinv_mod(#{$p_def_vect})
               elauu(#{$p_def_vect},idof1,idof1) = elauu(#{$p_def_vect},idof1,idof1) + fact0(#{$p_def_vect})
               elauu(#{$p_def_vect},idof2,idof2) = elauu(#{$p_def_vect},idof2,idof2) + fact0(#{$p_def_vect})
               elauu(#{$p_def_vect},idof3,idof3) = elauu(#{$p_def_vect},idof3,idof3) + fact0(#{$p_def_vect})
               do idime = 1,ndime
                  elrbu(#{$p_def_vect},idime,inode) = elrbu(#{$p_def_vect},idime,inode) - fact0(#{$p_def_vect}) * gpveo(#{$p_def_vect},idime)
                  elrbu(#{$p_def_vect},idime,inode) = elrbu(#{$p_def_vect},idime,inode) + fact0(#{$p_def_vect}) * elvel(#{$p_def_vect},idime,inode,2)
               end do
               do jnode = 1,pnode 
                  jdof1                       = 3*jnode-2
                  jdof2                       = 3*jnode-1
                  jdof3                       = 3*jnode
                  elauu(#{$p_def_vect},idof1,jdof1) = elauu(#{$p_def_vect},idof1,jdof1) - fact0(#{$p_def_vect}) * gpsha(#{$p_def_vect},jnode,igaus) 
                  elauu(#{$p_def_vect},idof2,jdof2) = elauu(#{$p_def_vect},idof2,jdof2) - fact0(#{$p_def_vect}) * gpsha(#{$p_def_vect},jnode,igaus) 
                  elauu(#{$p_def_vect},idof3,jdof3) = elauu(#{$p_def_vect},idof3,jdof3) - fact0(#{$p_def_vect}) * gpsha(#{$p_def_vect},jnode,igaus) 
               end do
            end do
         end do
      end if
   
   else if( kfl_lumped == 2 ) then 
      do igaus = 1,pgaus
         fact0(#{$p_def_vect}) = gpvol(#{$p_def_vect},igaus) * gpden(#{$p_def_vect},igaus) * dtinv_mod(#{$p_def_vect})
         do inode = 1, pnode
            fact1(#{$p_def_vect}) = fact0(#{$p_def_vect}) * gpsha(#{$p_def_vect},inode,igaus)
            do idime = 1,ndime
               idof1                       = (inode-1) * ndime + idime
               elauu(#{$p_def_vect},idof1,idof1) = elauu(#{$p_def_vect},idof1,idof1) + fact1(#{$p_def_vect})
               elrbu(#{$p_def_vect},idime,inode) = elrbu(#{$p_def_vect},idime,inode) + fact1(#{$p_def_vect}) * elvel(#{$p_def_vect},idime,inode,2)
            end do
         end do
      end do
   end if
EOF
   nests.push <<EOF
   if( ndime == 2 ) then
      do igaus = 1,pgaus
         do inode = 1,pnode
            idof1 = 2*inode-1
            idof2 = 2*inode
            do jnode = 1,pnode
               fact0(#{$p_def_vect})             = gpvol(#{$p_def_vect},igaus)       * gpsha(#{$p_def_vect},jnode,igaus) 
               fact1(#{$p_def_vect})             = fact0(#{$p_def_vect})             * gpcar(#{$p_def_vect},1,inode,igaus)
               fact2(#{$p_def_vect})             = fact0(#{$p_def_vect})             * gpcar(#{$p_def_vect},2,inode,igaus)
               elapu(#{$p_def_vect},jnode,idof1) = elapu(#{$p_def_vect},jnode,idof1) + fact1(#{$p_def_vect})
               elapu(#{$p_def_vect},jnode,idof2) = elapu(#{$p_def_vect},jnode,idof2) + fact2(#{$p_def_vect})
               elaup(#{$p_def_vect},idof1,jnode) = elaup(#{$p_def_vect},idof1,jnode) - fact1(#{$p_def_vect})
               elaup(#{$p_def_vect},idof2,jnode) = elaup(#{$p_def_vect},idof2,jnode) - fact2(#{$p_def_vect})
            end do
         end do
      end do
   else
      do igaus = 1,pgaus
         do inode = 1,pnode
            idof1 = 3*inode-2
            idof2 = 3*inode-1
            idof3 = 3*inode
            do jnode = 1,pnode
               fact0(#{$p_def_vect})             = gpvol(#{$p_def_vect},igaus)       * gpsha(#{$p_def_vect},jnode,igaus) 
               fact1(#{$p_def_vect})             = fact0(#{$p_def_vect})             * gpcar(#{$p_def_vect},1,inode,igaus)
               fact2(#{$p_def_vect})             = fact0(#{$p_def_vect})             * gpcar(#{$p_def_vect},2,inode,igaus)
               fact3(#{$p_def_vect})             = fact0(#{$p_def_vect})             * gpcar(#{$p_def_vect},3,inode,igaus)
               elapu(#{$p_def_vect},jnode,idof1) = elapu(#{$p_def_vect},jnode,idof1) + fact1(#{$p_def_vect})
               elapu(#{$p_def_vect},jnode,idof2) = elapu(#{$p_def_vect},jnode,idof2) + fact2(#{$p_def_vect})
               elapu(#{$p_def_vect},jnode,idof3) = elapu(#{$p_def_vect},jnode,idof3) + fact3(#{$p_def_vect}) 
               elaup(#{$p_def_vect},idof1,jnode) = elaup(#{$p_def_vect},idof1,jnode) - fact1(#{$p_def_vect})
               elaup(#{$p_def_vect},idof2,jnode) = elaup(#{$p_def_vect},idof2,jnode) - fact2(#{$p_def_vect})
               elaup(#{$p_def_vect},idof3,jnode) = elaup(#{$p_def_vect},idof3,jnode) - fact3(#{$p_def_vect})
            end do
         end do
      end do
   end if
EOF
   nests.push <<EOF
   if( kfl_stabi_nsi /= -1 ) then
      do igaus = 1,pgaus
         do inode = 1,pnode
            do jnode = inode+1,pnode
               fact1(#{$p_def_vect})             = gpsp1_p(#{$p_def_vect},igaus) * wgrgr(#{$p_def_vect},jnode,inode,igaus) * gpvol(#{$p_def_vect},igaus)
               elapp(#{$p_def_vect},jnode,inode) = elapp(#{$p_def_vect},jnode,inode) + fact1(#{$p_def_vect})
               elapp(#{$p_def_vect},inode,jnode) = elapp(#{$p_def_vect},inode,jnode) + fact1(#{$p_def_vect})
            end do
            fact1(#{$p_def_vect})             = gpsp1_p(#{$p_def_vect},igaus) * wgrgr(#{$p_def_vect},inode,inode,igaus) * gpvol(#{$p_def_vect},igaus)
            elapp(#{$p_def_vect},inode,inode) = elapp(#{$p_def_vect},inode,inode) + fact1(#{$p_def_vect})
         end do
      end do
    end if
EOF
   nests.push <<EOF
   do igaus = 1,pgaus
     fact1(#{$p_def_vect}) = penal_nsi * gpvol(#{$p_def_vect},igaus)
     do inode = 1,pnode
       elapp(#{$p_def_vect},inode,inode) = elapp(#{$p_def_vect},inode,inode) + fact1(#{$p_def_vect}) * gpsha(#{$p_def_vect},inode,igaus)
       elrbp(#{$p_def_vect},inode)       = elrbp(#{$p_def_vect},inode)       + fact1(#{$p_def_vect}) * gpsha(#{$p_def_vect},inode,igaus) * elpre(#{$p_def_vect},inode,1) 
     end do
   end do
EOF
   nests.push <<EOF
   if( kfl_stabi_nsi == -1 ) then
      gpvep(#{$p_def_vect},:,:) = 0.0_rp 
   else if( kfl_limit_nsi == -1 ) then
      gpvep(#{$p_def_vect},:,:) = 0.0_rp
   else if( kfl_limit_nsi > 0 ) then
      do igaus = 1,pgaus
         c1(#{$p_def_vect}) = 0.0_rp
         c2(#{$p_def_vect}) = 0.0_rp
         c3(#{$p_def_vect}) = 0.0_rp
         do idime = 1,ndime
            c4(#{$p_def_vect}) = 0.0_rp
            do inode = 1,pnode
               c4(#{$p_def_vect}) = c4(#{$p_def_vect}) + agrau(#{$p_def_vect},inode,igaus) * elvel(#{$p_def_vect},idime,inode,1)
            end do
            c4(#{$p_def_vect}) = gpsp1(#{$p_def_vect},igaus) * c4(#{$p_def_vect})
            c1(#{$p_def_vect}) = c1(#{$p_def_vect}) + ( gpvep(#{$p_def_vect},idime,igaus) - c4(#{$p_def_vect}) )**2
            c3(#{$p_def_vect}) = c3(#{$p_def_vect}) + gpvep(#{$p_def_vect},idime,igaus) * gpvep(#{$p_def_vect},idime,igaus)
            c2(#{$p_def_vect}) = c2(#{$p_def_vect}) + c4(#{$p_def_vect}) * c4(#{$p_def_vect})
         end do
         c3(#{$p_def_vect})   = sqrt( c2(#{$p_def_vect}) ) + sqrt( c3(#{$p_def_vect}) )
         c1(#{$p_def_vect})   = sqrt( c1(#{$p_def_vect}) )
         beta(#{$p_def_vect}) = c1(#{$p_def_vect}) / ( c3(#{$p_def_vect}) + epsilon(1.0_rp) )
         if( kfl_limit_nsi == 1 ) then
            alpha(#{$p_def_vect}) = min(1.0_rp,2.0_rp*(1.0_rp-beta(#{$p_def_vect})))
         else if( kfl_limit_nsi == 2 ) then
            alpha(#{$p_def_vect}) = 0.5_rp*(tanh(20.0_rp*(beta(#{$p_def_vect})-0.8_rp))+1.0_rp)
         end if
         do idime = 1,ndime
            gpvep(#{$p_def_vect},idime,igaus) = alpha(#{$p_def_vect}) * gpvep(#{$p_def_vect},idime,igaus)
         end do
      end do
   end if
EOF
   nests.push <<EOF
   if( kfl_stabi_nsi == -1 ) then
      gpgrp(#{$p_def_vect},:,:) = 0.0_rp
   else
      do igaus = 1,pgaus
         do idime = 1,ndime
            gpgrp(#{$p_def_vect},idime,igaus) = gpgrp(#{$p_def_vect},idime,igaus) + gpsp1_p(#{$p_def_vect},igaus) * gprhs(#{$p_def_vect},idime,igaus)
         end do
      end do
      if( kfl_sgsti_nsi == 1 ) then
         do igaus = 1,pgaus 
            fact1(#{$p_def_vect})    = gpden(#{$p_def_vect},igaus) * dtsgs(#{$p_def_vect}) * gpsp1_v(#{$p_def_vect},igaus)
            fact1_p (#{$p_def_vect}) = gpden(#{$p_def_vect},igaus) * dtsgs(#{$p_def_vect}) * gpsp1_p(#{$p_def_vect},igaus)
            do idime = 1,ndime
               gpvep(#{$p_def_vect},idime,igaus) = gpvep(#{$p_def_vect},idime,igaus) + fact1(#{$p_def_vect})   * gpsgs(#{$p_def_vect},idime,igaus,2)
               gpgrp(#{$p_def_vect},idime,igaus) = gpgrp(#{$p_def_vect},idime,igaus) + fact1_p(#{$p_def_vect}) * gpsgs(#{$p_def_vect},idime,igaus,2)
            end do
         end do
      end if
   end if
EOF
   nests.push <<EOF
   do igaus = 1,pgaus
      fact4(#{$p_def_vect}) = gpden(#{$p_def_vect},igaus) * dtinv_mod(#{$p_def_vect})
      do itime = 2,nbdfp_nsi
         do idime = 1,ndime
            gprhs(#{$p_def_vect},idime,igaus) = gprhs(#{$p_def_vect},idime,igaus) - pabdf_nsi(itime) * fact4(#{$p_def_vect}) * gpvel(#{$p_def_vect},idime,igaus,itime)
         end do
      end do
      do inode = 1,pnode
         fact1(#{$p_def_vect}) = gpvol(#{$p_def_vect},igaus) * gpsha(#{$p_def_vect},inode,igaus)  ! ( f + rho*u^n/dt , v )
         fact3(#{$p_def_vect}) = gpvol(#{$p_def_vect},igaus) * agrau(#{$p_def_vect},inode,igaus)  ! ( rho * a.grad(v) , P1' ) 
         do idime = 1,ndime
            elrbu(#{$p_def_vect},idime,inode) = elrbu(#{$p_def_vect},idime,inode) + fact1(#{$p_def_vect}) * gprhs(#{$p_def_vect},idime,igaus) &
                 &                                                    + fact3(#{$p_def_vect}) * gpvep(#{$p_def_vect},idime,igaus)              
         end do
         elrbp(#{$p_def_vect},inode) = elrbp(#{$p_def_vect},inode) + gpvol(#{$p_def_vect},igaus) * gpsha(#{$p_def_vect},inode,igaus) * gprhc(#{$p_def_vect},igaus)                ! ( rhs, q )
         do idime = 1,ndime
            elrbp(#{$p_def_vect},inode) = elrbp(#{$p_def_vect},inode) + gpvol(#{$p_def_vect},igaus) * gpcar(#{$p_def_vect},idime,inode,igaus) * gpgrp(#{$p_def_vect},idime,igaus) ! ( P2' , grad(q) ) 
         end do
      end do
   end do
EOF
   nests.push <<EOF
   if( maxval(pbubl) == 1 ) then
      if( kfl_stabi_nsi /= -1 ) then
         write(6,*) 'BUBBLE NOT CODED FOR SPLIT OSS'
         stop
      end if
   
      elauq = 0.0_rp
      elapq = 0.0_rp
      elaqu = 0.0_rp
      elaqp = 0.0_rp
      elaqq = 0.0_rp
      elrbq = 0.0_rp
   
      if( kfl_press_nsi == 1 ) then
         do igaus = 1,pgaus
            fact1(#{$p_def_vect}) = gpvol(#{$p_def_vect},igaus) * gpsha_bub(#{$p_def_vect},igaus)
            do inode = 1,pnode
               do idime = 1,ndime
                  idofv = (inode-1)*ndime + idime
                  elauq(#{$p_def_vect},idofv,1) = elauq(#{$p_def_vect},idofv,1) - fact1(#{$p_def_vect}) * gpcar(#{$p_def_vect},idime,inode,igaus)
                  elaqu(#{$p_def_vect},1,idofv) = elaqu(#{$p_def_vect},1,idofv) + fact1(#{$p_def_vect}) * gpcar(#{$p_def_vect},idime,inode,igaus) 
               end do
            end do
         end do
      else
         do igaus = 1,pgaus
            fact1(#{$p_def_vect}) = gpvol(#{$p_def_vect},igaus) * gpsha_bub(#{$p_def_vect},igaus)
            do inode = 1,pnode
               do idime = 1,ndime
                  idofv = (inode-1)*ndime + idime
                  elauq(#{$p_def_vect},idofv,1) = elauq(#{$p_def_vect},idofv,1) + gpvol(#{$p_def_vect},igaus) * gpsha(#{$p_def_vect},inode,igaus) * gpcar_bub(#{$p_def_vect},idime,igaus)
                  elaqu(#{$p_def_vect},1,idofv) = elaqu(#{$p_def_vect},1,idofv) + fact1(#{$p_def_vect}) * gpcar(#{$p_def_vect},idime,inode,igaus) 
               end do
            end do
         end do
      end if
         
      do igaus = 1,pgaus
         elaqq(#{$p_def_vect},1,1) = elaqq(#{$p_def_vect},1,1) + gpvol(#{$p_def_vect},igaus) * gpsha_bub(#{$p_def_vect},igaus) * penal_nsi
         elrbq(#{$p_def_vect},1)   = elrbq(#{$p_def_vect},1)   + gpvol(#{$p_def_vect},igaus) * gpsha_bub(#{$p_def_vect},igaus) * penal_nsi * elbub(#{$p_def_vect}) 
         elrbq(#{$p_def_vect},1)   = elrbq(#{$p_def_vect},1)   + gpvol(#{$p_def_vect},igaus) * gpsha_bub(#{$p_def_vect},igaus) * gprhc(#{$p_def_vect},igaus) 
      end do
   
   end if
EOF
   set_lang(FORTRAN)
   @kernel = CKernel::new(:includes => "immintrin.h")
   @kernel.procedure = declare_procedure("nsi_element_assembly_split_oss_ref")
   get_output.print macros
   get_output.print generate_mocks
   get_output.print generate_ref_declaration
   get_output.print generate_ref_initialization

   @opts[:nests].each{|n|
     get_output.print nests[n-1]
   }

   get_output.print "end subroutine nsi_element_assembly_split_oss_ref"
   return @kernel
 end  
end
