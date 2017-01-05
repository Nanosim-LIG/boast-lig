require './KASGS.rb'

class KASGSRef < KASGS
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

  def generate
    nests = []
    
    macro =<<EOF
#define VECTOR_SIZE #{@opts[:vector_length]}
#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
EOF

    subroutine_prototype =<<EOF
    subroutine nsi_element_assembly_asgs_oss(ndime,     &
       pnode,pgaus,gpden,gpvis,gppor,gpgvi,gpsp1,       &
       gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
       gpadv,gpvep,gprhs,gprhc,rmomu,rcont,elpre,elbub, &
       elauu,elaup,elapp,elapu,elrbu,elrbp,rmom2,gpst1, &
       gpgve,gprh2,gprhs_sgs,elvel,ellum,dtinv_loc,     &
       pbubl,gpsha_bub,gpcar_bub,gppre,elauq,elapq,     &
       elaqu,elaqp,elaqq,elrbq,&
       ! Original global variables
       kfl_lumped,&
       mnode,ntens,&
       kfl_duatss,&
       fact_duatss,&
       kfl_stabi_nsi,&
       fvins_nsi,fcons_nsi,bemol_nsi,kfl_regim_nsi,&
       fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi,&
       kfl_p1ve2_nsi,kfl_linea_nsi,pabdf_nsi,&
       kfl_confi_nsi,nbdfp_nsi,kfl_sgsti_nsi,&
       kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi,&
       penal_nsi,&
       kfl_bubbl_nsi)

    #{gen_def_parameters}

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
    real(rp),    intent(in)    :: pabdf_nsi(*)
    integer(ip), intent(in)    :: kfl_confi_nsi
    integer(ip), intent(in)    :: nbdfp_nsi
    integer(ip), intent(in)    :: kfl_sgsti_nsi
    integer(ip), intent(in)    :: kfl_nota1_nsi
    integer(ip), intent(in)    :: kfl_limit_nsi
    integer(ip), intent(in)    :: kfl_penal_nsi
    real(rp),    intent(in)    :: penal_nsi
    integer(ip), intent(in)    :: kfl_bubbl_nsi

    integer(ip), intent(in)    :: pnode,pgaus
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgvi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gplap(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gprh2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: rmomu(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    ! Element matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: gpst1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(out)   :: ellum(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)                    ! Ne
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)              ! dNe/dxi
    real(rp),    intent(in)    :: gppre(VECTOR_SIZE,pgaus) 
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)
    ! Indices and local variables
    integer(ip)                :: inode,jnode,idofv
    integer(ip)                :: idof1,jdime,ivect
    integer(ip)                :: igaus,idime,jdofv,kdime
    real(rp)                   :: xvis2,xvisc,one_rp
    ! Vectorized local arrays
    real(rp)                   :: p1ve2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp)                   :: p1vec(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: p2vec(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp)                   :: p2sca(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)                   :: wgrvi(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: fact0(VECTOR_SIZE),fact1(VECTOR_SIZE)
    real(rp)                   :: fact2(VECTOR_SIZE),fact3(VECTOR_SIZE)
    real(rp)                   :: fact4(VECTOR_SIZE),fact5(VECTOR_SIZE)
    real(rp)                   :: fact6(VECTOR_SIZE),fact7(VECTOR_SIZE)
    real(rp)                   :: fact8(VECTOR_SIZE),factz(VECTOR_SIZE)
    real(rp)                   :: factx(VECTOR_SIZE),facty(VECTOR_SIZE)
    real(rp)                   :: facx1(VECTOR_SIZE),facy1(VECTOR_SIZE)
    real(rp)                   :: ugraN(VECTOR_SIZE)
    real(rp)                   :: gramugraN(VECTOR_SIZE)
    real(rp)                   :: penal(VECTOR_SIZE)
    real(rp)                   :: gprhh(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: taupr(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpveo(VECTOR_SIZE,ndime)
    real(rp)                   :: gpcar1ji(VECTOR_SIZE)
    real(rp)                   :: gpcar2ji(VECTOR_SIZE)
    real(rp)                   :: gpcar3ji(VECTOR_SIZE)
    real(rp)                   :: p2sca_bub(VECTOR_SIZE,pgaus)
    real(rp)                   :: p2vec_bub(VECTOR_SIZE,ndime,pgaus)

    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu11,elauu21,elauu31
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu12,elauu22,elauu32
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu13,elauu23,elauu33
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4))       :: factvec1,factvec2,factvec3,factvec4
EOF

    locals =<<EOF
#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !----------------------------------------------------------------------
    !
    ! Local variables
    !
    !----------------------------------------------------------------------
    !
    ! Viscous term factor
    !
    xvis2 = 0.0_rp
    if( fvins_nsi > 0.9_rp ) then
       xvisc = 1.0_rp
       if( fvins_nsi > 1.9_rp ) xvis2 = 2.0_rp/3.0_rp
    else
       xvisc = 0.0_rp
    end if
    !
    ! Do not stabilize convection, reaction, coriolis
    !
    if( kfl_nota1_nsi == 1 ) then
       one_rp = 0.0_rp 
    else
       one_rp = 1.0_rp
    end if
EOF

    init =<<EOF
    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elaup(DEF_VECT,:,:)   = 0.0_rp
    elapu(DEF_VECT,:,:)   = 0.0_rp
    elapp(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    elrbp(DEF_VECT,:)     = 0.0_rp
 
    elauu11(DEF_VECT,:,:) = 0.0_rp
    elauu21(DEF_VECT,:,:) = 0.0_rp
    elauu12(DEF_VECT,:,:) = 0.0_rp
    elauu22(DEF_VECT,:,:) = 0.0_rp
    elauu31(DEF_VECT,:,:) = 0.0_rp
    elauu32(DEF_VECT,:,:) = 0.0_rp
    elauu13(DEF_VECT,:,:) = 0.0_rp
    elauu23(DEF_VECT,:,:) = 0.0_rp
    elauu33(DEF_VECT,:,:) = 0.0_rp
    !
    ! Incompressible vs low Mach
    !
    if( kfl_regim_nsi == 3 ) then
       taupr(DEF_VECT,1:pgaus) = 1.0_rp
    else
       taupr(DEF_VECT,1:pgaus) = 1.0_rp / max(gpst1(DEF_VECT,1:pgaus),zeror)
    end if

EOF

    nests[0] = <<EOF
    do igaus = 1,pgaus

       fact1(DEF_VECT) = gpsp1(DEF_VECT,igaus) * one_rp * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       fact2(DEF_VECT) = (gptt1(DEF_VECT,igaus)-gpsp1(DEF_VECT,igaus) * one_rp * gppor(DEF_VECT,igaus)) * gpvol(DEF_VECT,igaus)
       fact3(DEF_VECT) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       fact4(DEF_VECT) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode

          ugraN(DEF_VECT)     = 0.0_rp
          gramugraN(DEF_VECT) = 0.0_rp

          do idime = 1,ndime
             ugraN(DEF_VECT)                   = ugraN(DEF_VECT)     + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             gramugraN(DEF_VECT)               = gramugraN(DEF_VECT) + gpgvi(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             p2vec(DEF_VECT,idime,inode,igaus) = fact3(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
          end do

          p1vec(DEF_VECT,inode,igaus)          = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) + fact1(DEF_VECT) * ugraN(DEF_VECT)
          p2sca(DEF_VECT,inode,igaus)          = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
          wgrvi(DEF_VECT,inode,igaus)          = gramugraN(DEF_VECT) + gpvis(DEF_VECT,igaus) * gplap(DEF_VECT,inode,igaus)

          do jnode = 1,pnode           
             wgrgr(DEF_VECT,inode,jnode,igaus) = 0.0_rp
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do

       end do
    end do
EOF

    nests[1] = <<EOF
    if( kfl_p1ve2_nsi /= 0 ) then
       !
       ! tau1'*2*rho*(w x v)
       ! x-equation: v = (v,0,0) => w x v = (    0 , wz*v , -wy*v)     
       ! y-equation: v = (0,v,0) => w x v = (-wz*v ,    0 ,  wx*v)     
       ! z-equation: v = (0,0,v) => w x v = ( wy*v ,-wx*v ,     0)     
       ! 
       p1ve2 = 0.0_rp
       if( ndime == 2 ) then
          do igaus = 1,pgaus           
             fact3(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * fvela_nsi(3)
             factz(DEF_VECT) = gpsp1(DEF_VECT,igaus)* fact3(DEF_VECT) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode 
                fact1(DEF_VECT)                 =   factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                p1ve2(DEF_VECT,1,2,inode,igaus) =   fact1(DEF_VECT)
                p1ve2(DEF_VECT,2,1,inode,igaus) = - fact1(DEF_VECT)
             end do
          end do
       else
          do igaus = 1,pgaus   
             fact8(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * gpsp1(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
             factx(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(1)
             facty(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(2)
             factz(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(3)
             do inode = 1,pnode
                p1ve2(DEF_VECT,1,2,inode,igaus) =  factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                p1ve2(DEF_VECT,1,3,inode,igaus) = -facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                p1ve2(DEF_VECT,2,1,inode,igaus) = -factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                p1ve2(DEF_VECT,2,3,inode,igaus) =  factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                p1ve2(DEF_VECT,3,1,inode,igaus) =  facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation
                p1ve2(DEF_VECT,3,2,inode,igaus) = -factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation    
             end do
          end do
       end if
    end if
EOF

    nests[2] =<<EOF
    if( ndime == 2 ) then

       do igaus = 1,pgaus

          fact0(DEF_VECT) = ( gpsp2(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)                                      ! (1) + (5)
          fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

          do inode = 1,pnode
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) 
             fact2(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)
             fact5(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)

             do jnode = 1,pnode              
                fact4(DEF_VECT)               =   p1vec(DEF_VECT,inode,igaus) * rmomu(DEF_VECT,jnode,igaus) &                                          ! (2):       ( rmomu(DEF_VECT,u) , p1vec(DEF_VECT,v) ) 
                     &                          + fact5(DEF_VECT) * wgrvi(DEF_VECT,jnode,igaus)              &                                         ! (3) + (6): ( mu d^2ui/dxj^2 + dmu/dxj dui/dxj, vi )
                     &                          + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)                                                  ! (7):       ( mu dui/dxj , dvi/dxj ) 

                elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * rcont(DEF_VECT,1,jnode,igaus) + fact4(DEF_VECT)      ! Auu_xx
                elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * rcont(DEF_VECT,1,jnode,igaus)                        ! Auu_yx
                elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * rcont(DEF_VECT,2,jnode,igaus)                        ! Auu_xy
                elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * rcont(DEF_VECT,2,jnode,igaus) + fact4(DEF_VECT)      ! Auu_yy
             end do

          end do

       end do

    else

       do igaus = 1,pgaus    

          fact0(DEF_VECT) = ( gpsp2(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)                                          ! (1) + (5)
          fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

          factvec1(DEF_VECT,1:pnode) = gpcar(DEF_VECT,1,1:pnode,igaus)
          factvec2(DEF_VECT,1:pnode) = gpcar(DEF_VECT,2,1:pnode,igaus)
          factvec3(DEF_VECT,1:pnode) = gpcar(DEF_VECT,3,1:pnode,igaus)

          do jnode = 1,pnode              
             fact1(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,1,jnode,igaus) 
             fact2(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,2,jnode,igaus)
             fact3(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,3,jnode,igaus)
             fact5(DEF_VECT) = wgrvi(DEF_VECT,jnode,igaus)*gpvol(DEF_VECT,igaus)

             do inode = 1,pnode
                factvec4(DEF_VECT,inode) = p1vec(DEF_VECT,inode,igaus) * rmomu(DEF_VECT,jnode,igaus) + fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) &
                     + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)
             end do

             do inode = 1,pnode
                elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec1(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_xx
                elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec1(DEF_VECT,inode)                             ! Auu_xy
                elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec1(DEF_VECT,inode)                             ! Auu_xz

                elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec2(DEF_VECT,inode)                             ! Auu_yx
                elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec2(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_yy
                elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec2(DEF_VECT,inode)                             ! Auu_yz

                elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec3(DEF_VECT,inode)                             ! Auu_zx
                elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec3(DEF_VECT,inode)                             ! Auu_zy
                elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec3(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_zz
             end do

          end do

       end do

    end if

EOF

    nests[3] =<<EOF
    if( fvins_nsi > 0.9_rp ) then

       if( ndime == 2 ) then

          do igaus = 1,pgaus

             fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
             fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

             do inode = 1,pnode
                fact0(DEF_VECT) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) *(1.0_rp-xvis2)
                fact1(DEF_VECT) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                fact2(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                facx1(DEF_VECT) = fact1(DEF_VECT) * xvis2
                facy1(DEF_VECT) = fact2(DEF_VECT) * xvis2    
                do jnode = 1,pnode              
                   !
                   !                                         (4)  ( mu duj/dxi , dvi/dxj )
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                   ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                   ! Auu_xy
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                   ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                   ! Auu_yy
                   !
                   !                                         (10) - 2/3 * B * ( dmu/dxi (div u) , vi )
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - facx1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                   ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) - facx1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                   ! Auu_xy
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) - facy1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                   ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - facy1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                   ! Auu_yy
                   !
                   !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                   !                                               + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,jnode,igaus) ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_xy
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,jnode,igaus) ! Auu_yy
                   !
                   !                                        

                end do

             end do
          end do

       else

          do igaus = 1,pgaus

             fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
             fact6(DEF_VECT) = gpgvi(DEF_VECT,3,igaus) * gpvol(DEF_VECT,igaus)
             fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
             !          factx(DEF_VECT) = xvis2 * fact4(DEF_VECT)
             !          facty(DEF_VECT) = xvis2 * fact5(DEF_VECT)
             !          factz(DEF_VECT) = xvis2 * fact6(DEF_VECT)

             !yor! [Originally line 451] reverted loops order, refactored the factors, with an effort to vectorize each instruction 
             !
             !                                (4)  ( mu duj/dxi , dvi/dxj )
             factvec1(DEF_VECT,1:pnode) = gpcar(DEF_VECT,1,1:pnode,igaus)
             factvec2(DEF_VECT,1:pnode) = gpcar(DEF_VECT,2,1:pnode,igaus)
             factvec3(DEF_VECT,1:pnode) = gpcar(DEF_VECT,3,1:pnode,igaus)

             do jnode = 1,pnode
                gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)
                !
                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xx
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yx
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zx

                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xy
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yy
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zy

                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xz
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yz
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zz
                end do
             end do
             !
             !                     (10) - 2/3 * B * ( dmu/dxi (div u) , vi )

             do inode = 1,pnode
                factvec1(DEF_VECT,inode) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facx1) this copy is for pure esthaetical purpose
                factvec2(DEF_VECT,inode) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facy1)
                factvec3(DEF_VECT,inode) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facz1)  
             end do

             do jnode = 1,pnode
                gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)

                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_xx 
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_xy
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_xz

                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_yy
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_yz

                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_zx
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_zy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_zz
                end do
             end do

             !                                                                   
             !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
             !                                                + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi ) 
             do inode = 1,pnode
                factvec4(DEF_VECT,inode) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)*(1.0_rp-xvis2) !yor!(fact0(DEF_VECT))
                factvec1(DEF_VECT,inode) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact1(DEF_VECT))
                factvec2(DEF_VECT,inode) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact2(DEF_VECT))
                factvec3(DEF_VECT,inode) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact3(DEF_VECT))
             end do

             do jnode = 1,pnode
                gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)

                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,1,jnode,igaus) ! Auu_xx
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,4,jnode,igaus) ! Auu_yx
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,5,jnode,igaus) ! Auu_zx

                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,4,jnode,igaus) ! Auu_xy
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,2,jnode,igaus) ! Auu_yy
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,6,jnode,igaus) ! Auu_zy

                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,5,jnode,igaus) ! Auu_xz
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,6,jnode,igaus) ! Auu_yz
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_zz
                end do
             end do
          end do
       end if
    end if
EOF

    nests[4] =<<EOF
    if( fcons_nsi > 0.1_rp ) then   
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                do jnode = 1,pnode     
                   fact1(DEF_VECT) = 0.0_rp
                   fact2(DEF_VECT) = 0.0_rp
                   do idime = 1,2
                      fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                      fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) 
                   end do
                   fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                   fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT) 

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_yy
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                do jnode = 1,pnode
                   fact1(DEF_VECT) = 0.0_rp
                   fact2(DEF_VECT) = 0.0_rp
                   do idime = 1,3
                      fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                      fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) 
                   end do
                   fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                   fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT) 

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_yy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_zz
                end do
             end do
          end do
       end if
    end if
EOF
    
    nests[5] =<<EOF
    if( kfl_rmom2_nsi /= 0 ) then
       !
       ! p1vec * rmom2(DEF_VECT,
       !
       !yor! [Originally loop on line 577] Reordered loops to solve indirections problems 
       if (ndime == 2) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,1,1,jnode,igaus)
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,2,1,jnode,igaus)
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,1,2,jnode,igaus)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,2,2,jnode,igaus)            
                end do
             end do
          end do
       else !ndime==3
          do igaus = 1,pgaus
             factvec1(DEF_VECT,1:pnode) = p1vec(DEF_VECT,1:pnode,igaus)
             do jnode = 1,pnode
                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,1,jnode,igaus)
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,1,jnode,igaus)
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,1,jnode,igaus)
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,2,jnode,igaus)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,2,jnode,igaus)
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,2,jnode,igaus)
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,3,jnode,igaus)
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,3,jnode,igaus)
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,3,jnode,igaus)
                end do
             end do
          end do
       end if

       if( kfl_p1ve2_nsi /= 0 ) then
          !
          ! p1ve2 * rmomu
          ! p1ve2 * rmom2(DEF_VECT,
          !
          !yor! optim : [originally line 597] reverted loops, refactored 
          if(ndime==2)then   
             do igaus=1,pgaus
                do jnode=1,pnode
                   do inode=1,pnode
                      elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      do kdime = 1,ndime
                         elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                      end do
                   end do
                end do
             end do
          else ! ndime==3
             do igaus=1,pgaus
                do jnode=1,pnode
                   do inode=1,pnode
                      elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu31(DEF_VECT,inode,jnode)=elauu31(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu32(DEF_VECT,inode,jnode)=elauu32(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu13(DEF_VECT,inode,jnode)=elauu13(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu23(DEF_VECT,inode,jnode)=elauu23(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      elauu33(DEF_VECT,inode,jnode)=elauu33(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                      do kdime = 1,ndime
                         elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu31(DEF_VECT,inode,jnode)=elauu31(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                         elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu32(DEF_VECT,inode,jnode)=elauu32(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         elauu13(DEF_VECT,inode,jnode)=elauu13(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                         elauu23(DEF_VECT,inode,jnode)=elauu23(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                         elauu33(DEF_VECT,inode,jnode)=elauu33(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                      end do
                   end do
                end do
             end do
          end if
       end if
    end if
EOF

    nests[6] =<<EOF
if( abs(bemol_nsi).gt.1.0e-9_rp ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) )
             do inode = 1,pnode
                do jnode = 1,pnode              
                   !
                   ! fact2(DEF_VECT) = - b * ( v , rho * (a.grad) u )
                   ! fact3(DEF_VECT) = - b * ( u , rho * (a.grad) v ) 
                   ! fact4(DEF_VECT) = - b * ( rho * div(a) u , v ) 
                   !
                   fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) ) * gpsha(DEF_VECT,jnode,igaus)
                   fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) ) * gpsha(DEF_VECT,inode,igaus)
                   fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus) 
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) + gpgve(DEF_VECT,3,3,igaus))
             do inode = 1,pnode
                do jnode = 1,pnode
                   !
                   ! fact2(DEF_VECT) = - b * ( v , rho * (a.grad) u )
                   ! fact3(DEF_VECT) = - b * ( u , rho * (a.grad) v ) 
                   ! fact4(DEF_VECT) = - b * ( rho * div(a) u , v ) 
                   !
                   fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,inode,igaus)) * gpsha(DEF_VECT,jnode,igaus)
                   fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,jnode,igaus)) * gpsha(DEF_VECT,inode,igaus)
                   fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                end do
             end do
          end do
       end if

    end if
EOF

    nests[7] =<<EOF
    if( kfl_linea_nsi == 2 ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                do jnode = 1,pnode
                   fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   ! Auu_1i
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                   ! Auu_2i
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                do jnode = 1,pnode
                   fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   ! Auu_1i
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,1,igaus) ! rho * uz * dux^i-1/dz
                   ! Auu_2i
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,2,igaus) ! rho * uz * duy^i-1/dz
                   ! Auu_3i
                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,3,igaus) ! rho * ux * duz^i-1/dx
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,3,igaus) ! rho * uy * duz^i-1/dy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,3,igaus) ! rho * uz * duz^i-1/dz

                end do
             end do
          end do
       end if
    end if
EOF

    nests[8] =<<EOF
    if( kfl_lumped == 1 ) then 
       !
       ! Remove Galerkin term and add lumped term 
       ! 
       if( ndime == 2 ) then
          stop
       else
          do igaus =1, pgaus
             gpveo(DEF_VECT,1:3) = 0.0_rp
             do inode = 1,pnode
                gpveo(DEF_VECT,1) = gpveo(DEF_VECT,1) + elvel(DEF_VECT,1, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                gpveo(DEF_VECT,2) = gpveo(DEF_VECT,2) + elvel(DEF_VECT,2, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                gpveo(DEF_VECT,3) = gpveo(DEF_VECT,3) + elvel(DEF_VECT,3, inode, 2)* gpsha(DEF_VECT,inode, igaus)
             end do
             do inode =1, pnode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)*dtinv_loc(DEF_VECT)
                elauu11(DEF_VECT,inode,inode) =  elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                elauu22(DEF_VECT,inode,inode) =  elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                elauu33(DEF_VECT,inode,inode) =  elauu33(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode)     =  elrbu(DEF_VECT,idime,inode)   - fact0(DEF_VECT)*gpveo(DEF_VECT,idime)
                   elrbu(DEF_VECT,idime,inode)     =  elrbu(DEF_VECT,idime,inode)   + fact0(DEF_VECT)*elvel(DEF_VECT,idime, inode, 2)
                end do
                do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                   elauu11(DEF_VECT,inode,jnode) =  elauu11(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus) 
                   elauu22(DEF_VECT,inode,jnode) =  elauu22(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus) 
                   elauu33(DEF_VECT,inode,jnode) =  elauu33(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus) 
                end do
             end do
          end do
       end if
    else if( kfl_lumped == 2 ) then 
       !
       ! No time term have been added up to now: add Galerkin term
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * dtinv_loc(DEF_VECT)
          if( ndime == 2 ) then
             do inode = 1, pnode
                elauu11(DEF_VECT,inode,inode) = elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elauu22(DEF_VECT,inode,inode) = elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             end do
          else ! ndime==3
             do inode = 1, pnode
                elauu11(DEF_VECT,inode,inode) = elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elauu22(DEF_VECT,inode,inode) = elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elauu33(DEF_VECT,inode,inode) = elauu33(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             end do
          end if
          do inode = 1, pnode
             do idime = 1,ndime   
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elvel(DEF_VECT,idime,inode,2)
             end do
          end do
       end do
    end if
EOF

    nests[9] =<<EOF
    if( kfl_duatss == 1 ) then
       ellum = 0.0_rp 
       do igaus =1, pgaus
          do inode =1, pnode
             fact0(DEF_VECT)       = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_loc(DEF_VECT) * real(fact_duatss-1,rp)
             ellum(DEF_VECT,inode) = ellum(DEF_VECT,inode) + fact0(DEF_VECT) 
          end do
       end do
    end if
EOF

    nests[10] =<<EOF
    if( ndime == 2 ) then  
       do jnode = 1,pnode
          do inode = 1,pnode
             elauu(DEF_VECT,2*inode-1,2*jnode-1)=elauu11(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,2*inode  ,2*jnode-1)=elauu21(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,2*inode-1,2*jnode  )=elauu12(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,2*inode  ,2*jnode  )=elauu22(DEF_VECT,inode,jnode)
          end do
       end do
    else ! ndime==3
       do jnode = 1,pnode
          do inode = 1,pnode
             elauu(DEF_VECT,3*inode-2,3*jnode-2)=elauu11(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode-2)=elauu21(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode-2)=elauu31(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-2,3*jnode-1)=elauu12(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode-1)=elauu22(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode-1)=elauu32(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-2,3*jnode  )=elauu13(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode-1,3*jnode  )=elauu23(DEF_VECT,inode,jnode)
             elauu(DEF_VECT,3*inode  ,3*jnode  )=elauu33(DEF_VECT,inode,jnode)
          end do
       end do
    end if
EOF

    nests[11] =<<EOF
    if( kfl_press_nsi == 1 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,jnode,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,inode,igaus)+p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                  + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    else
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + p1vec(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    end if
EOF

    nests[12] =<<EOF
    if( kfl_p1ve2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + &
                           gpcar(DEF_VECT,jdime,jnode,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if

EOF

    nests[13] =<<EOF
    if (kfl_regim_nsi == 3 .and. kfl_confi_nsi == 1) then ! Penalization of pressure, never imposse pressure
       do igaus =1, pgaus
          penal(DEF_VECT) = 1.0e-4_rp*gpden(DEF_VECT,igaus)/gpvis(DEF_VECT,igaus)
          do inode =1, pnode
             fact0(DEF_VECT) = penal(DEF_VECT) * p2sca(DEF_VECT,inode,igaus)
             do jnode =inode +1, pnode
                fact1(DEF_VECT) = fact0(DEF_VECT)*gpsha(DEF_VECT,jnode,igaus)
                elapp(DEF_VECT,inode,jnode) = elapp(DEF_VECT,inode,jnode) + fact1(DEF_VECT)
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + fact1(DEF_VECT)
             end do
             elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode, igaus)
             elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact0(DEF_VECT) * gppre(DEF_VECT,igaus)
          end do
       end do
    end if

EOF

    nests[14] =<<EOF
#ifdef matiaslma

    if (ndime==2) then
       do igaus =1, pgaus
          do inode =1,pnode
             p2vec(DEF_VECT,1,inode,igaus)= p2vec(DEF_VECT,1,inode,igaus)*gpden(DEF_VECT,igaus)
             p2vec(DEF_VECT,2,inode,igaus)= p2vec(DEF_VECT,2,inode,igaus)*gpden(DEF_VECT,igaus)
             p2sca(DEF_VECT,inode, igaus) = p2sca(DEF_VECT,inode, igaus)*gpden(DEF_VECT,igaus)
          end do
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
          do inode =1, pnode           
             idof1 = 2*inode-1
             idof2 = idof1+1
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
             do jnode =1, pnode
                elapu(DEF_VECT,jnode, idof1) = elapu(DEF_VECT,jnode, idof1) -fact1(DEF_VECT)             * gpcar(DEF_VECT,1,jnode,igaus)  &
                     &                                  + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                elapu(DEF_VECT,jnode, idof2) = elapu(DEF_VECT,jnode, idof2) -fact1(DEF_VECT)             * gpcar(DEF_VECT,2,jnode,igaus)  &
                     &                                  + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
             end do
          end do
       end do
    else !ndime =3
       do igaus =1, pgaus
          do inode =1,pnode
             p2vec(DEF_VECT,1:ndime,inode,igaus) = p2vec(DEF_VECT,1:ndime,inode,igaus)*gpden(DEF_VECT,igaus)
             p2sca(DEF_VECT,inode, igaus)        = p2sca(DEF_VECT,inode, igaus) *gpden(DEF_VECT,igaus)
          end do
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
          do inode =1, pnode
             idof1 = 3*inode - 2 
             idof2 = idof1+1
             idof3 = idof2+1
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
             do jnode =1, pnode 
                elapu(DEF_VECT,jnode, idof1) = elapu(DEF_VECT,jnode, idof1) - fact1(DEF_VECT)              * gpcar(DEF_VECT,1,jnode,igaus) &
                     &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                elapu(DEF_VECT,jnode, idof2) = elapu(DEF_VECT,jnode, idof2) - fact1(DEF_VECT)              * gpcar(DEF_VECT,2,jnode,igaus) &
                     &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
                elapu(DEF_VECT,jnode, idof3) = elapu(DEF_VECT,jnode, idof3) - fact1(DEF_VECT)              * gpcar(DEF_VECT,3,jnode,igaus) &
                     &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,3,jnode,igaus)
             end do
          end do
       end do
    end if

#else
    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca(DEF_VECT,jnode,igaus)&
                     &                                  + rmomu(DEF_VECT,inode,igaus)       * p2vec(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

#endif
EOF

    nests[15] =<<EOF
    if( kfl_rmom2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do jdime = 1,ndime
                jdofv = (jnode-1)*ndime+jdime
                do inode = 1,pnode
                   do kdime = 1,ndime
                      elapu(DEF_VECT,inode,jdofv) = elapu(DEF_VECT,inode,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec(DEF_VECT,kdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if
EOF

    nests[16] =<<EOF
    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + p2vec(DEF_VECT,idime,jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
       end do
    end do
EOF

    nests[17] =<<EOF
    do igaus = 1,pgaus
       fact1(DEF_VECT) = penal_nsi * gpvol(DEF_VECT,igaus)
       do inode = 1,pnode
          elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
          elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elpre(DEF_VECT,inode,1) 
       end do
    end do
EOF

    nests[18] =<<EOF
    if( ndime == 2 ) then
       do igaus = 1,pgaus
          gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
          gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
          gprhh(DEF_VECT,1,igaus)     = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
          gprhh(DEF_VECT,2,igaus)     = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
          do inode = 1,pnode
             elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
             elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
             elrbp(DEF_VECT,inode)   = elrbp(DEF_VECT,inode)   + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                  &                                            + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus) 
          end do
       end do
    else
       do igaus = 1,pgaus
          gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
          gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
          gprhs_sgs(DEF_VECT,3,igaus) = gprhs_sgs(DEF_VECT,3,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,3,igaus)
          gprhh(DEF_VECT,1, igaus)    = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
          gprhh(DEF_VECT,2, igaus)    = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
          gprhh(DEF_VECT,3, igaus)    = gprhs(DEF_VECT,3,igaus)     + gprhs_sgs(DEF_VECT,3,igaus)
          do inode = 1,pnode
             elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
             elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
             elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,3,igaus)
             elrbp(DEF_VECT,inode)   = elrbp(DEF_VECT,inode)   + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                  &                                            + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus) &
                  &                                            + p2vec(DEF_VECT,3,inode,igaus) * gprhh(DEF_VECT,3,igaus)
          end do
       end do
    end if

EOF

    nests[19] =<<EOF
    if( kfl_p1ve2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                do kdime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        + p1ve2(DEF_VECT,idime,kdime,inode,igaus) * gprhh(DEF_VECT,kdime,igaus)
                end do
             end do
          end do
       end do
    end if
EOF

    nests[20] =<<EOF
    do igaus = 1,pgaus
       fact0(DEF_VECT) = gprhc(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       fact1(DEF_VECT) = gpsp2(DEF_VECT,igaus) * gprh2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       do inode = 1,pnode
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT)        * gpcar(DEF_VECT,idime,inode,igaus) ! ( rhs, tau2' * div(v) )
          end do
          elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)                             ! ( rhs , q)
       end do
    end do
EOF

    nests[21] =<<EOF
    if( kfl_linea_nsi == 2 ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) =  gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
             do inode = 1,pnode
                fact1(DEF_VECT)         = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) )
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) )
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) =  gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
             do inode = 1,pnode
                fact1(DEF_VECT)         = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) & 
                     + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) &
                     + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,1,igaus) )
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) &
                     + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) &
                     + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,2,igaus) )
                elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) &
                     + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,3,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,3,igaus) &
                     + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,3,igaus) )
             end do
          end do
       end if
    end if
EOF

    nests[22] =<<EOF
    if( kfl_stabi_nsi == 1 ) then

       do igaus = 1, pgaus
          do inode = 1,pnode
             fact0(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
             do idime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - gprhs_sgs(DEF_VECT,idime,igaus) * fact0(DEF_VECT)
             end do
          end do
       end do

    end if

#ifdef OPENACC
    end do
#endif
EOF
    
    bubble_nests = []

    bubble_nests[0] =<<EOF
       do igaus = 1,pgaus
          p2sca_bub(DEF_VECT,igaus) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
          do idime = 1,ndime
             p2vec_bub(DEF_VECT,idime,igaus) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
          end do
       end do
EOF

    bubble_nests[1] =<<EOF
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) + p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                            + fact1(DEF_VECT) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + p1vec(DEF_VECT,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       end if
EOF

    bubble_nests[2] =<<EOF
       if( kfl_p1ve2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + &
                           gpcar_bub(DEF_VECT,jdime,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end if
EOF

    bubble_nests[3] =<<EOF
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elaqu(DEF_VECT,1,idof1) = elaqu(DEF_VECT,1,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca_bub(DEF_VECT,igaus) &
                     &                                            + rmomu(DEF_VECT,inode,igaus)       * p2vec_bub(DEF_VECT,idime,igaus)
             end do
          end do
       end do
       if( kfl_rmom2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   do kdime = 1,ndime
                      elaqu(DEF_VECT,1,jdofv) = elaqu(DEF_VECT,1,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec_bub(DEF_VECT,kdime,igaus)
                   end do
                end do
             end do
          end do
       end if
EOF

    bubble_nests[4] =<<EOF
       do igaus = 1,pgaus 
          do inode = 1,pnode
             do idime = 1,ndime
                elapq(DEF_VECT,inode,1) = elapq(DEF_VECT,inode,1) + p2vec(DEF_VECT,idime,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                elaqp(DEF_VECT,1,inode) = elaqp(DEF_VECT,1,inode) + p2vec_bub(DEF_VECT,idime,igaus)   * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
          do idime = 1,ndime
             elaqq(DEF_VECT,1,1)             = elaqq(DEF_VECT,1,1)             + p2vec_bub(DEF_VECT,idime,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
             gprhs_sgs(DEF_VECT,idime,igaus) = gprhs_sgs(DEF_VECT,idime,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,idime,igaus)
             gprhh(DEF_VECT,idime,igaus)     = gprhs(DEF_VECT,idime,igaus)     + gprhs_sgs(DEF_VECT,idime,igaus)
             elrbq(DEF_VECT,1)               = elrbq(DEF_VECT,1)               + p2vec_bub(DEF_VECT,idime,igaus) * gprhh(DEF_VECT,idime,igaus)
          end do
       end do
       do igaus = 1,pgaus
          elrbq(DEF_VECT,1) = elrbq(DEF_VECT,1) + gprhc(DEF_VECT,igaus) * p2sca_bub(DEF_VECT,igaus) ! ( rhs , q_bub)
       end do
EOF

    bubble_nests[5] =<<EOF
       do igaus = 1,pgaus
          elaqq(DEF_VECT,1,1) = elaqq(DEF_VECT,1,1) + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) 
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * elbub(DEF_VECT)
       end do
EOF

    pressure_bubble =<<EOF
    if( maxval(pbubl) == 1 ) then
       !
       ! Initialization
       !
       elauq = 0.0_rp
       elapq = 0.0_rp
       elaqu = 0.0_rp
       elaqp = 0.0_rp
       elaqq = 0.0_rp
       elrbq = 0.0_rp

       #{bubble_nests.join}
    end if
EOF

    set_lang(FORTRAN)
    @kernel = CKernel::new(:includes => "immintrin.h")
    @kernel.procedure = declare_procedure
    get_output.print macro
    get_output.print subroutine_prototype
    get_output.print locals
    get_output.print <<EOF
#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif
EOF
    get_output.print init
    
    @opts[:nests].each{|n|
      get_output.print nests[n-1]
    }

    get_output.print <<EOF
#ifdef OPENACC
    end do
#endif
EOF
    
    get_output.print pressure_bubble

    get_output.print "end subroutine nsi_element_assembly_asgs_oss"
  end 
end
