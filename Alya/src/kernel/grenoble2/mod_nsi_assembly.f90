module mod_nsi_assembly 

  use def_parameters
  use mod_nsi_element_assembly, only : nsi_element_assembly_asgs_oss
  use mod_nsi_element_assembly, only : nsi_element_assembly_split_oss
  implicit none

  private

  public :: nsi_assembly

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    6/10/2016
  !> @brief   Call to element assembly
  !> @details Call to element assembly
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_assembly(iassembly,mnode,pnode,pgaus,ntens,ncomp_nsi,amatr,rhsid)
    !
    ! Dimensions
    !
    integer(ip), intent(in)  :: iassembly                                         !< Assembly type
    integer(ip), intent(in)  :: mnode                                             !< Max number of nodes
    integer(ip), intent(in)  :: pnode                                             !< Number of nodes
    integer(ip), intent(in)  :: pgaus                                             !< Number of Gauss points
    integer(ip), intent(in)  :: ntens                                             !< Hessian size
    integer(ip), intent(in)  :: ncomp_nsi                                         !< Number of temporal components

    real(rp),    intent(out) :: amatr(*)                                          !< Simulated output
    real(rp),    intent(out) :: rhsid(*)                                          !< Simulated output
    ! 
    ! Element matrices and vectors 
    !
    real(rp)                 :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                 :: elaup(VECTOR_SIZE,pnode*ndime,pnode)              ! Aup
    real(rp)                 :: elapp(VECTOR_SIZE,pnode,pnode)                    ! App
    real(rp)                 :: elapu(VECTOR_SIZE,pnode,pnode*ndime)              ! Apu
    real(rp)                 :: elrbu(VECTOR_SIZE,ndime,pnode)                    ! bu
    real(rp)                 :: elrbp(VECTOR_SIZE,pnode)                          ! bp
    real(rp)                 :: elrhs(VECTOR_SIZE,6*pnode)                        ! Generic RHS  
    real(rp)                 :: ellum(VECTOR_SIZE,pnode)                          ! Lumped mass matrix
    !
    ! Bubble matrices
    !
    real(rp)                 :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp)                 :: elapq(VECTOR_SIZE,pnode,1)
    real(rp)                 :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp)                 :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp)                 :: elaqq(VECTOR_SIZE,1,1)
    real(rp)                 :: elrbq(VECTOR_SIZE,1)
    !
    ! Perturbation and residuals
    !
    real(rp)                 :: rmomu(VECTOR_SIZE,pnode,pgaus)                    ! Residual velocity in momentum
    real(rp)                 :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)        ! Residual velocity in momentum
    real(rp)                 :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)              ! Residual velocity in continuity
    real(rp)                 :: p1vec(VECTOR_SIZE,ndime,ndime,pnode,pgaus)        ! Test funct. velocity in momentum
    real(rp)                 :: p1ve2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)        ! Test funct. velocity in momentum
    real(rp)                 :: p2vec(VECTOR_SIZE,ndime,pnode,pgaus)              ! Test funct. velocity in continuity
    real(rp)                 :: p2sca(VECTOR_SIZE,pnode,pgaus)                    ! Test function pressure in continuity
    !
    ! Gather 
    !
    real(rp)                 :: elvel(VECTOR_SIZE,ndime,pnode,ncomp_nsi)          ! u
    real(rp)                 :: elpre(VECTOR_SIZE,pnode,ncomp_nsi-1)              ! p
    real(rp)                 :: elbub(VECTOR_SIZE)                                ! Element bubble
    !
    ! Gauss point values
    !
    real(rp)                 :: gpsha(VECTOR_SIZE,pnode,pgaus)                    ! N
    real(rp)                 :: gpder(VECTOR_SIZE,ndime,pnode,pgaus)              ! dN/dsi
    real(rp)                 :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)                 :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)              ! d2N/dxidxj
    real(rp)                 :: gplap(VECTOR_SIZE,pnode,pgaus)                    ! Lapl(N)
    real(rp)                 :: gpvol(VECTOR_SIZE,pgaus)                          ! w*|J|, |J|
    real(rp)                 :: gpvis(VECTOR_SIZE,pgaus)                          ! Viscosity 
    real(rp)                 :: gpgvi(VECTOR_SIZE,ndime,pgaus)                    ! Viscosity gradients
    real(rp)                 :: gppor(VECTOR_SIZE,pgaus)                          ! Porosity
    real(rp)                 :: gpden(VECTOR_SIZE,pgaus)                          ! Density
    real(rp)                 :: gpst1(VECTOR_SIZE,pgaus)                          ! tau1
    real(rp)                 :: gpsp1(VECTOR_SIZE,pgaus)                          ! tau1'
    real(rp)                 :: gpsp2(VECTOR_SIZE,pgaus)                          ! tau2'
    real(rp)                 :: gptt1(VECTOR_SIZE,pgaus)                          ! tau1'/tau1
    real(rp)                 :: gptt2(VECTOR_SIZE,pgaus)                          ! tau2'/tau2
    real(rp)                 :: gpadv(VECTOR_SIZE,ndime,pgaus)                    ! u+u'
    real(rp)                 :: gprhs(VECTOR_SIZE,ndime,pgaus)                    ! RHS
    real(rp)                 :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)                ! RHS due to subscales
    real(rp)                 :: gprhc(VECTOR_SIZE,pgaus)                          ! RHS for the continuity equation (Low Mach)
    real(rp)                 :: gprh2(VECTOR_SIZE,pgaus)                          ! RHS for the residual of continuity equation (Low Mach)
    real(rp)                 :: gpsgs(VECTOR_SIZE,ndime,pgaus,2)                  ! u'
    real(rp)                 :: gppre(VECTOR_SIZE,pgaus,ncomp_nsi-1)              ! p
    real(rp)                 :: gpvel(VECTOR_SIZE,ndime,pgaus,ncomp_nsi-1)        ! u
    real(rp)                 :: gpgrp(VECTOR_SIZE,ndime,pgaus)                    ! tau1'*( grad(p) - rho*f )
    real(rp)                 :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                 :: gpvep(VECTOR_SIZE,ndime,pgaus)                    ! -tau1'*R(u) or tau1*rho*(a.grad)u
    !
    ! Enrichement
    !
    real(rp)                 :: gpsha_bub(VECTOR_SIZE,pgaus)                      ! Ne
    real(rp)                 :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)                ! dNe/dxi

    integer(ip)              :: pbubl(VECTOR_SIZE)
    real(rp)                 :: dtinv_loc(VECTOR_SIZE)
    real(rp)                 :: dtsgs_loc(VECTOR_SIZE)
    !
    ! Originally passed by a module
    !
    integer(ip)    :: kfl_lumped
    integer(ip)    :: kfl_duatss
    integer(ip)    :: fact_duatss
    integer(ip)    :: kfl_stabi_nsi
    real(rp)       :: fvins_nsi
    real(rp)       :: fcons_nsi
    real(rp)       :: bemol_nsi
    integer(ip)    :: kfl_regim_nsi
    real(rp)       :: fvela_nsi(3)
    integer(ip)    :: kfl_rmom2_nsi
    integer(ip)    :: kfl_press_nsi
    integer(ip)    :: kfl_p1ve2_nsi
    integer(ip)    :: kfl_linea_nsi
    real(rp)       :: pabdf_nsi(10)
    integer(ip)    :: kfl_confi_nsi
    integer(ip)    :: nbdfp_nsi
    integer(ip)    :: kfl_sgsti_nsi
    integer(ip)    :: kfl_nota1_nsi
    integer(ip)    :: kfl_limit_nsi
    integer(ip)    :: kfl_penal_nsi
    real(rp)       :: penal_nsi
    integer(ip)    :: kfl_bubbl_nsi
    !
    ! Variables from original module
    !
    kfl_lumped    = 0
    kfl_duatss    = 0
    fact_duatss   = 0.0_rp
    kfl_stabi_nsi = 2
    fvins_nsi     = 1.0_rp

    fcons_nsi     = 0.0_rp
    bemol_nsi     = 0.0_rp
    kfl_regim_nsi = 0
    fvela_nsi     = 0.0_rp
    kfl_rmom2_nsi = 0
    kfl_press_nsi = 1
    kfl_p1ve2_nsi = 0
    kfl_linea_nsi = 1
    pabdf_nsi     = 1.0_rp
    kfl_confi_nsi = 0
    nbdfp_nsi     = ncomp_nsi-1
    kfl_sgsti_nsi = 1
    kfl_nota1_nsi = 0
    kfl_limit_nsi = 0
    kfl_penal_nsi = 0
    penal_nsi     = 0.0_rp
    kfl_bubbl_nsi = 0

    gpden     = 1.0_rp
    gpvis     = 1.0_rp
    gppor     = 1.0_rp
    gpgvi     = 1.0_rp
    gpsp1     = 1.0_rp
    gptt1     = 1.0_rp
    gpsp2     = 1.0_rp
    gptt2     = 1.0_rp
    gpvol     = 1.0_rp
    gpsha     = 1.0_rp
    gpcar     = 1.0_rp
    gplap     = 1.0_rp
    gphes     = 1.0_rp
    gpadv     = 1.0_rp
    gpvep     = 1.0_rp
    gprhs     = 1.0_rp
    gprhc     = 1.0_rp
    rmomu     = 1.0_rp
    rcont     = 1.0_rp
    elpre     = 1.0_rp
    elbub     = 1.0_rp
    rmom2     = 1.0_rp
    gpst1     = 1.0_rp
    gpgve     = 1.0_rp
    gprh2     = 1.0_rp
    gprhs_sgs = 1.0_rp
    elvel     = 1.0_rp
    ellum     = 1.0_rp
    dtinv_loc = 1.0_rp
    dtsgs_loc = 1.0_rp

    pbubl     = 0
    gpsha_bub = 1.0_rp
    gpcar_bub = 1.0_rp
    gppre     = 1.0_rp
    gpvel     = 1.0_rp
    gpgrp     = 0.0_rp

    elauu     = 0.0_rp
    elaup     = 0.0_rp
    elapp     = 0.0_rp
    elapu     = 0.0_rp
    elrbu     = 0.0_rp
    elrbp     = 0.0_rp

    elauq     = 0.0_rp 
    elapq     = 0.0_rp
    elaqu     = 0.0_rp
    elaqp     = 0.0_rp
    elaqq     = 0.0_rp
    elrbq     = 0.0_rp
    !
    ! Elemental assembly
    !
    if( iassembly == 1 ) then
       call nsi_element_assembly_split_oss(                  &
            pnode,pgaus,gpden,gpvis,gppor,gpsp1,gpsp2,gpvol, &
            gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,gprhc,gpvel, &
            gpsgs,elvel,elpre,elbub,elauu,elaup,elapp,elapu, &
            elrbu,elrbp,dtinv_loc,dtsgs_loc,pbubl,gpsha_bub, &
            gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq,elrbq,   &
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
    else
       call nsi_element_assembly_asgs_oss(                   &
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
    end if
    !
    ! Simulated global assembly
    !
    call nsi_global_assembly(&
         pnode,elauu,elaup,elapu,elapp,elrbu,elrbp,amatr,rhsid)

  end subroutine nsi_assembly


  subroutine nsi_global_assembly(&
       pnode,elauu,elaup,elapu,elapp,elrbu,elrbp,amatr,rhsid)
    
    integer(ip)                :: pnode
    real(rp),    intent(in)    :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(in)    :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(in)    :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(in)    :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(in)    :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(inout) :: amatr(*)
    real(rp),    intent(inout) :: rhsid(*)
    integer(ip)                :: inode,jnode,idime,jdime
    integer(ip)                :: idofn,jdofn,itotv,ivect
    
    do ivect = 1,VECTOR_SIZE
       itotv = 0
       do inode = 1,pnode
          do idime = 1,ndime
             idofn = (inode-1)*ndime+idime
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdofn = (jnode-1)*ndime+jdime
                   itotv = itotv + 1
                   amatr(itotv) = amatr(itotv) + elauu(ivect,idofn,jdofn)
                end do
             end do
             amatr(itotv) = amatr(itotv) + elaup(ivect,idofn,jnode)
             amatr(itotv) = amatr(itotv) + elapu(ivect,inode,jnode)
             amatr(itotv) = amatr(itotv) + elapp(ivect,inode,jnode)
             rhsid(idofn) = rhsid(idofn) + elrbu(ivect,idime,inode)
          end do
          rhsid(idofn) = rhsid(idofn) + elrbp(ivect,inode)
       end do
    end do

  end subroutine nsi_global_assembly

end module mod_nsi_assembly
