  integer, parameter :: event_number = 1
  integer, dimension(event_number) :: events
  data events /&
!                PAPI_L1_DCM,&
!                PAPI_L1_ICM,&
!                 PAPI_TLB_DM,&
                 PAPI_TLB_IM&
!                PAPI_HW_INT,&
!                PAPI_BR_MSP,&
!                PAPI_TOT_IIS,&
!                PAPI_TOT_INS,&
!                PAPI_FP_INS,&
!                PAPI_LD_INS,&
!                PAPI_SR_INS,&
!                PAPI_BR_INS,&
!                PAPI_VEC_INS,&
!                PAPI_TOT_CYC,&
!                PAPI_L1_DCA&
                /
  integer(kind=8), dimension(event_number) :: counters

 character(*), parameter, dimension(event_number) :: event_name = (/&
!              'PAPI_L1_DCM ',&
!               'PAPI_L1_ICM ',&
!               'PAPI_TLB_DM ',&
               'PAPI_TLB_IM '&
!               'PAPI_HW_INT ',&
!               'PAPI_BR_MSP ',&
!                'PAPI_TOT_IIS',&
!                'PAPI_TOT_INS'&
!               'PAPI_FP_INS ',&
!               'PAPI_LD_INS ',&
!               'PAPI_SR_INS ',&
!               'PAPI_BR_INS ',&
!               'PAPI_VEC_INS',&
!               'PAPI_TOT_CYC',&
!               'PAPI_L1_DCA '&
                /)
