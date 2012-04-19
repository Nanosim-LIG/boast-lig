! 1 Conv
           print *,'1 Conv'
           call magicfilter_per_1(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_1(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 2 Conv
           print *,'2 Conv'
           call magicfilter_per_2(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_2(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 3 Conv
           print *,'3 Conv'
           call magicfilter_per_3(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_3(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 4 Conv
           print *,'4 Conv'
           call magicfilter_per_4(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_4(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 5 Conv
           print *,'5 Conv'
           call magicfilter_per_5(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_5(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 6 Conv
           print *,'6 Conv'
           call magicfilter_per_6(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_6(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 7 Conv
           print *,'7 Conv'
           call magicfilter_per_7(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_7(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 8 Conv
           print *,'8 Conv'
           call magicfilter_per_8(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_8(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 9 Conv
           print *,'9 Conv'
           call magicfilter_per_9(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_9(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 10 Conv
           print *,'10 Conv'
           call magicfilter_per_10(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_10(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 11 Conv
           print *,'11 Conv'
           call magicfilter_per_11(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_11(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
! 12 Conv
           print *,'12 Conv'
           call magicfilter_per_12(n1,ndat,psi_in,psi_out)

           call start_counters(events,event_number,ierror)
           call nanosec(tsc0);
           call magicfilter_per_12(n1,ndat,psi_in,psi_out)
           call nanosec(tsc1);

           call read_print_stop_counters(counters,event_name,event_number,ierror)
           CPUtime=real(tsc1-tsc0,kind=8)*1d-9
           call compare_2D_results(ndat, n1, psi_out, psi_cuda, maxdiff, 3.d-7)

           call print_time(CPUtime,n1*ndat,32,ntimes)
