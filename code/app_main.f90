!==============================================================================
!
! app_main.f90
! Part of the cart_pole_f90 example.
!
! This file contains the 'app_main' function, called by Smarties, where the
! training takes place.
!
! For clarity:
!
! C++     Interface        Fortran
! double  real(c_double)   double precision = real*8 = real(kind=8)
! bool    logical(c_bool)  logical
! int     integer(c_int)   integer
! *       type(c_ptr)      Fortran does not really like explicit pointers
!
!
! Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland. All rights reserved.
! Distributed under the terms of the MIT license.
!
!==============================================================================

include 'smarties.f90'
module app_main_module
  
  implicit none

  include 'mpif.h'

  ! 'app_main' is called from the outside (by the interface function
  ! 'app_main_interface' located in 'main.cpp')
  public app_main

contains

  function app_main(smarties_comm, f_mpicomm) result(result_value) &
    bind(c, name='app_main')
    ! This is the main function called from Smarties when the training starts.
  
    use, intrinsic :: iso_c_binding
    use smarties
    ! Modules
    use iso_fortran_env, Only : error_unit, Int32, Int64
    use global
    use input_output
    use initialization
    use time_integration
    use monitor
    use statistics
    use finalization
    use smarties_stat

    implicit none
  
    type(c_ptr),    intent(in), value :: smarties_comm ! this is the pointer to the Smarties communicator
    integer(c_int), intent(in), value :: f_mpicomm ! this is the MPI_COMM_WORLD handle initialized by Smarties
    !
    integer(c_int) :: result_value

    integer :: rank, numProcs, mpiIerr
    real   (c_double), dimension(NUM_ACTIONS), target :: upper_action_bound, lower_action_bound
     
    logical(c_int)  , dimension(STATE_SIZE),  target :: b_observable
    real   (c_double), dimension(STATE_SIZE),  target :: upper_state_bound, lower_state_bound
    logical(kind=c_bool)   :: bounded
    Integer(Int32)    :: i
  

    write(*,*) 'Fortran side begins'

    ! initialize MPI ranks and check that things are working
    call mpi_comm_rank(f_mpicomm, rank, mpiIerr)
    call mpi_comm_size(f_mpicomm, numProcs, mpiIerr)
    write(*,*) 'rank #', rank, ' of ', numProcs, ' is alive in Fortran'

    nagents = NUM_AGENTS/numProcs/2
    call smarties_setNumAgents(smarties_comm, nagents*2)

    ! inform Smarties about the size of the state and the number of actions it can take
    call smarties_setStateActionDims(smarties_comm, STATE_SIZE, NUM_ACTIONS, 0)
  
    ! OPTIONAL: aciton bounds
    bounded = .true.
    upper_action_bound = (/0.9/)
    lower_action_bound = (/1.1/)
    call smarties_setActionScales(smarties_comm, &
          c_loc(upper_action_bound), c_loc(lower_action_bound), &
          bounded, NUM_ACTIONS, 0)
  
    ! OPTIONAL: hide state variables.
    do i = 1,2
       b_observable(i) = 1 
    end do
    do i = 3,6
       b_observable(i) = 0 
    end do
    call smarties_setStateObservable(smarties_comm, c_loc(b_observable), STATE_SIZE, 0)

    ! train loop
    do while (.true.)
  
      Call initialize(f_mpicomm)

      ! small summary of input parameters
      Call summary
      
      ! temporal loop
      Do istep = 1, nsteps
     
        ! send state/reward and recv action from smarties
        Call send_recv_state_action(smarties_comm)

        ! compute dt based on CFL
        Call compute_dt

        ! time step
        Call compute_time_step_RK3

        ! compute a few statistics
        Call compute_statistics 

        ! output some key values
        Call output_monitor

        ! write snapshot if needed
        Call output_data

      End Do

      ! finalize stuff
      Call finalize
  
    end do ! train loop
  
    result_value = 0
    write(*,*) 'Fortran side ends'

  end function app_main

end module app_main_module
