!----------------------------------------------!
!     Module for interfacing with smarties     !
!----------------------------------------------!
Include 'smarties.f90'

Module smarties_stat

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use boundary_conditions
  Use subgrid
  Use smarties
  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------!
  !                                              !
  !         communicate with smarties            !
  !                                              !
  !----------------------------------------------!
  Subroutine send_recv_state_action(smarties_comm) bind(c, name='send_recv_state_action')
    Use, intrinsic :: iso_c_binding
    type(c_ptr),    intent(in), value :: smarties_comm 
    Integer(Int32) :: i, j, na, iwall, y_ind
    Real(Int64)    :: v_sca, l_sca, v_sca_2, Uloc, dUloc, y_rand, aa

    real   (c_double), dimension(NUM_ACTIONS), target :: action
    real   (c_double), dimension(STATE_SIZE),  target :: state
    real   (c_double) :: reward
    integer :: agent_id 

    If (istep .eq. 1) Then 
       Call compute_eddy_viscosity(U,V,W,nu_t)    
    End If
    ! sum over for each processor
    Do j=1,nyg
      if (myid < nprocs-1) then
        Umean  (j) = Sum( U(2:nx-1, j,2:nzg-1) )
      else
        Umean  (j) = Sum( U(2:nx-1, j,2:nzg-2) )
      end if
    End Do

    ! reduce statatistics between processors      
    If ( myid==0 ) Then
      Call MPI_Reduce(MPI_IN_PLACE,Umean,nyg,MPI_real8,MPI_sum,0,comm,ierr)
    Else
      Call MPI_Reduce(Umean,0,nyg,MPI_real8,MPI_sum,0,comm,ierr)
    End If

    ! compute mean and broadcast 
    Umean  = Umean/Real( ( nx_global-2)*(nzg_global-3), 8)
    Call Mpi_bcast ( Umean, nyg, MPI_real8, 0, comm, ierr )
    Ucl = maxval(Umean)
    v_sca   = Ucl

    Do i = 1,nagents

      Call RANDOM_NUMBER(y_rand)
      If (y_rand < 1d0/2d0) Then
         y_ind = 3
      Else
         y_ind = 4
      End If
      Call RANDOM_NUMBER(y_rand)

      utau = sqrt((nu+0.25d0*(nu_t((i-1)*nx_nagents+3,2,2) + &
                              nu_t((i-1)*nx_nagents+2,2,2) + &
                              nu_t((i-1)*nx_nagents+3,1,2) + &
                              nu_t((i-1)*nx_nagents+2,1,2))) * &
              (U((i-1)*nx_nagents+2,2,2)-U((i-1)*nx_nagents+2,1,2)) / &
              (yg(    2)-yg(  1)))
      v_sca_2 = utau
      l_sca   = y_rand*yg(y_ind)+(1d0-y_rand)*yg(y_ind+1)

      Uloc = 0.2d0 * (     y_rand *U((i-1)*nx_nagents+2,y_ind  ,2)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+2,y_ind+1,2))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+1,y_ind  ,2)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+1,y_ind+1,2))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+3,y_ind  ,2)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+3,y_ind+1,2))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+2,y_ind  ,1)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+2,y_ind+1,1))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+2,y_ind  ,3)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+2,y_ind+1,3)))

      dUloc = 0.2d0 * (U((i-1)*nx_nagents+2,y_ind+1,2) - &
                       U((i-1)*nx_nagents+2,y_ind  ,2) + &
                       U((i-1)*nx_nagents+1,y_ind+1,2) - &
                       U((i-1)*nx_nagents+1,y_ind  ,2) + &
                       U((i-1)*nx_nagents+3,y_ind+1,2) - &
                       U((i-1)*nx_nagents+3,y_ind  ,2) + &
                       U((i-1)*nx_nagents+2,y_ind+1,1) - &
                       U((i-1)*nx_nagents+2,y_ind  ,1) + &
                       U((i-1)*nx_nagents+2,y_ind+1,3) - &
                       U((i-1)*nx_nagents+2,y_ind  ,3)) / &
                       (yg(y_ind+1) - yg(y_ind))

      states(1,i,1) = (dUloc*l_sca)/v_sca_2
      states(1,i,2) = Uloc/v_sca_2 - (dUloc*l_sca)/v_sca_2*log(l_sca*v_sca_2/nu) 

      ! hidden states
      states(1,i,3) = Uloc
      states(1,i,4) = l_sca
      states(1,i,5) = v_sca_2
      states(1,i,6) = dPdx_ref 


      Call RANDOM_NUMBER(y_rand)
      If (y_rand < 1d0/2d0) Then
         y_ind = 3
      Else
         y_ind = 4
      End If
      Call RANDOM_NUMBER(y_rand)

      utau = sqrt((nu+0.25d0*(nu_t((i-1)*nx_nagents+3,nyg  ,2) + &
                              nu_t((i-1)*nx_nagents+2,nyg  ,2) + &
                              nu_t((i-1)*nx_nagents+3,nyg-1,2) + &
                              nu_t((i-1)*nx_nagents+2,nyg-1,2))) * &
              (U((i-1)*nx_nagents+2,nyg-1,2) - U((i-1)*nx_nagents+2,nyg,2))/&
              (yg(2)-yg(1)))
      v_sca_2 = utau
      l_sca   = 2d0-(y_rand*yg(nyg-y_ind+1)+(1d0-y_rand)*yg(nyg-y_ind  ))

      Uloc = 0.2d0 * (     y_rand *U((i-1)*nx_nagents+2,nyg-y_ind+1,2)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+2,nyg-y_ind  ,2))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+1,nyg-y_ind+1,2)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+1,nyg-y_ind  ,2))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+3,nyg-y_ind+1,2)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+3,nyg-y_ind  ,2))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+2,nyg-y_ind+1,1)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+2,nyg-y_ind  ,1))) + &
             0.2d0 * (     y_rand *U((i-1)*nx_nagents+2,nyg-y_ind+1,3)    + &
                     ((1d0-y_rand)*U((i-1)*nx_nagents+2,nyg-y_ind  ,3))) 
      dUloc = 0.2d0 * (U((i-1)*nx_nagents+2,nyg-y_ind  ,2) - &
                       U((i-1)*nx_nagents+2,nyg-y_ind+1,2) + &
                       U((i-1)*nx_nagents+1,nyg-y_ind  ,2) - &
                       U((i-1)*nx_nagents+1,nyg-y_ind+1,2) + &
                       U((i-1)*nx_nagents+3,nyg-y_ind  ,2) - &
                       U((i-1)*nx_nagents+3,nyg-y_ind+1,2) + &
                       U((i-1)*nx_nagents+2,nyg-y_ind  ,1) - &
                       U((i-1)*nx_nagents+2,nyg-y_ind+1,1) + &
                       U((i-1)*nx_nagents+2,nyg-y_ind  ,3) - &
                       U((i-1)*nx_nagents+2,nyg-y_ind+1,3)) / &
                       (yg(nyg-y_ind+1) - yg(nyg-y_ind))

      states(2,i,1) = (dUloc*l_sca)/v_sca_2
      states(2,i,2) = Uloc/v_sca_2 - (dUloc*l_sca)/v_sca_2*log(l_sca*v_sca_2/nu) 

      ! hidden states
      states(2,i,3) = Uloc
      states(2,i,4) = l_sca
      states(2,i,5) = v_sca_2
      states(2,i,6) = dPdx_ref 

      rewards(1,i) = 100d0 - 5000d0* &
                     abs(((nu+0.25*(nu_t((i-1)*nx_nagents+2,1,2)+ &
                                    nu_t((i-1)*nx_nagents+3,1,2)+ &
                                    nu_t((i-1)*nx_nagents+2,2,2)+ &
                                    nu_t((i-1)*nx_nagents+3,2,2))) * &
                         (U((i-1)*nx_nagents+2,2,2)-&
                          U((i-1)*nx_nagents+2,1,2))/&
                          (yg(2)-yg(1)))**0.5d0 &
                         - dPdx_ref**0.5d0) / dPdx_ref**0.5d0 
      rewards(2,i) = 100d0 - 5000d0* &
                     abs(((nu+0.25*(nu_t((i-1)*nx_nagents+2,nyg-1,2)+&
                                    nu_t((i-1)*nx_nagents+3,nyg-1,2)+ &
                                    nu_t((i-1)*nx_nagents+2,nyg  ,2)+ &
                                    nu_t((i-1)*nx_nagents+3,nyg  ,2))) * &
                         (U((i-1)*nx_nagents+2,nyg-1,2)-&
                          U((i-1)*nx_nagents+2,nyg,2))/&
                          (yg(nyg)-yg(nyg-1)))**0.5d0 &
                         - dPdx_ref**0.5d0) / dPdx_ref**0.5d0  
    End Do

    If (istep .lt. 1000) Then
      If (istep .eq. 1) Then
         alpha_x(:,1,:) = dPdx_ref 
         alpha_x(:,2,:) = dPdx_ref 
         If (myid.eq.0) then 
            Call RANDOM_NUMBER(aa)
         End If
         Call Mpi_bcast ( aa, 1, MPI_real8, 0, comm, ierr )
         alpha_x = alpha_x * (1d0+(aa-0.5d0))
      End If
    ElseIf (istep .eq. 1000) Then 
      Do iwall = 1,2
      Do na = 1,nagents
        state = states(iwall,na,:) 
        agent_id = 2*(na-1)+iwall-1
        call smarties_sendInitState(smarties_comm, c_loc(state) ,  STATE_SIZE, agent_id) 
        call smarties_recvAction   (smarties_comm, c_loc(action), NUM_ACTIONS, agent_id) 
        actions(iwall,na,:) = alpha_x((na-1)*nx_nagents+2,iwall,2)*action(1)
      End Do
      End Do
      rewards_old = rewards
    Elseif ((istep .lt. nsteps) .and. (mod(istep,200) .eq. 0)) Then 
      Do iwall = 1,2
      Do na = 1,nagents
        state = states(iwall,na,:) 
        reward = rewards(iwall,na)-rewards_old(iwall,na)
        if (rewards(iwall,na) > 0d0) then
           reward = reward + rewards(iwall,na)
        end if
        agent_id = 2*(na-1)+iwall-1

        call smarties_sendState (smarties_comm, c_loc(state) ,  STATE_SIZE, reward, agent_id)
        call smarties_recvAction(smarties_comm, c_loc(action), NUM_ACTIONS, agent_id) 
        actions(iwall,na,:) = alpha_x((na-1)*nx_nagents+2,iwall,2)*action(1)
      End Do
      End Do
      rewards_old = rewards
    ElseIf (istep .eq. nsteps) Then
      Do iwall = 1,2
      Do na = 1,nagents
        state = states(iwall,na,:) 
        reward = rewards(iwall,na)-rewards_old(iwall,na)
        if (rewards(iwall,na) > 0d0) then
           reward = reward + rewards(iwall,na)
        end if
        agent_id = 2*(na-1)+iwall-1
        call smarties_sendLastState(smarties_comm, c_loc(state), STATE_SIZE, reward, agent_id) 
      End Do
      End Do
    End If

    If ((istep .gt. 999).and.(istep .lt. nsteps)) Then
      ! interpolate along x-direction
      Do j = 1,nx_nagents
        Do i = 1,nagents-1
          alpha_x((i-1)*nx_nagents+j+1,1,:) = actions(1,i,1)*real(nx_nagents-j+1)/real(nx_nagents) + &
                                            actions(1,i+1,1)*real(j-1)/real(nx_nagents)
          alpha_x((i-1)*nx_nagents+j+1,2,:) = actions(2,i,1)*real(nx_nagents-j+1)/real(nx_nagents) + &
                                            actions(2,i+1,1)*real(j-1)/real(nx_nagents)
        End do
        i = nagents
        alpha_x((i-1)*nx_nagents+j+1,1,:) = actions(1,i,1)*real(nx_nagents-j+1)/real(nx_nagents) + &
                                          actions(1,1,1)*real(j-1)/real(nx_nagents)
        alpha_x((i-1)*nx_nagents+j+1,2,:) = actions(2,i,1)*real(nx_nagents-j+1)/real(nx_nagents) + &
                                          actions(2,1,1)*real(j-1)/real(nx_nagents)
      end do
      alpha_z = 0d0

      call apply_periodic_bc_x(alpha_x,1)
      call apply_periodic_bc_z(alpha_x,5)
      call update_ghost_interior_planes(alpha_x,5)

      ! interpolate along z-direction
      do i = 1,nz_nagents
         alpha_x(:,:,i+1) = alpha_x(:,:,  2)*real(nz_nagents-i+1)/real(nz_nagents) + &
                            alpha_x(:,:,nzg)*real(i-1)/real(nz_nagents) 
      end do

    End If

  End Subroutine send_recv_state_action
End Module smarties_stat
