!-----------------------------------------!
!      Module to finalize FFT and etc     !
!-----------------------------------------!
Module finalization
  
  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  
  ! prevent implicit typing
  Implicit None
  
Contains

  !-------------------------------------------!
  !        finalize FFTW plans and MPI        !
  !-------------------------------------------!
  Subroutine finalize
  
    ! finalize FFTW
    Call fftw_free(cplane_fft)
    Call fftw_mpi_cleanup()
    
    ! deallocate
    deallocate (  k1_global,  k2_global )
    deallocate ( kg1_global, kg2_global )
    deallocate ( x_global ,  y_global ,  z_global   )
    deallocate ( xm_global,  ym_global,  zm_global  )
    deallocate ( xg_global,  yg_global,  zg_global  )

    deallocate (  x ,  y ,  z  )
    deallocate (  xm,  ym,  zm )
    deallocate (  xg,  yg,  zg )

    deallocate ( yg_m  )
    deallocate ( yg_mm )
    
    ! global interior + boundary + ghost points
    !deallocate (U  , V  , W  , P  ) 
    deallocate (U, V, W, P)
    deallocate (Uo , Vo , Wo , Po )
    !deallocate (Uoo, Voo, Woo, Poo)
    deallocate (Uoo, Voo, Woo, Poo)

    ! Auxiliary arrays
    deallocate ( term   ) 
    deallocate ( term_1 ) 
    deallocate ( term_2 )
    deallocate ( term_3 )
    deallocate ( term_4 )

    ! RHS: interior points only
    deallocate ( rhs_uo ) 
    deallocate ( rhs_vo )
    deallocate ( rhs_wo )
    deallocate ( rhs_p  )


    deallocate ( buffer_ui, buffer_vi, buffer_wi, buffer_ci, buffer_ai )
    ! local velocity, ending  z-planes
    deallocate ( buffer_ue, buffer_ve, buffer_we, buffer_ce, buffer_ae )
    ! local pressure z-plane
    deallocate ( buffer_p ) 

    !------------------------Interior communications------------------------!
    deallocate ( buffer_us, buffer_ur )
    deallocate ( buffer_vs, buffer_vr )
    deallocate ( buffer_ws, buffer_wr )
    deallocate ( buffer_ps, buffer_pr ) 
    deallocate ( buffer_as, buffer_ar ) 

    deallocate ( rhs_p_hat ) 
    deallocate ( rhs_aux   ) 

    deallocate ( kxx, kzz ) 
    deallocate ( kmode_map ) 

    deallocate ( imode_map_fft ) 
    deallocate ( kmode_map_fft ) 

    deallocate ( pivot )
    deallocate ( Dyy, M )
    deallocate ( D, DL, DU )

    deallocate ( bc_1, bc_2 )
    deallocate ( bc_1_hat, bc_2_hat )

    deallocate ( weight_y_0, weight_y_1 )

    deallocate (  Umean,  Vmean,  Wmean )
    deallocate ( U2mean, V2mean, W2mean, UVmean)

    deallocate( rk_coef, rk_t )

    deallocate ( Fu1, Fu2, Fu3 )
    deallocate ( Fv1, Fv2, Fv3 )
    deallocate ( Fw1, Fw2, Fw3 )

    If (LES_model>0) Then
       deallocate( Lij, Mij, Sij, Rij, S )
       deallocate( ten_buf   ) 
       deallocate( ten_buf_2 )
       
       ! Filtered velocities
       deallocate( Uf  , Vf  , Wf   )
    End If

    ! Eddy viscosity
    deallocate( nu_t, nu_t_hat )  
    deallocate( alpha_x, alpha_z )

    ! SMARTIES
    deallocate ( actions, rewards, rewards_old, states )

    If ( myid==0 ) Then
       Write(*,*) 'Done!'
    End If
    call MPI_Barrier(comm)
    
  End Subroutine finalize
  
End Module finalization
