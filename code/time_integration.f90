!---------------------------------------------!
!     Module for temporal integration         !
!---------------------------------------------!
Module time_integration

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use equations
  Use projection
  Use boundary_conditions
  Use mass_flow
  Use subgrid
  Use wallmodel

  ! prevent implicit typing
  Implicit None

Contains
  
  !-----------------------------------------------!
  !          Explicit Runge-Kutta 3 steps         !
  !-----------------------------------------------!
  Subroutine compute_time_step_RK3

    Real(Int64) :: to

    ! save previous state
    to = t
    Uo = U
    Vo = V
    Wo = W

    ! step 1
    rk_step = 1
    Call compute_eddy_viscosity(U,V,W,nu_t)    
    Call compute_wall_model(U,V,W,nu_t) ! uses nu_t from compute_eddy_viscosity
    Call compute_rhs_u(U,V,W,Fu1)
    Call compute_rhs_v(U,V,W,Fv1)
    Call compute_rhs_w(U,V,W,Fw1)

    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fu1
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rk_coef(1,1)*Fv1
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rk_coef(1,1)*Fw1
    t = to + rk_t(rk_step)*dt

    Call apply_boundary_conditions
    Call compute_projection_step
    Call apply_boundary_conditions

    ! step 2
    rk_step = 2
    Call compute_eddy_viscosity(U,V,W,nu_t)    
    Call compute_wall_model(U,V,W,nu_t) ! uses nu_t from compute_eddy_viscosity
    Call compute_rhs_u(U,V,W,Fu2)
    Call compute_rhs_v(U,V,W,Fv2)
    Call compute_rhs_w(U,V,W,Fw2)

    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*( rk_coef(2,1)*Fu1 + rk_coef(2,2)*Fu2 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*( rk_coef(2,1)*Fv1 + rk_coef(2,2)*Fv2 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*( rk_coef(2,1)*Fw1 + rk_coef(2,2)*Fw2 )
    t = to + rk_t(rk_step)*dt

    Call apply_boundary_conditions
    Call compute_projection_step
    Call apply_boundary_conditions

    ! step 3
    rk_step = 3
    Call compute_eddy_viscosity(U,V,W,nu_t)    
    Call compute_wall_model(U,V,W,nu_t) ! uses nu_t from compute_eddy_viscosity
    Call compute_rhs_u(U,V,W,Fu3)
    Call compute_rhs_v(U,V,W,Fv3)
    Call compute_rhs_w(U,V,W,Fw3)

    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fu1 + rk_coef(3,2)*Fu2 + rk_coef(3,3)*Fu3 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fv1 + rk_coef(3,2)*Fv2 + rk_coef(3,3)*Fv3 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + &
         dt*( rk_coef(3,1)*Fw1 + rk_coef(3,2)*Fw2 + rk_coef(3,3)*Fw3 )
    t = to + rk_t(rk_step)*dt

    Call apply_boundary_conditions
    Call compute_projection_step
    Call apply_boundary_conditions

    ! compute mean pressure gradient for constant mass flow in x
    If ( x_mass_cte == 1 ) Then
       Call compute_dPx_for_constant_mass_flow(U,dPdx)
       U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) + dPdx
       dPdx = 0d0 !dPdx/dt ! to be used later by rhs_*
       Call apply_boundary_conditions
    End If

  End Subroutine compute_time_step_RK3

  !-----------------------------------------------!
  !            compute dt based on CFL            !
  !-----------------------------------------------!
  ! NOTE: add rotating and eddy viscosity CFL
  Subroutine compute_dt

    Integer(Int32) :: i, j, k
    Real   (Int64) :: lUmax, lVmax, lWmax, dt_local
    Real   (Int64) :: dt_conv_u, dt_conv_v, dt_conv_w, dt_conv
    Real   (Int64) :: dt_vis_u, dt_vis_v, dt_vis_w, dt_vis
    Real   (Int64) :: dt_max
    
    ! convective time step
    lUmax = 0d0
    lVmax = 0d0
    lWmax = 0d0
    Do i=2,nxg-1
       Do j=2,nyg-1
          Do k=2,nzg-1
             lUmax = Max( lUmax,(xg(i+1)-xg(i))/Abs(U(i,j,k)) )
             lVmax = Max( lVmax,(yg(j+1)-yg(j))/Abs(V(i,j,k)) )
             lWmax = Max( lWmax,(zg(k+1)-zg(k))/Abs(W(i,j,k)) )
          End Do
       End Do
    End Do
    
    dt_conv_u = CFL*lUmax 
    dt_conv_v = CFL*lVmax
    dt_conv_w = CFL*lWmax
    
    dt_conv = Minval( (/dt_conv_u,dt_conv_v,dt_conv_w/) )
    
    ! viscous time step
    dt_vis_u = CFL*dxmin**2d0/nu
    dt_vis_v = CFL*dymin**2d0/nu
    dt_vis_w = CFL*dzmin**2d0/nu
    
    dt_vis = Minval( (/dt_vis_u,dt_vis_v,dt_vis_w/) )
    
    ! time step
    dt_local = Min ( dt_conv,dt_vis )

    ! compute global minimum and communicate results to all processors
    Call MPI_Allreduce(dt_local,dt,1,MPI_real8,MPI_min,comm,ierr)

    ! time step limiter
    dt_max = 1d-3
    dt     = Min( dt, dt_max )
    If ( CFL<0 ) Then
       dt = -CFL
    End If
     
   End Subroutine compute_dt

End Module time_integration
