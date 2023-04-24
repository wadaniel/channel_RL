!------------------------------------!
!     Module for LES wall-models     !
!------------------------------------!
Module wallmodel

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use interpolation
  Use subgrid
  Use boundary_conditions

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------!
  !                                              !
  !              Select wall model               !
  !                                              !
  !----------------------------------------------!
  Subroutine compute_wall_model(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx, nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nxg,nyg,nzg), Intent(InOut) :: nu_t_

    Integer(Int32) :: i, j

    term_2(1:nx,1,:) = alpha_x(1:nx,1,:)/&
            ((U_(:,    2,:) - U_(:,  1,:))/(yg(2)-yg(1))) - nu
    term_2(1:nx,2,:) = alpha_x(1:nx,2,:)/&
            ((U_(:,nyg-1,:) - U_(:,nyg,:))/(yg(2)-yg(1))) - nu

    Call interpolate_x(term_2(1:nx,:,:),term(2:nx-1,:,:),2)

    nu_t_(:,    1,:) = term(:,1,:)
    nu_t_(:,    2,:) = term(:,1,:)
    nu_t_(:,nyg-1,:) = term(:,2,:)
    nu_t_(:,nyg  ,:) = term(:,2,:)

    Call apply_periodic_bc_x(nu_t_,2) 
    
  End Subroutine compute_wall_model

  !-------------------------------------------------------------!
  !                                                             !
  !            Compute Neumann boundary conditions              !
  !    for pseudo-pressure when slip-wall model is active       !
  !                                                             !
  ! This has to be called every sub-step                        !
  !                                                             !
  ! Conditions bottom wall:                                     !
  !                                                             !
  !     V1 = V1*  -  (p2-p1)/(yg(2)-yg(1))                      !
  !     V2 = V2*  -  (p3-p2)/(yg(3)-yg(2))                      !
  !                                                             !
  !     => p1 = p_b2*p2 + p_bc3*p3                              !
  !        p_bc2   = 1 + beta*Delta_r                           !
  !        p_bc3   =   - beta*Delta_r                           !
  !        Delta_r = ( yg(2)-yg(1) )/( yg(3)-yg(2) )            !
  !                                                             !
  !     V* -> velocity without pressure                         !
  !                                                             !
  ! Equation for first interior points:                         !
  !                                                             !
  !    (a + c*p_bc3)*p3 + (b + c*p_bc2)*p2 = rhs_p2             !
  !                                                             !
  !                                                             !
  ! Conditions top wall: (n->ny, ng->nyg)                       !
  !                                                             !
  !     V(n)   = V(  n)* -(p(ng)  -p(ng-1))/(yg(ng)  -yg(ng-1)) !
  !     V(n-1) = V(n-1)* -(p(ng-1)-p(ng-2))/(yg(ng-1)-yg(ng-2)) !
  !                                                             !
  !     => p(ng)    = p_bn1*p(ng-2) + p_bcn*p(ng-1)             !
  !        p_bcn    = 1d0 + beta*Delta_r                        !
  !        p_bcn1   =     - beta*Delta_r                        !
  !        Delta_r  = (yg(ng)-yg(ng-1))/(yg(ng-1)-yg(ng-2))     !
  !                                                             !
  ! Equation for last interior points:                          !
  !                                                             !
  !   (b+a*p_bcn)*p(ng-1) + (c+a*p_bcn1)*p(ng-2) = rhs_p2(ng-1) !
  !                                                             !
  !                                                             !
  !-------------------------------------------------------------!
  Subroutine compute_pseudo_pressure_bc_for_robin_bc

    ! local variables
    Real   (Int64) :: a, b, c
    Real   (Int64) :: beta, Delta_r, alphad
    Real   (Int64) :: p_bc2, p_bc3, p_bcn, p_bcn1
    Integer(Int64) :: j    

    ! bottom wall
    j        = 2 
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    Delta_r  = ( yg(2)-yg(1) )/( yg(3)-yg(2) )
    alphad   = 0 
    beta     = alphad/(alphad + y(2)-y(1) )
    p_bc2    = 1d0 + beta*Delta_r
    p_bc3    =     - beta*Delta_r

    Dyy(2,2) = b + c*p_bc2
    Dyy(2,3) = a + c*p_bc3
    
    ! top wall
    j        = nyg-1
    a        = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b        = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c        = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    alphad   = 0 
    Delta_r  = ( yg(nyg) - yg(nyg-1) )/( yg(nyg-1) - yg(nyg-2) )
    beta     = alphad/( alphad - (y(ny)-y(ny-1)) )
    p_bcn    = 1d0 + beta*Delta_r
    p_bcn1   =     - beta*Delta_r
    
    Dyy(nyg-1,nyg-1) = b + a*p_bcn
    Dyy(nyg-1,nyg-2) = c + a*p_bcn1
    
  End Subroutine compute_pseudo_pressure_bc_for_robin_bc
End Module wallmodel

