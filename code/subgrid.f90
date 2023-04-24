!----------------------------------------!
!    Module for subgrid-scale models     ! 
!----------------------------------------!
Module subgrid

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global 
  Use interpolation
  Use boundary_conditions
  Use mpi

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------------------!
  !                                                      !
  !             Select eddy viscosity model              !
  !                                                      ! 
  !------------------------------------------------------!
  Subroutine compute_eddy_viscosity(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,nzm+2), Intent(Out) :: nu_t_
    If     ( LES_model==1 ) Then
       ! constant coef. Smagorinsky
       Call sgs_Smagorinsky(U_,V_,W_,nu_t_)   
    Elseif ( LES_model==2 ) Then
       ! dynamic Smagorinsky
       Call sgs_dynamic_Smagorinsky(U_,V_,W_,nu_t_) 
    Elseif ( LES_model==3 ) Then
       ! AMD 
       Call sgs_AMD(U_,V_,W_,nu_t_) 
    Elseif ( LES_model==4 ) Then
       ! AMD 
       Call sgs_Vreman(U_,V_,W_,nu_t_) 
    Else
       ! No explicit model
       nu_t_      = 0d0
    End If

  End Subroutine compute_eddy_viscosity

  !--------------------------------------------------------!
  !      Filter field in homogeneous directions xz         !
  !                                                        !
  ! Simpson's rule (4th order):                            !
  !                                                        !
  !   Uf(0,0) = sum( fil(i,j)*U(i,j) )                     !
  !                                                        !
  ! Weights fil:                                           !
  !   1/36  1/9  1/36                                      !
  !   1/9   4/9  1/9                                       !
  !   1/36  1/9  1/36                                      !
  !                                                        !
  ! Note:                                                  !
  !   Must be called as filter_xz(U,Uf(2:n1-1,:,2:n3-1))   !
  !                                                        !
  !                                                        !
  ! Input:  U   (original flow field)                      !
  ! Output: Uf  (filtered flow field)                      !
  !                                                        ! 
  !--------------------------------------------------------!
  Subroutine filter_xz(U_,Uf_)

    Real(Int64), Intent(In)  :: U_  (:,:,:)
    Real(Int64), Intent(Out) :: Uf_ (:,:,:) 

    ! local variables
    Real   (Int64) :: fil(-1:1,-1:1)
    Integer(Int64) :: n(3), n1, n2, n3, i, j, k, ii, kk

    ! local size
    n  = Shape(U_)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ! filter size
    fil_size = 1.5874d0

    ! Set up filter coefficients
    fil(-1,-1) = 1d0/36d0
    fil( 1, 1) = 1d0/36d0
    fil(-1, 1) = 1d0/36d0
    fil( 1,-1) = 1d0/36d0

    fil( 0, 1) = 1d0/9d0
    fil( 1, 0) = 1d0/9d0
    fil( 0,-1) = 1d0/9d0
    fil(-1, 0) = 1d0/9d0

    fil( 0, 0) = 4d0/9d0

    ! Apply filter in homogeneous directions xz
    Uf_ = 0d0 
    ! loop for Uf(ii,:,kk)
    Do ii = 2, n1-1
       Do kk = 2, n3-1
          ! loop for filter
          Do i = -1, 1
             Do k = -1, 1
                ! the -1 shift is because filter_xz is called as filter_xz(U,Uf(2:end-1,:,2:end-1))
                Uf_(ii-1,1:n2,kk-1) = Uf_(ii-1,1:n2,kk-1) + U_(ii+i,1:n2,kk+k)*fil(i,k)
             End do
          End do
       End do
    End do
    
  End Subroutine filter_xz

  !--------------------------------------------------------!
  !             Filter field in xyz directions             !
  !                                                        !
  ! Simpson's rule (4th order):                            !
  !                                                        !
  !   Uf(0,0) = sum( fil(i,j)*U(i,j) )                     !
  !                                                        !
  ! Weights fil in x and z:                                !
  !   1/36  1/9  1/36                                      !
  !   1/9   4/9  1/9                                       !
  !   1/36  1/9  1/36                                      !
  !                                                        !
  ! Weights fil in y:                                      !
  !   1/6   2/3  1/6                                       !
  !                                                        !
  ! Note:                                                  !
  !   Must be called as filter_xzy(U,Uf(2:n1-1,:,2:n3-1))  !
  !                                                        !
  !                                                        !
  ! Input:  U   (original flow field)                      !
  ! Output: Uf  (filtered flow field)                      !
  !                                                        ! 
  !--------------------------------------------------------!
  Subroutine filter_xzy(U_,Uf_)

    Real(Int64), Intent(In)  :: U_  (:,:,:)
    Real(Int64), Intent(Out) :: Uf_ (:,:,:)

    ! local variables
    Real   (Int64) :: fil(-1:1,-1:1)
    Integer(Int64) :: n(3), n1, n2, n3, i, j, k, ii, kk

    ! local size
    n  = Shape(U_)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ! filter size
    fil_size = 2d0

    ! Set up filter coefficients 
    fil(-1,-1) = 1d0/36d0
    fil( 1, 1) = 1d0/36d0
    fil(-1, 1) = 1d0/36d0
    fil( 1,-1) = 1d0/36d0

    fil( 0, 1) = 1d0/9d0
    fil( 1, 0) = 1d0/9d0
    fil( 0,-1) = 1d0/9d0
    fil(-1, 0) = 1d0/9d0

    fil( 0, 0) = 4d0/9d0

    !Apply filter in homogeneous directions
    Uf_ = 0d0
    Do ii = 2, n1-1
       Do i = -1, 1
          Do kk = 2, n3-1
             Do k = -1,1
                Do j = 2, n2-1
                   Uf_(ii-1,j,kk-1) = Uf_(ii-1,j,kk-1) + (1d0/6d0 * U_(ii+i,j-1,kk+k) &
                           + 2d0/3d0 * U_(ii+i,j,kk+k) + 1d0/6d0 * U_(ii+i,j+1,kk+k)) * fil(i,k)
                End Do
                Uf_(ii-1,1,kk-1) = Uf_(ii-1,1,kk-1) + (2d0/3d0 * U_(ii+i,1,kk+k) + &
                        1d0/3d0 * U_(ii+i,2,kk+k)) * fil(i,k)
                Uf_(ii-1,n2,kk-1) = Uf_(ii-1,n2,kk-1) + (1d0/3d0 * U_(ii+i,n2-1,kk+k) + &
                        2d0/3d0 * U_(ii+i,n2,kk+k)) * fil(i,k)
             End do
          End do
       End do
    End do

  End Subroutine filter_xzy

  !--------------------------------------------------------!
  !             Filter field in xyz directions             !
  !                                                        !
  ! Trapezoidal rule (2th order):                          !
  !                                                        !
  !   Uf(0,0) = sum( fil(i,j)*U(i,j) )                     !
  !                                                        !
  ! Note:                                                  !
  !   Must be called as filter_xzy(U,Uf(2:n1-1,:,2:n3-1))  !
  !                                                        !
  !                                                        !
  ! Input:  U   (original flow field)                      !
  ! Output: Uf  (filtered flow field)                      !
  !                                                        ! 
  ! NOTE: this is setting velocities to zero at the wall   !
  !--------------------------------------------------------!
  Subroutine filter_xzy_2nd(U_,Uf_)

    Real(Int64), Intent(In)  :: U_  (:,:,:)
    Real(Int64), Intent(Out) :: Uf_ (:,:,:)

    ! local variables
    Real   (Int64) :: Umean_center(8)
    Integer(Int64) :: n(3), n1, n2, n3, i, j, k, ii, kk, jj

    ! local size
    n  = Shape(U_)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ! filter size
    fil_size = 2d0

    ! apply filter 
    Uf_ = 0d0
    Do ii = 2, n1-1
       Do jj = 2, n2-1
          Do kk = 2, n3-1
             Umean_center = 0d0 
             ! (1,1,1)
             Do i = 0,1
                Do j = 0,1
                   Do k = 0,1
                      Umean_center(1) = Umean_center(1) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (-1,1,1)
             Do i = -1,0
                Do j = 0,1
                   Do k = 0,1
                      Umean_center(2) = Umean_center(2) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (1,-1,1)
             Do i = 0,1
                Do j = -1,0
                   Do k = 0,1
                      Umean_center(3) = Umean_center(3) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (1,1,-1)
             Do i = 0,1
                Do j = 0,1
                   Do k = -1,0
                      Umean_center(4) = Umean_center(4) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (1,-1,-1)
             Do i = 0,1
                Do j = -1,0
                   Do k = -1,0
                      Umean_center(5) = Umean_center(5) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (-1,1,-1)
             Do i = -1,0
                Do j = 0,1
                   Do k = -1,0
                      Umean_center(6) = Umean_center(6) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (-1,-1,1)
             Do i = -1,0
                Do j = -1,0
                   Do k = 0,1
                      Umean_center(7) = Umean_center(7) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             ! (-1,-1,-1)
             Do i = -1,0
                Do j = -1,0
                   Do k = -1,0
                      Umean_center(8) = Umean_center(8) + U_(ii+i,jj+j,kk+k) 
                   End Do
                End Do
             End Do
             !
             Umean_center      = Umean_center/8d0
             Uf_(ii-1,jj,kk-1) = Sum(Umean_center)/8d0
             !
          End do
       End do
    End do

  End Subroutine filter_xzy_2nd

  !------------------------------------------------------!
  !      Filter tensor in homogeneous directions xz      !
  !                                                      !
  ! Input:  T  (4th component is tensor position)        !
  ! Output: Tf                                           !
  !                                                      !
  !------------------------------------------------------!
  Subroutine filter_tensor_xz(T_,Tf_)

    Real(Int64), Intent(In)  :: T_  (:,:,:,:)
    Real(Int64), Intent(Out) :: Tf_ (:,:,:,:) 

    ! local variables
    Integer(Int64) :: n(4), n4, i

    ! local size
    n  = Shape(T_)
    n4 = n(4)

    ! filtering
    Do i = 1, n4
      Call filter_xz(T_(:,:,:,i),Tf_(:,:,:,i))
    End do 

  End Subroutine filter_tensor_xz

  !--------------------------------------------------------!
  !      Filter tensor in homogeneous directions xz and y  !
  !                                                        !
  ! Input:  T  (4th component is tensor position)          !
  ! Output: Tf                                             !
  !                                                        !
  !--------------------------------------------------------!
  Subroutine filter_tensor_xzy(T_,Tf_)

    Real(Int64), Intent(In)  :: T_  (:,:,:,:)
    Real(Int64), Intent(Out) :: Tf_ (:,:,:,:)

    ! local variables
    Integer(Int64) :: n(4), n4, i

    ! local size
    n  = Shape(T_)
    n4 = n(4)

    ! filtering
    Do i = 1, n4
      Call filter_xzy(T_(:,:,:,i),Tf_(:,:,:,i))
    End do

  End Subroutine filter_tensor_xzy

  !------------------------------------------------------------!
  !                                                            !
  !        Compute dynamic Smagorinsky eddy-viscosity          !
  !                     Lilly's approach                       !
  !                                                            !
  ! nu_t = -0.5 * <(Lij * Mij> / <Mij * Mij>*S  (clipping)     !
  ! Lij  = hat(ui*uj) - hat(ui)*hat(uj)                        !
  ! Mij  = fil_size^2 * |hat(S)| * hat(Sij) - hat(|S| * Sij)   !
  ! |S|  = sqrt(2 Sij Sij)                                     !
  ! fil_size = (2*dx*2*dz*2*dy)^(1/3)/(dx*dz*dy)^(1/3)         !
  !                                                            !
  ! Tensors are organized in arrays as                         !
  !   ( 1 4 5 )                                                !
  !   ( 4 2 6 )                                                !
  !   ( 5 6 3 )                                                !
  ! where the number is the 4th component of the array         !
  !                                                            !
  ! Input:  U_,V_,W_ (velocities)                              !
  ! Output: nu_t_                                              !
  !                                                            !
  !------------------------------------------------------------!
  Subroutine sgs_dynamic_Smagorinsky(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,nzm+2), Intent(Out) :: nu_t_

    ! local variables
    Real   (Int64), Dimension(1,nyg,1) :: LijMij_local, LijMij, MijMij
    Integer(Int32) :: i, j, k

    !------------------------------------------------------------------!
    ! Part 1: Compute Leonard term Lij = hat(ui*uj) - hat(ui)*hat(uj)
    ! interpolate velocity to cell centers (faces to centers)
    Call interpolate_x(U_,term  (2:nxg-1,:,:),1) 
    Call interpolate_y(V_,term_1(:,2:nyg-1,:),1) 
    Call interpolate_z(W_,term_2(:,:,2:nzg),1) 

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(term,  2)
    Call apply_periodic_bc_z(term_2,4)
    Call update_ghost_interior_planes(term_2,4)

    ! Lij (Sij as placeholder) = ui*uj at cell centers (why to nzg?)
    Sij(:,2:nyg-1,:,1) = term  (2:nxg,2:nyg-1,2:nzg) * term  (2:nxg,2:nyg-1,2:nzg)  ! u^2
    Sij(:,2:nyg-1,:,2) = term_1(2:nxg,2:nyg-1,2:nzg) * term_1(2:nxg,2:nyg-1,2:nzg)  ! v^2
    Sij(:,2:nyg-1,:,3) = term_2(2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! w^2
    Sij(:,2:nyg-1,:,4) = term  (2:nxg,2:nyg-1,2:nzg) * term_1(2:nxg,2:nyg-1,2:nzg)  ! uv
    Sij(:,2:nyg-1,:,5) = term  (2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! uw
    Sij(:,2:nyg-1,:,6) = term_1(2:nxg,2:nyg-1,2:nzg) * term_2(2:nxg,2:nyg-1,2:nzg)  ! vw

    ten_buf(2:nxg,2:nyg,2:nzg,:) = Sij;
    Do i = 1,6
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do
    ! Lij = hat(ui*uj)
    Call filter_tensor_xzy(ten_buf(1:nxg,2:nyg-1,1:nzg,:),Lij(2:nxg-1,2:nyg-1,2:nzg-1,:)) ! this stores filtered values in Lij(3:nxg-1,2:nyg-1,3:nzg-1)

    ! filter interpolated velocity
    Call filter_xzy(term  (1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,1)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,1)
    Call filter_xzy(term_1(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)
    Call filter_xzy(term_2(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,3)) ! this stores filtered values in Sij(2:nxg-1,2:nyg-1,2:nzg-1,2)

    ! Lij = hat(ui*uj) - hat(ui)*hat(uj)
    ! Only valid on (2:nxg-1,2:nyg-1,2:nzg-1)
    Lij(:,:,:,1) = Lij(:,:,:,1) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,1) ! hat(u^2) - u^2
    Lij(:,:,:,2) = Lij(:,:,:,2) - Sij(2:nxg,2:nyg-1,2:nzg,2) * Sij(2:nxg,2:nyg-1,2:nzg,2) ! hat(v^2) - v^2
    Lij(:,:,:,3) = Lij(:,:,:,3) - Sij(2:nxg,2:nyg-1,2:nzg,3) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat(w^2) - w^2
    Lij(:,:,:,4) = Lij(:,:,:,4) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,2) ! hat( uv) - uv
    Lij(:,:,:,5) = Lij(:,:,:,5) - Sij(2:nxg,2:nyg-1,2:nzg,1) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat( uw) - uw
    Lij(:,:,:,6) = Lij(:,:,:,6) - Sij(2:nxg,2:nyg-1,2:nzg,2) * Sij(2:nxg,2:nyg-1,2:nzg,3) ! hat(u^2) - u^2

    !------------------------------------------------------------------!
    ! Part 2: Compute Mij = fil_size^2*|hat(S)|*hat(Sij) - hat(|S|*Sij) 

    ! filter velocity Same indices as U_
    Call filter_xzy(U_(1:nx,1:nyg,1:nzg),Uf(2:nx-1,1:nyg,2:nzg-1)) ! this stores filtered values in Uf(2:nx-1,2:nyg-1,2:nzg-1)
    Call filter_xzy(V_(1:nxg,1:ny,1:nzg),Vf(2:nxg-1,1:ny,2:nzg-1)) ! this stores filtered values in Vf(2:nxg-1,1:ny,2:nzg-1)
    Call filter_xzy(W_(1:nxg,1:nyg,1:nz),Wf(2:nxg-1,1:nyg,2:nz-1)) ! this stores filtered values in Wf(2:nxg-1,2:nyg-1,2:nz-1)

    ! fill in missing values (periodicity)
    Call apply_periodic_bc_x(Uf,1)
    Call apply_periodic_bc_z(Uf,1)
    Call update_ghost_interior_planes(Uf,1)

    Call apply_periodic_bc_x(Vf,2)
    Call apply_periodic_bc_z(Vf,2)
    Call update_ghost_interior_planes(Vf,2)

    Call apply_periodic_bc_x(Wf,2)
    Call apply_periodic_bc_z(Wf,3)
    Call update_ghost_interior_planes(Wf,3)

    ! Compute hat(Sij) with filtered velocities
    Call compute_Sij(Uf,Vf,Wf,Sij,S)

    ! Mij = fil_size^2 * |hat(S)| * hat(Sij)
    Do i = 1,6
      Mij(:,:,:,i) = fil_size**2d0 * S * Sij(2:nxg-1,2:nyg-1,2:nzg-1,i)
    End Do

    ! Compute Sij with unfiltered velocities
    Call compute_Sij(U_,V_,W_,Sij,S)

    ! Compute hat(|S|Sij) and subtract from Mij
    Do i = 1,6
      term(2:nxg-1,2:nyg-1,2:nzg-1) = S * Sij(2:nxg-1,2:nyg-1,2:nzg-1,i)
      call apply_periodic_bc_x(term,2)
      call apply_periodic_bc_z(term,4)
      call update_ghost_interior_planes(term,4)
      Call filter_xzy(term(1:nxg,2:nyg-1,1:nzg),Sij(2:nxg-1,2:nyg-1,2:nzg-1,i))
    End Do

    ! Mij = fil_size^2*|hat(S)|*hat(Sij) - hat(S*Sij)
    Mij = Mij - Sij(2:nxg-1,2:nyg-1,2:nzg-1,:)

    !------------------------------------------------------------------!
    ! Part 3: Compute eddy viscosity nu_t = Mij*Lij / Mij*Mij

    ! compute Lij*Mij -> stored in nu_t_
    nu_t_ = 0d0
    nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = ( Lij(2:nxg-1,:,2:nzg-1,1)*Mij(:,:,:,1)  +      &
                                       Lij(2:nxg-1,:,2:nzg-1,2)*Mij(:,:,:,2)  +      &
                                       Lij(2:nxg-1,:,2:nzg-1,3)*Mij(:,:,:,3)  +      &
                                 2d0*( Lij(2:nxg-1,:,2:nzg-1,4)*Mij(:,:,:,4)  +      &
                                       Lij(2:nxg-1,:,2:nzg-1,5)*Mij(:,:,:,5)  +      &
                                       Lij(2:nxg-1,:,2:nzg-1,6)*Mij(:,:,:,6) ) )

    ! compute <Lij*Mij> -> LijMij_local
    LijMij_local = 0d0
    Do j = 2,nyg-1
       LijMij_local(1,j,1) = Sum( nu_t_(2:nxg-1,j,2:nzg-1) ) 
    End Do
    Call MPI_Allreduce(LijMij_local,LijMij,nyg,MPI_real8,MPI_sum,comm,ierr)

    ! compute Mij*Mij
    nu_t_(2:nxg-1,2:nyg-1,2:nzg-1) = ( Mij(:,:,:,1)*Mij(:,:,:,1)  +      &
                                       Mij(:,:,:,2)*Mij(:,:,:,2)  +      & 
                                       Mij(:,:,:,3)*Mij(:,:,:,3)  +      & 
                                 2d0*( Mij(:,:,:,4)*Mij(:,:,:,4)  +      & 
                                       Mij(:,:,:,5)*Mij(:,:,:,5)  +      & 
                                       Mij(:,:,:,6)*Mij(:,:,:,6) ) )

    ! compute <Mij*Mij>  -> LijMij_local
    LijMij_local = 0d0
    Do j = 2,nyg-1
       LijMij_local(1,j,1) = Sum( nu_t_(2:nxg-1,j,2:nzg-1) )
    End Do
    Call MPI_Allreduce(LijMij_local,MijMij,nyg,MPI_real8,MPI_sum,comm,ierr)

    ! compute <Lij*Mij>/<Mij*Mij> -> LijMij_local 
    Do j = 2,nyg-1
       LijMij_local(1,j,1) = LijMij(1,j,1) / MijMij(1,j,1) 
    End Do

    ! compute -0.5*<Lij*Mij>/<Mij*Mij>*|S|
    Do j = 2,nyg-1
       nu_t_(2:nxg-1,j,2:nzg-1) = -0.5d0 *LijMij_local(1,j,1) * S(2:nxg-1,j,2:nzg-1) 
    End Do

    ! boundary conditions 
    call apply_periodic_bc_x(nu_t_,2)
    call apply_periodic_bc_z(nu_t_,4)
    Call update_ghost_interior_planes(nu_t_,4)

    ! clipping negative values
    nu_t_ = Max( nu_t_, 0d0 )

    ! boundary conditions ghost cell (must be done after clipping)
    If ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(nu_t_    ,2)
    Else
       Call apply_Neumann_bc_y  (nu_t_    ,2)
    End If

  End Subroutine sgs_dynamic_Smagorinsky

  !-----------------------------------------------------------!
  !                                                           !
  !  Compute constant coefficient Smagorinsky eddy-viscosity  !
  !       with van Driest damping function at the wall        !
  !                                                           !
  ! nu_t  = (Cs*Delta*f)^2 |S|                                !
  ! Cs    = 0.11 (usual range 0.1-0.2 )                       !
  ! Delta = (Cell V_lume)^1/3                                 !
  ! |S|   = sqrt(2 Sij Sij)                                   !
  !                                                           !
  ! van Driest damping function:                              !
  ! f = 1 - exp(-y+/25)                                       !
  ! + is wall-units computed from pressure gradient           !
  !                                                           !
  ! Tensors are organized in arrays as                        !
  !   ( 1 4 5 )                                               !
  !   ( 4 2 6 )                                               !
  !   ( 5 6 3 )                                               !
  ! where the number denotes the 4th component of the array   !
  !                                                           !
  ! Input:  U_,V_,W_ (velocities)                             !
  ! Output: nu_t,                                             !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine sgs_Smagorinsky(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,nzm+2), Intent(Out) :: nu_t_

    ! local variables
    Real   (Int64) :: dum, maxerr, Cs, Delta, f, utau_t
    Integer(Int32) :: i, j, k

    Stop 'ERROR: NEED TO PARALLELIZED Smagorinsky constant coefficient'

    ! Smagorinsky constant coefficient
    Cs = 0.11d0
    
    ! utau from pressure gradient
    utau_t = dPdx**0.5d0

    !------------------------------------------------------------------!
    ! Part 1: Compute rate-of-strain tensor Sij

    ! Compute rate-of-strain and S at cell centers
    Call compute_Sij(U_,V_,W_,Sij,S)

    !------------------------------------------------------------------!
    ! Part 2: Compute eddy viscosity nu_t = (Cs*Delta*f)^2|S| at cell centers
    Do j = 2, nyg-1
       If ( j<nyg/2 ) Then
          f = 1d0 - dexp( -( yg(j)-y(1)  )/25d0/(nu/utau_t) )
       Else
          f = 1d0 - dexp( -( y(ny)-yg(j) )/25d0/(nu/utau_t) )
       End if
       If ( Dirichlet_nu_t == 0 ) f = 1d0
       Delta                   = ( dx*dz*( y(j) - y(j-1) ) )**(1d0/3d0)
       nu_t_(2:nxg-1,j,2:nzg-1) = ( (Cs*Delta*f)**2d0 ) * S(2:nxg-1,j,2:nzg-1)
    End Do
 
    ! clipping negative values
    nu_t_ = Max( nu_t_, 0d0 )

    ! boundary conditions ghost cell (must be done after clipping)
    If ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(nu_t_,2)
    Else
       Call apply_Neumann_bc_y  (nu_t_,2)
    End If

  End Subroutine sgs_Smagorinsky

  !------------------------------------------------------------!
  !                                                            !
  !  Compute annisotropic minimum dissipation model            ! 
  !                                                            !
  ! nu_t_  = C * -(d_hat_i u_j)^2*S_ij / (d_k u_l)*(d_k u_l)   !
  ! C      = 0.3 , the Poincare constant for 2nd order methods !
  ! d_hat_i= the grid weighted derivative                      !
  !                                                            !
  ! Tensors are organized in arrays as                         !
  !   ( 1 4 5 )                                                !
  !   ( 4 2 6 )                                                !
  !   ( 5 6 3 )                                                !
  ! where the number denotes the 4th component of the array    !
  !                                                            !
  ! Input:  U_,V_,W_ (velocities)                              !
  ! Output: nu_t_,                                             !
  !                                                            !
  !------------------------------------------------------------!
  Subroutine sgs_AMD(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,nzg), Intent(Out) :: nu_t_

    Real(Int64)    :: C
    Integer(Int32) :: i, j, k

    ! Poincare constant
    C = 0.3d0   
 
    !------------------------------------------------------------------!
    ! Part 1: Compute rate-of-strain tensor Sij

    ! Compute rate-of-strain and S at cell centers
    Call compute_Sij(U_,V_,W_,Sij,S)

    !------------------------------------------------------------------!
    ! Part 2: Compute the numerator and denominator
    ! nu_t_  = (d_hat_i U) * (d_hat_i U) * S11
    ! term_4 = (d_i     U) * (d_i     U)
    Do i = 2,nx
      term_2(i,:,:) = (U_(i,:,:)-U_(i-1,:,:))
      term_4(i,:,:) = (U_(i,:,:)-U_(i-1,:,:))**2d0 / (x(i)-x(i-1))**2d0
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,1) 

    Do j = 2,nyg
      term  (1:nx,j,:) = (U_(:,j,:)-U_(:,j-1,:))
      term_1(1:nx,j,:) = term(1:nx,j,:) / (yg(j) - yg(j-1))
    End Do
    Call interpolate_x(term_1(1:nx,2:nyg,:),term_2(2:nx,2:nyg  ,:),1)
    Call interpolate_y(term_2(2:nx,2:nyg,:),term_3(2:nx,2:nyg-1,:),1) 
    term_4 = term_4 + term_3**2d0
    Call interpolate_x(term  (1:nx,2:nyg,:),term_1(2:nx,2:nyg  ,:),1)
    Call interpolate_y(term_1(2:nx,2:nyg,:),term_2(2:nx,2:nyg-1,:),1) 
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,1) 

    Do k = 2,nzg
      term  (1:nx,:,k) = (U_(:,:,k)-U_(:,:,k-1))
      term_1(1:nx,:,k) = term(1:nx,:,k) / (zg(k)-zg(k-1))
    End Do
    Call interpolate_x(term_1(1:nx,:,2:nzg),term_2(2:nx,:,2:nzg  ),1)
    Call interpolate_z(term_2(2:nx,:,2:nzg),term_3(2:nx,:,2:nzg-1),1) 
    term_4 = term_4 + term_3**2d0

    Call interpolate_x(term  (1:nx,:,2:nzg),term_1(2:nx,:,2:nzg  ),1)
    Call interpolate_z(term_1(2:nx,:,2:nzg),term_2(2:nx,:,2:nzg-1),1) 
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,1) 
    
    ! nu_t_  <- nu_t_  + (d_hat_i V) * (d_hat_i V) * S22
    ! term_4 <- term_4 + (d_i     V) * (d_i     V) 
    Do j = 2,ny
      term_2(:,j,:) = (V_(:,j,:)-V_(:,j-1,:))
      term_4(:,j,:) = term_4(:,j,:) + term_2(:,j,:)**2d0 / (y(j)-y(j-1))**2d0 
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,2) 

    Do i = 2,nxg
      term  (i,1:ny,:) = (V_(i,:,:)-V_(i-1,:,:))
      term_1(i,1:ny,:) = term(i,1:ny,:) / (xg(i) - xg(i-1))
    End Do

    Call interpolate_y(term_1(2:nxg,1:ny,:),term_2(2:nxg  ,2:ny,:),1)
    Call interpolate_x(term_2(2:nxg,2:ny,:),term_3(2:nxg-1,2:ny,:),1) 
    term_4 = term_4 + term_3**2d0

    Call interpolate_y(term  (2:nxg,1:ny,:),term_1(2:nxg  ,2:ny,:),1)
    Call interpolate_x(term_1(2:nxg,2:ny,:),term_2(2:nxg-1,2:ny,:),1) 
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,2) 

    Do k = 2,nzg
      term  (:,1:ny,k) = (V_(:,:,k)-V_(:,:,k-1))
      term_1(:,1:ny,k) = term(:,1:ny,k) / (zg(k) - zg(k-1))
    End Do
    Call interpolate_y(term_1(:,1:ny,2:nzg),term_2(:,2:ny,2:nzg  ),1)
    Call interpolate_z(term_2(:,2:ny,2:nzg),term_3(:,2:ny,2:nzg-1),1) 
    term_4 = term_4 + term_3**2d0

    Call interpolate_y(term  (:,1:ny,2:nzg),term_1(:,2:ny,2:nzg  ),1)
    Call interpolate_z(term_1(:,2:ny,2:nzg),term_2(:,2:ny,2:nzg-1),1) 
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,2) 
    
    ! nu_t_  <- nu_t_  + (d_hat_i W) * (d_hat_i W) * S33
    ! term_4 <- term_4 + (d_i     W) * (d_i     W) * S33
    Do k = 2,nz
      term_2(:,:,k) = (W_(:,:,k)-W_(:,:,k-1))**2d0
      term_4(:,:,k) = term_4(:,:,k) + term_2(:,:,k) / (z(k) - z(k-1))**2d0 
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,3) 

    Do i = 2,nxg
      term  (i,:,1:nz) = (W_(i,:,:)-W_(i-1,:,:))
      term_1(i,:,1:nz) = term(i,:,1:nz) / (xg(i) - xg(i-1))
    End Do
    Call interpolate_z(term_1(2:nxg,:,1:nz),term_2(2:nxg  ,:,2:nz),1)
    Call interpolate_x(term_2(2:nxg,:,2:nz),term_3(2:nxg-1,:,2:nz),1) 
    term_4 = term_4 + term_3**2d0

    Call interpolate_z(term  (2:nxg,:,1:nz),term_1(2:nxg  ,:,2:nz),1)
    Call interpolate_x(term_1(2:nxg,:,2:nz),term_2(2:nxg-1,:,2:nz),1)
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,3) 

    Do j = 2,nyg
      term  (:,j,1:nz) = (W_(:,j,:)-W_(:,j-1,:))
      term_1(:,j,1:nz) = term(:,j,1:nz) / (yg(j) - yg(j-1)) 
    End Do
    Call interpolate_z(term_1(:,2:nyg,1:nz),term_2(:,2:nyg  ,2:nz),1)
    Call interpolate_y(term_2(:,2:nyg,2:nz),term_3(:,2:nyg-1,2:nz),1) 
    term_4 = term_4 + term_3**2d0

    Call interpolate_z(term  (:,2:nyg,1:nz),term_1(:,2:nyg  ,2:nz),1)
    Call interpolate_y(term_1(:,2:nyg,2:nz),term_2(:,2:nyg-1,2:nz),1)
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + term_2(2:nxg,2:nyg,2:nzg)**2d0 * Sij(:,:,:,3) 

    ! nu_t_  <- nu_t_  + 2 * (d_hat_i U) * (d_hat_i V) * S12
    Do i = 2,nxg
      term(i,1:ny,:) = (V_(i,:,:)-V_(i-1,:,:))
    End Do
    Call interpolate_y(term  (2:nxg,1:ny,:),term_1(2:nxg  ,2:ny,:),1)
    Call interpolate_x(term_1(2:nxg,2:ny,:),term_2(2:nxg-1,2:ny,:),1) 
    Do i = 2,nx
      term_1(i,2:nyg-1,:) = term_2(i,2:nyg-1,:) * (U_(i,2:nyg-1,:)-U_(i-1,2:nyg-1,:))
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,4) 

    Do j = 2,nyg
      term(1:nx,j,:) = (U_(:,j,:)-U_(:,j-1,:))
    End Do
    Call interpolate_x(term  (1:nx,2:nyg,:),term_1(2:nx,2:nyg  ,:),1)
    Call interpolate_y(term_1(2:nx,2:nyg,:),term_2(2:nx,2:nyg-1,:),1) 
    Do j = 2,ny
      term_1(2:nxg-1,j,:) = term_2(2:nxg-1,j,:) * (V_(2:nxg-1,j,:)-V_(2:nxg-1,j-1,:))
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,4) 

    Do k = 2,nzg
      term(1:nx,:,k) = (U_(:,:,k)-U_(:,:,k-1))
    End Do
    Call interpolate_x(term(1:nx,:,2:nzg),term_1(2:nx,:,2:nzg),1)
    Do k = 2,nzg
      term(:,1:ny,k) = (V_(:,:,k)-V_(:,:,k-1))
    End Do
    Call interpolate_y(term(:,1:ny,2:nzg),term_2(:,2:ny,2:nzg),1)
    term_3 = term_1*term_2
    Call interpolate_z(term_3(2:nxg-1,2:nyg-1,2:nzg),term_1(2:nxg-1,2:nyg-1,2:nzg-1),1)
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,4) 

    ! nu_t_  <- nu_t_  + 2 * (d_hat_i U) * (d_hat_i W) * S13
    Do i = 2,nxg
      term(i,:,1:nz) = (W_(i,:,:)-W_(i-1,:,:))
    End Do
    Call interpolate_z(term  (2:nxg,:,1:nz),term_1(2:nxg  ,:,2:nz),1)
    Call interpolate_x(term_1(2:nxg,:,2:nz),term_2(2:nxg-1,:,2:nz),1) 
    Do i = 2,nx
      term_1(i,2:nyg-1,:) = term_2(i,2:nyg-1,:) * (U_(i,2:nyg-1,:)-U_(i-1,2:nyg-1,:))
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,5) 

    Do k = 2,nzg
      term(1:nx,:,k) = (U_(:,:,k)-U_(:,:,k-1))
    End Do
    Call interpolate_x(term  (1:nx,:,2:nzg),term_1(2:nx,:,2:nzg  ),1)
    Call interpolate_z(term_1(2:nx,:,2:nzg),term_2(2:nx,:,2:nzg-1),1) 
    Do k = 2,nz
      term_1(2:nxg-1,:,k) = term_2(2:nxg-1,:,k) * (W_(2:nxg-1,:,k)-W_(2:nxg-1,:,k-1))
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,5) 

    Do j = 2,nyg
      term(1:nx,j,:) = (U_(:,j,:)-U_(:,j-1,:))
    End Do
    Call interpolate_x(term(1:nx,2:nyg,:),term_1(2:nx,2:nyg,:),1)
    Do j = 2,nyg
      term(:,j,1:nz) = (W_(:,j,:)-W_(:,j-1,:))
    End Do
    Call interpolate_z(term(:,2:nyg,1:nz),term_2(:,2:nyg,2:nz),1)
    term_3 = term_1*term_2
    Call interpolate_y(term_3(2:nxg-1,2:nyg,2:nzg-1),term_1(2:nxg-1,2:nyg-1,2:nzg-1),1)
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,5) 

    ! nu_t_  <- nu_t_  + 2 * (d_hat_i V) * (d_hat_i W) * S23
    Do j = 2,nyg
      term(:,j,1:nz) = (W_(:,j,:)-W_(:,j-1,:))
    End Do
    Call interpolate_z(term  (:,2:nyg,1:nz),term_1(:,2:nyg  ,2:nz),1)
    Call interpolate_y(term_1(:,2:nyg,2:nz),term_2(:,2:nyg-1,2:nz),1) 
    Do j = 2,ny
      term_1(2:nxg-1,j,:) = term_2(2:nxg-1,j,:) * (V_(2:nxg-1,j,:)-V_(2:nxg-1,j-1,:))
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,6) 

    Do k = 2,nzg
      term(:,1:ny,k) = (V_(:,:,k)-V_(:,:,k-1))
    End Do
    Call interpolate_y(term  (:,1:ny,2:nzg),term_1(:,2:ny,2:nzg  ),1)
    Call interpolate_z(term_1(:,2:ny,2:nzg),term_2(:,2:ny,2:nzg-1),1) 
    Do k = 2,nz
      term_1(:,2:nyg-1,k) = term_2(:,2:nyg-1,k) * (W_(:,2:nyg-1,k)-W_(:,2:nyg-1,k-1))
    End Do
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,6) 

    Do i = 2,nxg
      term(i,1:ny,:) = (V_(i,:,:)-V_(i-1,:,:))
    End Do
    Call interpolate_y(term(2:nxg,1:ny,:),term_1(2:nxg,2:ny,:),1)
    Do i = 2,nxg
      term(i,:,1:nz) = (W_(i,:,:)-W_(i-1,:,:))
    End Do
    Call interpolate_z(term(2:nxg,:,1:nz),term_2(2:nxg,:,2:nz),1)
    term_3 = term_1*term_2
    Call interpolate_x(term_3(2:nxg,2:nyg-1,2:nzg-1),term_1(2:nxg-1,2:nyg-1,2:nzg-1),1)
    nu_t_(2:nxg,2:nyg,2:nzg) = nu_t_(2:nxg,2:nyg,2:nzg) + 2d0 * term_1(2:nxg,2:nyg,2:nzg) * Sij(:,:,:,6) 

    ! nu_t_ = numerator / denominator
    nu_t_ = - C * nu_t_ / term_4
    ! clipping negative values
    nu_t_ = Max( nu_t_, 0d0 )

    Call apply_periodic_bc_x(nu_t_,2)
    Call update_ghost_interior_planes(nu_t_,4)
    Call apply_periodic_bc_z(nu_t_,4)
 
    !!!!!!!!!!!!!!!!!!!!!
!    nu_t_(:,    2,:) = nu_t_(:,    3,:)
!    nu_t_(:,nyg-1,:) = nu_t_(:,nyg-2,:)
    !!!!!!!!!!!!!!!!!!!!!

    ! boundary conditions ghost cell (must be done after clipping)
    If      ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(nu_t_,2)
    Elseif  ( Dirichlet_nu_t == 0 ) Then
       Call apply_Neumann_bc_y  (nu_t_,2)
    Elseif  ( Dirichlet_nu_t == 2 ) Then
       nu_t_(:,  1,:) = -nu_t_(:,    2,:) + 2d0*0.38d0*dPdx**0.5d0*(y(2)-y(1))/2d0
       nu_t_(:,nyg,:) = -nu_t_(:,nyg-1,:) + 2d0*0.38d0*dPdx**0.5d0*(y(2)-y(1))/2d0
    End If


  End Subroutine sgs_AMD

  !------------------------------------------------------------!
  !                                                            !
  !  Compute Vreman model                                      ! 
  !                                                            !
  ! nu_t_  = C * sqrt( B_b / a_ij * a_ij )                     !
  ! C      = 2.5*C_s^2, C_s Smagorinsky constant               !
  ! a_ij   = d u_j / d x_i                                     ! 
  ! b_ij   = delta_m^2 * a_mj * a_mi                           ! 
  ! B_b    = b_11 * b_22 - b_12^2 + b_11*b_33 - b_13^2         !
  !          + b_22*b_33 - b_23^2                              ! 
  !                                                            !
  ! Tensors are organized in arrays as                         !
  !   ( 1 4 5 )                                                !
  !   ( 7 2 6 )                                                !
  !   ( 8 9 3 )                                                !
  ! where the number denotes the 4th component of the array    !
  !                                                            !
  ! Input:  U_,V_,W_ (velocities)                              !
  ! Output: nu_t_,                                             !
  !                                                            !
  !------------------------------------------------------------!
  Subroutine sgs_Vreman(U_,V_,W_,nu_t_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_

    Real(Int64), Dimension(nxg,nyg,nzg), Intent(Out) :: nu_t_

    Real(Int64)    :: C
    Integer(Int32) :: i, j, k
    Real(Int64), Dimension(nxg,nyg,nzg,9) :: a_ij 
    Real(Int64), Dimension(nxg,nyg,nzg,6) :: b_ij 
    Real(Int64), Dimension(nxg,nyg,nzg  ) :: B_b 
    Real(Int64), Dimension(nxg,nyg,nzg  ) :: A_b 

    ! Vreman constant
    C = 2.5d0*(0.1d0**2d0)
  
    ! Compute a_ij
    Do i = 2,nxg-1
       a_ij(i,:,:,1)       = ( U_(i,:,:) - U_(i-1,:,:) ) / (x(i) - x(i-1)) 
       a_ij(i,2:nyg-1,:,4) = 0.25d0*( V_(i,2:ny  ,:) - V_(i-1,2:ny  ,:) ) / (xg(i) - xg(i-1)) + &
                             0.25d0*( V_(i,1:ny-1,:) - V_(i-1,1:ny-1,:) ) / (xg(i) - xg(i-1)) + &
                             0.25d0*( V_(i+1,2:ny  ,:) - V_(i,2:ny  ,:) ) / (xg(i+1) - xg(i)) + &
                             0.25d0*( V_(i+1,1:ny-1,:) - V_(i,1:ny-1,:) ) / (xg(i+1) - xg(i)) 
       a_ij(i,:,2:nzg-1,5) = 0.25d0*( W_(i,:,2:nz  ) - W_(i-1,:,2:nz  ) ) / (xg(i) - xg(i-1)) + &
                             0.25d0*( W_(i,:,1:nz-1) - W_(i-1,:,1:nz-1) ) / (xg(i) - xg(i-1)) + &
                             0.25d0*( W_(i+1,:,2:nz  ) - W_(i,:,2:nz  ) ) / (xg(i+1) - xg(i)) + &
                             0.25d0*( W_(i+1,:,1:nz-1) - W_(i,:,1:nz-1) ) / (xg(i+1) - xg(i)) 
    End Do
    Do j = 2,nyg-1
       a_ij(:,j,:,2)       = ( V_(:,j,:) - V_(:,j-1,:) ) / (y(j) - y(j-1)) 
       a_ij(2:nxg-1,j,:,7) = 0.25d0*( U_(2:nx  ,j,:) - U_(2:nx  ,j-1,:) ) / (yg(j) - yg(j-1)) + &
                             0.25d0*( U_(1:nx-1,j,:) - U_(1:nx-1,j-1,:) ) / (yg(j) - yg(j-1)) + &
                             0.25d0*( U_(2:nx  ,j+1,:) - U_(2:ny  ,j,:) ) / (yg(j+1) - yg(j)) + &
                             0.25d0*( U_(1:nx-1,j+1,:) - U_(1:ny-1,j,:) ) / (yg(j+1) - yg(j)) 
       a_ij(:,j,2:nzg-1,6) = 0.25d0*( W_(:,j,2:nz  ) - W_(:,j-1,2:nz  ) ) / (yg(j) - yg(j-1)) + &
                             0.25d0*( W_(:,j,1:nz-1) - W_(:,j-1,1:nz-1) ) / (yg(j) - yg(j-1)) + &
                             0.25d0*( W_(:,j+1,2:nz  ) - W_(:,j,2:nz  ) ) / (yg(j+1) - yg(j)) + &
                             0.25d0*( W_(:,j+1,1:nz-1) - W_(:,j,1:nz-1) ) / (yg(j+1) - yg(j)) 
    End Do
    Do k = 2,nzg-1
       a_ij(:,:,k,3)       = ( W_(:,:,k) - W_(:,:,k-1) ) / (z(k) - z(k-1)) 
       a_ij(2:nxg-1,:,k,8) = 0.25d0*( U_(2:nx  ,:,k) - U_(2:nx  ,:,k-1) ) / (zg(k) - zg(k-1)) + &
                             0.25d0*( U_(1:nx-1,:,k) - U_(1:nx-1,:,k-1) ) / (zg(k) - zg(k-1)) + &
                             0.25d0*( U_(2:nx  ,:,k+1) - U_(2:ny  ,:,k) ) / (zg(k+1) - zg(k)) + &
                             0.25d0*( U_(1:nx-1,:,k+1) - U_(1:ny-1,:,k) ) / (zg(k+1) - zg(k)) 
       a_ij(:,2:nyg-1,k,9) = 0.25d0*( V_(:,2:ny  ,k) - V_(:,2:ny  ,k-1) ) / (zg(k) - zg(k-1)) + &
                             0.25d0*( V_(:,1:ny-1,k) - V_(:,1:ny-1,k-1) ) / (zg(k) - zg(k-1)) + &
                             0.25d0*( V_(:,2:ny  ,k+1) - V_(:,2:ny  ,k) ) / (zg(k+1) - zg(k)) + &
                             0.25d0*( V_(:,1:ny-1,k+1) - V_(:,1:ny-1,k) ) / (zg(k+1) - zg(k)) 
    End Do

    b_ij = 0d0
    ! Compute b_ij
    Do i = 2,nxg-1
      b_ij(i,:,:,1) = b_ij(i,:,:,1) + (x(i)-x(i-1))**2d0 * a_ij(i,:,:,1) * a_ij(i,:,:,1)   
      b_ij(i,:,:,2) = b_ij(i,:,:,2) + (x(i)-x(i-1))**2d0 * a_ij(i,:,:,4) * a_ij(i,:,:,4)   
      b_ij(i,:,:,3) = b_ij(i,:,:,3) + (x(i)-x(i-1))**2d0 * a_ij(i,:,:,5) * a_ij(i,:,:,5)   
      b_ij(i,:,:,4) = b_ij(i,:,:,4) + (x(i)-x(i-1))**2d0 * a_ij(i,:,:,1) * a_ij(i,:,:,4)   
      b_ij(i,:,:,5) = b_ij(i,:,:,5) + (x(i)-x(i-1))**2d0 * a_ij(i,:,:,1) * a_ij(i,:,:,5)   
      b_ij(i,:,:,6) = b_ij(i,:,:,6) + (x(i)-x(i-1))**2d0 * a_ij(i,:,:,4) * a_ij(i,:,:,5)   
    End Do

    Do j = 2,nyg-1
      b_ij(:,j,:,1) = b_ij(:,j,:,1) + (y(j)-y(j-1))**2d0 * a_ij(:,j,:,7) * a_ij(:,j,:,7)   
      b_ij(:,j,:,2) = b_ij(:,j,:,2) + (y(j)-y(j-1))**2d0 * a_ij(:,j,:,2) * a_ij(:,j,:,2)   
      b_ij(:,j,:,3) = b_ij(:,j,:,3) + (y(j)-y(j-1))**2d0 * a_ij(:,j,:,6) * a_ij(:,j,:,6)   
      b_ij(:,j,:,4) = b_ij(:,j,:,4) + (y(j)-y(j-1))**2d0 * a_ij(:,j,:,7) * a_ij(:,j,:,2)   
      b_ij(:,j,:,5) = b_ij(:,j,:,5) + (y(j)-y(j-1))**2d0 * a_ij(:,j,:,7) * a_ij(:,j,:,6)   
      b_ij(:,j,:,6) = b_ij(:,j,:,6) + (y(j)-y(j-1))**2d0 * a_ij(:,j,:,2) * a_ij(:,j,:,6)   
    End Do

    Do k = 2,nzg-1
      b_ij(:,:,k,1) = b_ij(:,:,k,1) + (z(k)-z(k-1))**2d0 * a_ij(:,:,k,8) * a_ij(:,:,k,8)   
      b_ij(:,:,k,2) = b_ij(:,:,k,2) + (z(k)-z(k-1))**2d0 * a_ij(:,:,k,9) * a_ij(:,:,k,9)   
      b_ij(:,:,k,3) = b_ij(:,:,k,3) + (z(k)-z(k-1))**2d0 * a_ij(:,:,k,3) * a_ij(:,:,k,3)   
      b_ij(:,:,k,4) = b_ij(:,:,k,4) + (z(k)-z(k-1))**2d0 * a_ij(:,:,k,8) * a_ij(:,:,k,9)   
      b_ij(:,:,k,5) = b_ij(:,:,k,5) + (z(k)-z(k-1))**2d0 * a_ij(:,:,k,8) * a_ij(:,:,k,3)   
      b_ij(:,:,k,6) = b_ij(:,:,k,6) + (z(k)-z(k-1))**2d0 * a_ij(:,:,k,9) * a_ij(:,:,k,3)   
    End Do

    ! Compute B_b
    B_b = b_ij(:,:,:,1) * b_ij(:,:,:,2) - b_ij(:,:,:,4)**2d0 + &
          b_ij(:,:,:,1) * b_ij(:,:,:,3) - b_ij(:,:,:,5)**2d0 + &
          b_ij(:,:,:,2) * b_ij(:,:,:,3) - b_ij(:,:,:,6)**2d0

    A_b =  a_ij(:,:,:,1)**2d0 + a_ij(:,:,:,2)**2d0 + a_ij(:,:,:,3)**2d0 + & 
           a_ij(:,:,:,4)**2d0 + a_ij(:,:,:,5)**2d0 + a_ij(:,:,:,6)**2d0 + & 
           a_ij(:,:,:,7)**2d0 + a_ij(:,:,:,8)**2d0 + a_ij(:,:,:,9)**2d0 
    
    nu_t_ = C * (B_b / A_b)**0.5d0 

    Call apply_periodic_bc_x(nu_t_,2)
    Call update_ghost_interior_planes(nu_t_,4)
    Call apply_periodic_bc_z(nu_t_,4)

  
    If      ( Dirichlet_nu_t == 1 ) Then
       Call apply_Dirichlet_bc_y(nu_t_,2)
    Elseif  ( Dirichlet_nu_t == 0 ) Then
       Call apply_Neumann_bc_y  (nu_t_,2)
    Elseif  ( Dirichlet_nu_t == 2 ) Then
       nu_t_(:,  1,:) = -nu_t_(:,    2,:) + 2d0*0.38d0*dPdx**0.5d0*(y(2)-y(1))/2d0
       nu_t_(:,nyg,:) = -nu_t_(:,nyg-1,:) + 2d0*0.38d0*dPdx**0.5d0*(y(2)-y(1))/2d0
    End If

  End Subroutine sgs_Vreman

  !-----------------------------------------------------------!
  !                                                           !
  !           Compute rate-of-strain at cell centers          !
  !                                                           !
  ! Input:  U, V, W, (flow velocities)                        !
  ! Output: Sij, S=sqrt(2Sij*Sij)                             !
  !                                                           !
  !-----------------------------------------------------------!
  Subroutine compute_Sij(U_,V_,W_,Sij_,S_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_    

    Real(Int64), Dimension(2:nxg,  2:nyg,  2:nzg,  6), Intent(Out) :: Sij_
    Real(Int64), Dimension(2:nxg-1,2:nyg-1,2:nzg-1),   Intent(Out) :: S_

    ! local variables
    Integer(Int32) :: i, j, k

    ! Compute diagonal Sij
    Do k = 2,nzg-1
      Do j = 2,nyg-1
         Do i = 2,nxg-1
          ! These are at cell center
          Sij_(i,j,k,1) = (U_(i,j,k) - U_(i-1,j,k))/(x(i)-x(i-1)) ! dU/dx
          Sij_(i,j,k,2) = (V_(i,j,k) - V_(i,j-1,k))/(y(j)-y(j-1)) ! dV/dy
          Sij_(i,j,k,3) = (W_(i,j,k) - W_(i,j,k-1))/(z(k)-z(k-1)) ! dW/dz
        End Do
      End Do
    End Do

    ! Compute off-diagonal Sij
    Do k = 2,nzg 
      Do j = 2,nyg
         Do i = 2,nxg
          ! These are at cell edges
          Sij_(i,j,k,4) = 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j-1,k))/(yg(j)-yg(j-1)) + (V_(i,j-1,k)-V_(i-1,j-1,k))/(xg(i)-xg(i-1))) ! 1/2*(dU/dy + dV/dx)
          Sij_(i,j,k,5) = 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j,k-1))/(zg(k)-zg(k-1)) + (W_(i,j,k-1)-W_(i-1,j,k-1))/(xg(i)-xg(i-1))) ! 1/2*(dU/dz + dW/dx)
          Sij_(i,j,k,6) = 0.5d0 * ((V_(i,j-1,k)-V_(i,j-1,k-1))/(zg(k)-zg(k-1)) + (W_(i,j,k-1)-W_(i,j-1,k-1))/(yg(j)-yg(j-1))) ! 1/2*(dV/dz + dW/dy)
        End Do
      End Do
    End Do

    ! Move values from cell edge to cell center
    Do k = 2,nzg-1
      Do j = 2,nyg-1
         Do i = 2,nxg-1       
          Sij_(i,j,k,4) = 0.25d0 * (Sij_(i,j,k,4) + Sij_(i+1,j,k,4) + Sij_(i,j+1,k,4) + Sij_(i+1,j+1,k,4)) 
          Sij_(i,j,k,5) = 0.25d0 * (Sij_(i,j,k,5) + Sij_(i+1,j,k,5) + Sij_(i,j,k+1,5) + Sij_(i+1,j,k+1,5))
          Sij_(i,j,k,6) = 0.25d0 * (Sij_(i,j,k,6) + Sij_(i,j+1,k,6) + Sij_(i,j,k+1,6) + Sij_(i,j+1,k+1,6))
        End Do
      End Do
    End Do

    ! Compute |S| = sqrt(2*Sij*Sij) (at cell centers)
    S_ = 0d0
    Do i = 1,6
      If (i .le. 3) Then
        S_ = S_ + 2d0 * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i) * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i)
      Else 
        S_ = S_ + 4d0 * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i) * Sij_(2:nxg-1,2:nyg-1,2:nzg-1,i)
      End If
    End Do
    S_ = S_**0.5d0

  End Subroutine compute_Sij

  !------------------------------------------------------------------------------!
  !                                                                              !
  !                    Compute rate-of-strain at V locations                     !
  !                     for the first two cell at the wall                       !
  !                                                                              !
  ! Assumed uniform mesh in y                                                    !
  !                                                                              !
  ! Input:  U, V, W, (flow velocities)                                           !
  ! Output: Sij(2:nxg,2:3,2:nzg) (bottom) and Sij(2:nxg,nyg-1:nyg,2:nzg) (top)   !
  !                                                                              !
  !------------------------------------------------------------------------------!
  Subroutine compute_Sij_at_V_location_wall(U_,V_,W_,Sij_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_    

    Real(Int64), Dimension(2:nxg,2:nyg,2:nzg,6), Intent(Out) :: Sij_

    ! local variables
    Integer(Int32) :: i, j, k, jv
    Integer(Int32) :: jindex(6), jj

    ! indices for diagonal part
    jindex(1) = 2
    jindex(2) = 3
    jindex(3) = 4
    jindex(4) = nyg-2
    jindex(5) = nyg-1
    jindex(6) = nyg

    ! Compute diagonal Sij
    Do k = 2, nzg-1
       Do jj = 1, 6
          j  = jindex(jj)
          Do i = 2, nxg-1
             ! These are at cell center
             ! dU/dx
             Sij_(i,j,k,1) = 0.5d0*( (U_(i,j,k) - U_(i-1,j,k))/(x(i)-x(i-1)) + (U_(i,j-1,k) - U_(i-1,j-1,k))/(x(i)-x(i-1)) )
             ! dV/dy first order (this could be 2nd order for jindex(2:3))
             jv = j
             If ( j>nyg/2 ) jv = j-1

             If ( j==2 .or. j==nyg) Then
                Sij_(i,j,k,2) = (V_(i,jv,k) - V_(i,jv-1,k))/(y(jv)-y(jv-1)) 
             Else If ( j < nyg/2 ) Then 
                Sij_(i,j,k,2) = (V_(i,jv,k) - V_(i,jv-2,k))/(y(jv+1)-y(jv-1)) 
             Else 
                Sij_(i,j,k,2) = (V_(i,jv+1,k) - V_(i,jv-1,k))/(y(jv+1)-y(jv-1)) 
             End If
             ! dW/dz
             Sij_(i,j,k,3) = 0.5d0*( (W_(i,j,k) - W_(i,j,k-1))/(z(k)-z(k-1)) + (W_(i,j-1,k) - W_(i,j-1,k-1))/(z(k)-z(k-1)) )
          End Do
       End Do
    End Do

    ! indices for off-diagonal part
    jindex(1) = 2
    jindex(2) = 3
    jindex(3) = 4
    jindex(4) = nyg-2
    jindex(5) = nyg-1
    jindex(6) = nyg

    ! Compute off-diagonal Sij
    Do k = 2, nzg-1
       Do jj = 1, 6
          j  = jindex(jj)
          jv = j
          !If ( j>nyg/2 ) jv = j-1
          Do i = 2, nxg-1
             ! These are at cell edges
             ! 1/2*(dU/dy + dV/dx)
                                     !(i,j)
             Sij_(i,j,k,4) = 0.5d0*( 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j-1,k))/(yg(j)-yg(j-1)) + &
                     (V_(i,jv-1,k)-V_(i-1,jv-1,k))/(xg(i)-xg(i-1))) + & 
                     !(i+1,j)
                     0.5d0 * ((U_(i-1+1,j,k)-U_(i-1+1,j-1,k))/(yg(j)-yg(j-1)) + &
                             (V_(i+1,jv-1,k)-V_(i-1+1,jv-1,k))/(xg(i+1)-xg(i-1+1))) )   

             ! 1/2*(dU/dz + dW/dx), average over planes j=1:2
             Sij_(i,j,k,5) = 1d0/8d0*( &
                  !-------------j
                  !(i,k)
                  0.5d0 * ((U_(i-1,j,k)-U_(i-1,j,k-1))/(zg(k)-zg(k-1)) + &
                          (W_(i,j,k-1)-W_(i-1,j,k-1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k)
                  0.5d0 * ((U_(i-1+1,j,k)-U_(i-1+1,j,k-1))/(zg(k)-zg(k-1)) + &
                          (W_(i+1,j,k-1)-W_(i-1+1,j,k-1))/(xg(i+1)-xg(i-1+1))) + &  
                  !(i,k+1)
                  0.5d0 * ((U_(i-1,j,k+1)-U_(i-1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) + &
                          (W_(i,j,k-1+1)-W_(i-1,j,k-1+1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k+1)
                  0.5d0 * ((U_(i-1+1,j,k+1)-U_(i-1+1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) + &
                          (W_(i+1,j,k-1+1)-W_(i-1+1,j,k-1+1))/(xg(i+1)-xg(i-1+1))) + &   
                  !-------------j-1
                  !(i,k)
                  0.5d0 * ((U_(i-1,j-1,k)-U_(i-1,j-1,k-1))/(zg(k)-zg(k-1)) + &
                          (W_(i,j-1,k-1)-W_(i-1,j-1,k-1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k)
                  0.5d0 * ((U_(i-1+1,j-1,k)-U_(i-1+1,j-1,k-1))/(zg(k)-zg(k-1)) + &
                          (W_(i+1,j-1,k-1)-W_(i-1+1,j-1,k-1))/(xg(i+1)-xg(i-1+1))) + &  
                  !(i,k+1)
                  0.5d0 * ((U_(i-1,j-1,k+1)-U_(i-1,j-1,k-1+1))/(zg(k+1)-zg(k-1+1)) + &
                          (W_(i,j-1,k-1+1)-W_(i-1,j-1,k-1+1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k+1)
                  0.5d0 * ((U_(i-1+1,j-1,k+1)-U_(i-1+1,j-1,k-1+1))/(zg(k+1)-zg(k-1+1)) + &
                          (W_(i+1,j-1,k-1+1)-W_(i-1+1,j-1,k-1+1))/(xg(i+1)-xg(i-1+1))) )

             ! 1/2*(dV/dz + dW/dy) 
                                     !(j,k)
             Sij_(i,j,k,6) = 0.5d0*( 0.5d0 * ((V_(i,jv-1,k)-V_(i,jv-1,k-1))/(zg(k)-zg(k-1)) + &
                     (W_(i,j,k-1)-W_(i,j-1,k-1))/(yg(j)-yg(j-1))) + & 
                                     !(j,k+1)
                                     0.5d0 * ((V_(i,jv-1,k+1)-V_(i,jv-1,k-1+1))/(zg(k+1)-zg(k-1+1)) + &
                                             (W_(i,j,k-1+1)-W_(i,j-1,k-1+1))/(yg(j)-yg(j-1))) ) 
          End Do
       End Do
    End Do

    ten_buf(2:nxg,2:nyg,2:nzg,:) = Sij
    Do i = 1,6
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do
    Sij = ten_buf(2:nxg,2:nyg,2:nzg,:)
    
  End Subroutine compute_Sij_at_V_location_wall
  
  !------------------------------------------------------------------------------!
  !                                                                              !
  !                    Compute rate-of-rotation at V locations                   !
  !                     for the first two cell at the wall                       !
  !                                                                              !
  ! Assumed uniform mesh in y                                                    !
  !                                                                              !
  ! Input:  U, V, W, (flow velocities)                                           !
  ! Output: Rij(2:nxg,2:3,2:nzg) (bottom) and Rij(2:nxg,nyg-1:nyg,2:nzg) (top)   !
  !                                                                              !
  !------------------------------------------------------------------------------!
  Subroutine compute_Rij_at_V_location_wall(U_,V_,W_,Rij_)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_    

    Real(Int64), Dimension(2:nxg,2:nyg,2:nzg,6), Intent(Out) :: Rij_

    ! local variables
    Integer(Int32) :: i, j, k, jv
    Integer(Int32) :: jindex(6), jj

    Rij_ = 0d0

    ! indices for off-diagonal part
    jindex(1) = 2
    jindex(2) = 3
    jindex(3) = 4
    jindex(4) = nyg-2
    jindex(5) = nyg-1
    jindex(6) = nyg

    ! Compute off-diagonal Sij
    Do k = 2, nzg-1
       Do jj = 1, 6
          j  = jindex(jj)
          jv = j
          !If ( j>nyg/2 ) jv = j-1
          Do i = 2, nxg-1
             ! These are at cell edges
             ! 1/2*(dU/dy - dV/dx)
                                     !(i,j)
             Rij_(i,j,k,4) = 0.5d0*( 0.5d0 * ((U_(i-1,j,k)-U_(i-1,j-1,k))/(yg(j)-yg(j-1)) - &
                     (V_(i,jv-1,k)-V_(i-1,jv-1,k))/(xg(i)-xg(i-1))) + & 
                                     !(i+1,j)
                                     0.5d0 * ((U_(i-1+1,j,k)-U_(i-1+1,j-1,k))/(yg(j)-yg(j-1)) - &
                                             (V_(i+1,jv-1,k)-V_(i-1+1,jv-1,k))/(xg(i+1)-xg(i-1+1))) )   

             ! 1/2*(dU/dz - dW/dx), average over planes j=1:2
             Rij_(i,j,k,5) = 1d0/8d0*( &
                  !-------------j
                  !(i,k)
                  0.5d0 * ((U_(i-1,j,k)-U_(i-1,j,k-1))/(zg(k)-zg(k-1)) - &
                          (W_(i,j,k-1)-W_(i-1,j,k-1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k)
                  0.5d0 * ((U_(i-1+1,j,k)-U_(i-1+1,j,k-1))/(zg(k)-zg(k-1)) - &
                          (W_(i+1,j,k-1)-W_(i-1+1,j,k-1))/(xg(i+1)-xg(i-1+1))) + &  
                  !(i,k+1)
                  0.5d0 * ((U_(i-1,j,k+1)-U_(i-1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) - &
                          (W_(i,j,k-1+1)-W_(i-1,j,k-1+1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k+1)
                  0.5d0 * ((U_(i-1+1,j,k+1)-U_(i-1+1,j,k-1+1))/(zg(k+1)-zg(k-1+1)) - &
                          (W_(i+1,j,k-1+1)-W_(i-1+1,j,k-1+1))/(xg(i+1)-xg(i-1+1))) + &   
                  !-------------j-1
                  !(i,k)
                  0.5d0 * ((U_(i-1,j-1,k)-U_(i-1,j-1,k-1))/(zg(k)-zg(k-1)) - &
                          (W_(i,j-1,k-1)-W_(i-1,j-1,k-1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k)
                  0.5d0 * ((U_(i-1+1,j-1,k)-U_(i-1+1,j-1,k-1))/(zg(k)-zg(k-1)) - &
                          (W_(i+1,j-1,k-1)-W_(i-1+1,j-1,k-1))/(xg(i+1)-xg(i-1+1))) + &  
                  !(i,k+1)
                  0.5d0 * ((U_(i-1,j-1,k+1)-U_(i-1,j-1,k-1+1))/(zg(k+1)-zg(k-1+1)) - &
                          (W_(i,j-1,k-1+1)-W_(i-1,j-1,k-1+1))/(xg(i)-xg(i-1))) + &  
                  !(i+1,k+1)
                  0.5d0 * ((U_(i-1+1,j-1,k+1)-U_(i-1+1,j-1,k-1+1))/(zg(k+1)-zg(k-1+1)) - &
                          (W_(i+1,j-1,k-1+1)-W_(i-1+1,j-1,k-1+1))/(xg(i+1)-xg(i-1+1))) )

             ! 1/2*(dV/dz - dW/dy) 
                                     !(j,k)
             Rij_(i,j,k,6) = 0.5d0*( 0.5d0 * ((V_(i,jv-1,k)-V_(i,jv-1,k-1))/(zg(k)-zg(k-1)) - &
                     (W_(i,j,k-1)-W_(i,j-1,k-1))/(yg(j)-yg(j-1))) + & 
                                     !(j,k+1)
                                     0.5d0 * ((V_(i,jv-1,k+1)-V_(i,jv-1,k-1+1))/(zg(k+1)-zg(k-1+1)) - &
                                             (W_(i,j,k-1+1)-W_(i,j-1,k-1+1))/(yg(j)-yg(j-1))) ) 
          End Do
       End Do
    End Do

    ten_buf(2:nxg,2:nyg,2:nzg,:) = Rij
    Do i = 1,6
       Call apply_periodic_bc_x(ten_buf(:,:,:,i),2)
       Call apply_periodic_bc_z(ten_buf(:,:,:,i),4)
       Call update_ghost_interior_planes(ten_buf(:,:,:,i),4)
    End Do
    Rij = ten_buf(2:nxg,2:nyg,2:nzg,:)
    
  End Subroutine compute_Rij_at_V_location_wall

End Module subgrid
