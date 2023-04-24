!-----------------------------------------!
! Module with all shared global variables !
!-----------------------------------------!
Module global

  ! General Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use, Intrinsic :: iso_c_binding

  ! prevent implicit typing
  Implicit None

  ! FFTW
  Include 'fftw3-mpi.f03'  

  !------------------Declarations-----------------!

  ! MPI Communicator
  Integer(Int32) :: comm

  ! step number
  Integer(Int32) :: istep, rk_step
  Real   (Int64) :: time1, time2

  ! constants
  Real(Int64) :: pi = 4d0*datan(1d0)  

  ! files
  Character(200) :: filein, fileout, filepar
  Integer(Int32) :: nsave, nmonitor

  ! set random initial condition
  Integer(Int32) :: random_init

  ! domain size
  Real(Int64) :: Lx, Lz, Ly, Lxp, Lzp

  ! steps
  Integer(Int32) :: nsteps, nstep_init
  Real   (Int64) :: dt, t

  ! viscosity
  Real(Int64) :: nu, delta

  ! global face points
  Integer(Int32) :: nx_global, ny_global, nz_global

  ! local face points
  Integer(Int32) :: nx, ny, nz

  ! global center points
  Integer(Int32) :: nxm_global, nym_global, nzm_global

  ! local center points
  Integer(Int32) :: nxm, nym, nzm

  ! global center points + ghost cells
  Integer(Int32) :: nxg_global, nyg_global, nzg_global

  ! local center points + ghost cells
  Integer(Int32) :: nxg, nyg, nzg

  ! global grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x_global, y_global, z_global

  ! local grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x, y, z

  ! global grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm_global, ym_global, zm_global

  ! local grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm, ym, zm

  ! global grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg_global, yg_global, zg_global

  ! local grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg, yg, zg

  ! middle points for yg->yg_m and yg_m->yg_mm
  Real(Int64), Allocatable, Dimension(:) :: yg_m, yg_mm

  ! local velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: U,V,W,P
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uo,Vo,Wo,Po
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uoo,Voo,Woo,Poo

  ! local auxiliary 
  Real(Int64), Allocatable, Dimension(:,:,:) :: term_1, term_2, term_3, term_4, term

  ! local rhs for velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:)   :: px_bottom, px_top
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_p
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uo, rhs_vo, rhs_wo
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uf, rhs_vf, rhs_wf

  ! local rhs for pressure in Fourier
  Complex(Int64), Dimension(:,:,:), Allocatable :: rhs_p_hat
  Complex(Int64), Dimension(:),     Allocatable :: rhs_aux

  ! local auxiliary arrays for MPI_sendrev boundary conditions
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ui, buffer_vi, buffer_ci, buffer_ai
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ue, buffer_ve, buffer_we, buffer_wi, buffer_ce, buffer_ae
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_p

  ! local auxiliary arrays for MPI_sendrev interior planes
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_us, buffer_ur
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_vs, buffer_vr
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_ws, buffer_wr
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_as, buffer_ar
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_ps, buffer_pr
  
  ! local auxiliary planes for FFTW
  Type(C_PTR) :: cplane_fft
  Complex(C_DOUBLE_COMPLEX), Pointer, Dimension(:,:) :: plane, plane_hat

  ! Fourier points and wave numbers 
  Integer(C_INTPTR_T) :: nxp_global, nzp_global, local_k_offset
  Integer(C_INTPTR_T) :: nxp, nzp
  Integer(C_INTPTR_T) :: mx_global, mz_global
  Integer(C_INTPTR_T) :: mx, mz
  Real   (Int64)      :: dx, dz
  Real   (Int64), Dimension(:), Allocatable :: kxx, kzz

  ! Mappings for fft modes
  Integer(Int64), Dimension(:),   Allocatable :: kmode_map
  Integer(Int64), Dimension(:,:), Allocatable :: imode_map_fft, kmode_map_fft
  
  ! FFTW plans
  Integer(C_INTPTR_T) :: alloc_local
  Type   (C_PTR)      :: plan_d, plan_i

  ! finite differences (second derivative)
  Real(Int64) :: ddx1, ddx2, ddx3
  Real(Int64) :: ddy1, ddy2, ddy3
  Real(Int64) :: ddz1, ddz2, ddz3

  ! linear solver
  Integer (Int32) :: nr, nrhs
  Integer (Int32), Dimension(:),   Allocatable :: pivot  
  Complex (Int64), Dimension(:),   Allocatable :: D, DL, DU
  Complex (Int64), Dimension(:,:), Allocatable :: M, Dyy

  ! pressure gradients
  Real(Int64) :: dPdx, dPdx_ref

  ! constant mass flow
  Real   (Int64) :: Qflow_x_0
  Integer(Int32) :: x_mass_cte
    
  ! CFL parameters
  Real(Int64) :: CFL, dxmin, dymin, dzmin

  ! interpolation weights 
  Integer(Int32) :: in1, in2
  Real(Int64), Dimension(:), Allocatable :: weight_y_0, weight_y_1

  ! actual pressure boundary conditions
  Real   (Int64) :: coef_bc_1, coef_bc_2
  Real   (Int64), Dimension(:,:), Allocatable :: bc_1,     bc_2
  Complex(Int64), Dimension(:,:), Allocatable :: bc_1_hat, bc_2_hat
  Logical(Int32) :: pressure_computed

  ! statistics
  Integer(Int32) :: nstats
  Real   (Int64) :: Retau, utau, Qflow_x, Qflow_y
  Real   (Int64), Dimension(:), Allocatable ::  Umean,  Vmean,  Wmean, Pmean
  Real   (Int64), Dimension(:), Allocatable :: U2mean, V2mean, W2mean, UVmean, P2mean

  ! Runge-Kutta 3 coefficients and buffers
  Real(Int64), Dimension(:),     Allocatable :: rk_t
  Real(Int64), Dimension(:,:),   Allocatable :: rk_coef
  Real(Int64), Dimension(:,:,:), Allocatable :: Fu1, Fu2, Fu3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fv1, Fv2, Fv3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fw1, Fw2, Fw3

  ! sgs model
  Integer(Int32) :: LES_model, Dirichlet_nu_t
  Real   (Int64), Allocatable, Dimension(:,:,:,:) :: Lij, Mij, Sij, Rij, ten_buf, ten_buf_2
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: Uf, Vf, Wf, S
  Real   (Int64), Allocatable, Dimension(:,:,:)   :: nu_t, nu_t_hat
  Real   (Int64) :: fil_size, fil_size_2

  ! wall-model
  Real   (Int64), Allocatable, Dimension(:,:,:) :: alpha_x, alpha_z

  ! smarites
  integer, parameter :: NUM_ACTIONS = 1
  integer, parameter :: STATE_SIZE  = 6 
  integer, parameter :: NUM_AGENTS  = 256
  integer :: nagents, nx_nagents, nz_nagents
  real :: Ucl

  
  real(Int64), Dimension(:,:,:), Allocatable :: actions
  real(Int64), Dimension(:,:),   Allocatable :: rewards, rewards_old
  real(Int64), Dimension(:,:,:), Allocatable :: states 
  
End Module global
