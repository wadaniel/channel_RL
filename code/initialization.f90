!------------------------------------------------!
! Module with initialization of global variables !
!------------------------------------------------!
Module initialization

  ! Modules
  Use, Intrinsic :: iso_c_binding
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use input_output
  Use mass_flow

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  !----------------------------------------!
  !         Initialize everything          !
  !----------------------------------------!
  Subroutine initialize(in_comm)
    
    Integer(Int32) :: in_comm
    Integer(Int32) :: i, j, k, kk, nzpe, pos, ipos, nze, nzme
    Real   (Int64) :: dy1, dy2, det, a, b, c, r, Qflow_ref
    Integer(Int32), Dimension(:,:), Allocatable :: A_kmodes, A_kmodes_local

    !---------------------first initialize MPI-------------------!
    comm = in_comm

    !call Mpi_init(ierr)
    call Mpi_comm_size(comm, nprocs, ierr)
    call Mpi_comm_rank(comm,   myid, ierr)

    If (myid==0) Then 
       Write(*,*) '----------------------------------------------------------------------'       
       Write(*,*) ' '
       Write(*,*) '           My channel ^^, parallel version 0.25                       '
       Write(*,*) ' '
       Write(*,*) '----------------------------------------------------------------------'
    End If

    !--------------read parameters from standard input-----------!
    If ( myid==0 ) Write(*,*) 'reading input paramaters...'
    Call RANDOM_NUMBER(r)
    if (r < 1d0/3d0) then
       filepar = '../case_folder/input_parameters_2000'
    elseif (r < 2d0/3d0) then
       filepar = '../case_folder/input_parameters_4200'
    else
       filepar = '../case_folder/input_parameters_8000'
    end if
    Call read_input_parameters

    ! time
    delta = 1d0
    t = 0d0
    
    !-------------------grid definitions-------------------------!
    Allocate (  k1_global(0:nprocs-1),  k2_global(0:nprocs-1) )
    Allocate ( kg1_global(0:nprocs-1), kg2_global(0:nprocs-1) )

    ! restrictions for FFTW mapping
    If ( Mod( nx_global   , 2     )/=0 ) Stop 'Error: nx must be even'
    If ( Mod( nz_global   , 2     )/=0 ) Stop 'Error: nz must be even'
    If ( Mod( nz_global-2 , nprocs)/=0 ) Stop 'nz-2 should be divisible by number of processors'

    ! number of interior z-planes per processor based on fftw decomposition
    nslices_z = Nint( Real((nz_global-2))/Real(nprocs) ) 
    
    ! restriction for MPI boundaries
    If ( nslices_z<2 ) Stop 'Error: nslices_z must be at least 2' 

    ! domain decomposition. Must be consistent with fftw
    Do i = 0, nprocs-1
       ! range index for faces in each processor
       k1_global(i)  = i*nslices_z  + 1
       k2_global(i)  = k1_global(i) + nslices_z + 1
       ! range index for centers in each processor
       kg1_global(i) = i*nslices_z   + 1
       kg2_global(i) = kg1_global(i) + nslices_z + 1
    End Do    

    ! remaining planes in last processor
    k2_global (nprocs-1) = nz_global 
    kg2_global(nprocs-1) = nz_global + 1

    ! face points
    nx = nx_global
    ny = ny_global
    nz = k2_global(myid) - k1_global(myid) + 1 

    ! middle points
    nxm_global = nx_global - 1
    nym_global = ny_global - 1
    nzm_global = nz_global - 1

    nxm = nx - 1
    nym = ny - 1
    nzm = kg2_global(myid) - kg1_global(myid) + 1 - 2  
    
    ! middle points + ghost cells
    nxg_global = nxm_global + 2
    nyg_global = nym_global + 2
    nzg_global = nzm_global + 2

    nxg = nxm + 2
    nyg = nym + 2
    nzg = kg2_global(myid) - kg1_global(myid) + 1 

    ! size for last proccesor nz and nzm -> nze and nzme
    nze  = nz
    nzme = nzm
    Call Mpi_bcast (  nze,1,MPI_integer,nprocs-1,comm,ierr )
    Call Mpi_bcast ( nzme,1,MPI_integer,nprocs-1,comm,ierr )
   
    ! Allocate main arrays
    If ( myid==0 ) Write(*,*) 'allocating main arrays...'
    Allocate ( x_global (  nx_global),  y_global (  ny_global),  z_global (  nz_global)  )
    Allocate ( xm_global( nxm_global),  ym_global( nym_global),  zm_global( nzm_global)  )
    Allocate ( xg_global(nxm_global+2), yg_global(nym_global+2), zg_global(nzm_global+2) )

    Allocate (  x (  nx),  y (  ny),  z (  nz) )
    Allocate (  xm( nxm),  ym( nym),  zm( nzm) )
    Allocate ( xg(nxm+2), yg(nym+2), zg(nzm+2) )

    Allocate ( yg_m (nyg-1) )
    Allocate ( yg_mm(nyg-2) )
    
    ! global interior + boundary + ghost points
    Allocate (U (    nx, nym+2, nzm+2) )
    Allocate (V ( nxm+2,    ny, nzm+2) )
    Allocate (W ( nxm+2, nym+2,    nz) )
    Allocate (P ( nxm+2, nym+2, nzm+2) )

    Allocate (Uo  (    nx,  nym+2, nzm+2) )
    Allocate (Vo  ( nxm+2,     ny, nzm+2) )
    Allocate (Wo  ( nxm+2,  nym+2,    nz) )
    Allocate (Po  ( nxm+2,  nym+2, nzm+2) )

    Allocate (Uoo (    nx,  nym+2, nzme+2) ) ! z-planes modified for I/O
    Allocate (Voo ( nxm+2,     ny, nzme+2) )
    Allocate (Woo ( nxm+2,  nym+2,    nze) )
    Allocate (Poo ( nxm+2,  nym+2, nzme+2) )

    ! Auxiliary arrays
    Allocate ( term   ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_1 ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_2 ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_3 ( nxg, nyg, nzm+2 ) ) 
    Allocate ( term_4 ( nxg, nyg, nzm+2 ) ) 

    ! RHS: interior points only
    Allocate ( rhs_uo ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) ) 
    Allocate ( rhs_vo ( 2:nxg-1, 2:ny-1,  2:nzg-1 ) )
    Allocate ( rhs_wo ( 2:nxg-1, 2:nyg-1, 2:nz-1  ) )
    Allocate ( rhs_p  ( 2:nxg-1, 2:nyg-1, 2:nzg   ) ) ! ONE EXTRA PLANE IN Z FOR GHOST CELL
    

    !--------------------------Boundary conditions--------------------------!
    ! local velocity, initial z-planes
    Allocate ( buffer_ui(nx,nyg,2:3), buffer_vi(nxg,ny,2:3), buffer_wi(nxg,nyg), buffer_ci(nxg,nyg,2:3) )
    Allocate ( buffer_ai(nx,2,2:3) )
    ! local velocity, ending  z-planes
    Allocate ( buffer_ue(nx,nyg),     buffer_ve(nxg,ny),     buffer_we(nxg,nyg), buffer_ce(nxg,nyg) )
    Allocate ( buffer_ae(nx,2) )
    ! local pressure z-plane
    Allocate ( buffer_p(2:nxg-1,2:nyg-1) ) 

    !------------------------Interior communications------------------------!
    Allocate ( buffer_us(nx ,nyg), buffer_ur(nx ,nyg) )
    Allocate ( buffer_vs(nxg, ny), buffer_vr(nxg, ny) )
    Allocate ( buffer_ws(nxg,nyg), buffer_wr(nxg,nyg) )
    Allocate ( buffer_as(nx,2), buffer_ar(nx,2) )
    Allocate ( buffer_ps(2:nxg-1,2:nyg-1), buffer_pr(2:nxg-1,2:nyg-1) ) 

    ! read data 
    If ( myid==0 ) Write(*,*) 'preparing initial condition...'
    If ( random_init==1 ) Then
       Call init_flow
    Else
       Call read_input_data
    End If

    ! definie global grids from x_global, y_global and z_global (face to centers)
    ! local faces
    x = x_global
    y = y_global
    z = z_global( k1_global(myid):k2_global(myid) )

    ! global interior centers
    Do i = 1, nxm_global
       xm_global(i) = 0.5d0*( x_global(i) + x_global(i+1) )
    End Do
    Do j=1,nym_global
       ym_global(j) = 0.5d0*( y_global(j) + y_global(j+1) )
    End Do
    Do k=1,nzm_global
       zm_global(k) = 0.5d0*( z_global(k) + z_global(k+1) )
    End Do

    ! local interior centers
    xm = xm_global
    ym = ym_global
    zm = zm_global( kg1_global(myid):kg2_global(myid)-2 )

    ! global 
    xg_global(2:nxm_global+1) = xm_global
    xg_global(1)              = xm_global(1)          - 2d0*(xm_global(1)-x_global(1))
    xg_global(nxm_global+2)   = xm_global(nxm_global) + 2d0*(x_global(nx_global)-xm_global(nxm_global))
    
    yg_global(2:nym_global+1) = ym_global
    yg_global(1)              = ym_global(1)          - 2d0*(ym_global(1)-y_global(1))
    yg_global(nym_global+2)   = ym_global(nym_global) + 2d0*(y_global(ny_global)-ym_global(nym_global))
    
    zg_global(2:nzm_global+1) = zm_global
    zg_global(1)              = zm_global(1)          - 2d0*(zm_global(1)-z_global(1))
    zg_global(nzm_global+2)   = zm_global(nzm_global) + 2d0*(z_global(nz_global)-zm_global(nzm_global))    

    xg = xg_global
    yg = yg_global
    zg = zg_global( kg1_global(myid):kg2_global(myid) )

    ! middle points for yg (.not. equal to y in general)
    yg_m = 0.5d0*( yg(2:nyg) + yg(1:nyg-1) )

    ! middle points for yg_m (.not. equal to ym in general)
    yg_mm = 0.5d0*( yg_m(2:nyg-1) + yg_m(1:nyg-2) )

    ! local minimum grid size for CFL
    dxmin = Minval ( xg_global(2:nxg_global) - xg_global(1:nxg_global-1) )
    dymin = Minval ( yg_global(2:nyg_global) - yg_global(1:nyg_global-1) )
    dzmin = Minval ( zg_global(2:nzg_global) - zg_global(1:nzg_global-1) )

    ! total domain size
    Lx = x_global(nx_global) - x_global(1)
    Ly = y_global(ny_global) - y_global(1)
    Lz = z_global(nz_global) - z_global(1)



    !---------------------------Fourier transform---------------------------!
    If ( myid==0 ) Write(*,*) 'initializing FFT...'
    ! initialize MPI FFTW
    Call fftw_mpi_init()

    ! Fourier constant grid spacing
    dx = dxmin
    dz = dzmin

    ! length for periodic domain
    Lxp = Lx - dx 
    Lzp = Lz - dz

    ! global points for periodic domain in physical space
    nxp_global = nxm_global - 1
    nzp_global = nzm_global - 1

    ! global indices for fourier modes starting from 0
    mx_global = nxp_global - 1
    mz_global = nzp_global - 1

    ! Get local sizes:
    ! local data size in x direction
    nxp = nxp_global
    mx  =  mx_global
    ! local data size in z direction (note dimension reversal)
    alloc_local = fftw_mpi_local_size_2d(nzp_global, nxp_global, comm, nzp, local_k_offset)
    mz  = nzp - 1

    ! sanity check and restrictions in fftw
    If ( (nzp/=nzm .And. myid/=nprocs-1) .Or. (nzp/=nzm-1 .And. myid==nprocs-1) ) Then 
       Write(*,*) nzp,nzm
       Stop 'Error: something wrong in FFTW size'
    End If

    ! allocate variables
    cplane_fft = fftw_alloc_complex(alloc_local)
    Call c_f_pointer(cplane_fft,plane,[nxp,nzp])
    plane_hat(0:,0:) => plane
    Allocate ( rhs_p_hat ( 0:mx, 2:nyg-1, 0:mz ) )
    Allocate ( rhs_aux   ( 2:nyg-1 ) )
   
    ! create MPI plan for forward DFT (note dimension reversal and transposed_out/in)
    plan_d = fftw_mpi_plan_dft_2d( nzp_global, nxp_global, plane, plane_hat,           & 
             comm,  FFTW_FORWARD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT) ) 
    plan_i = fftw_mpi_plan_dft_2d( nzp_global, nxp_global, plane_hat, plane,           & 
             comm, FFTW_BACKWARD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN)  ) 

    ! global Fourier coeficients with modified wave-number for the second derivative
    Allocate ( kxx(0:mx_global), kzz(0:mz_global) ) 
    kxx = 0d0
    kzz = 0d0
    Do i = 0, Ceiling( Real(nxp_global)/2d0 )
       kxx(i) = 2d0*( dcos(2d0*pi*Real(i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0  
    End do
    Do i = Ceiling( Real(nxp_global)/2d0 )+1, mx_global
       kxx(i) = 2d0*( dcos(2d0*pi*Real(-nxp_global+i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0
    End do

    Do k = 0, Ceiling( Real(nzp_global)/2d0 )
       kzz(k) = 2d0*( dcos(2d0*pi*Real(k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0  
    End do
    Do k = Ceiling( Real(nzp_global)/2d0 )+1, mz_global
       kzz(k) = 2d0*( dcos(2d0*pi*Real(-nzp_global+k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
    End do
    
    ! MPI mapping for z-modes: from local to global without transposed_out/in
    ! Not used in this version
    Allocate ( kmode_map(0:mz) ) 
    kmode_map = 0
    Do k = 0, mz
       kmode_map(k) = k + myid*nslices_z
    End Do

    ! FFTW+MPI mapping for x and z-modes when using FFTW with transposed_out/in
    ! from local to global
    ! this needs (mz_global+1)*(mx_global+1)/nprocs to be an integer 
    If ( Mod((mz_global+1)*(mx_global+1),nprocs)/=0 ) Stop 'Error: (mz_global+1)*(mx_global+1)/nprocs should be an integer'
    Allocate ( imode_map_fft(0:mx_global,0:mz) ) 
    Allocate ( kmode_map_fft(0:mx_global,0:mz) ) 
    Do i = 0, mx_global
       Do k = 0, mz            
          pos = i + (mx_global+1)*k + (mz_global+1)*(mx_global+1)/nprocs*myid
          imode_map_fft(i,k) = Floor( Real(pos/(mz_global+1)) )
          kmode_map_fft(i,k) = Mod  ( pos, mz_global+1 )
          ! sanity check
       end Do
    End Do

    ! Sanity check for FFTW mapping
    Allocate(A_kmodes      (0:mx_global,0:mz_global))
    Allocate(A_kmodes_local(0:mx_global,0:mz_global))
    A_kmodes       = 0
    A_kmodes_local = 0
    Do i = 0, mx_global
       Do k = 0, mz            
          A_kmodes_local( imode_map_fft(i,k), kmode_map_fft(i,k) ) =  A_kmodes_local(imode_map_fft(i,k), kmode_map_fft(i,k) ) + 1
       end Do
    End Do
    Call MPI_AllReduce(A_kmodes_local,A_kmodes,(mx_global+1)*(mz_global+1),MPI_integer,MPI_sum,comm,ierr)
    If ( Any(A_kmodes>1) .Or. Any(A_kmodes==0) ) Stop 'Error: wrong combination of nx, nz and processors'
    Deallocate(A_kmodes)
    Deallocate(A_kmodes_local)
     
    !------------------------Tridiagonal linear solver-------------------------!
    If ( myid==0 ) Write(*,*) 'initializing pressure solver...'
    Allocate ( pivot(nyg) )
    Allocate ( Dyy(2:nyg-1,2:nyg-1), M(2:nyg-1,2:nyg-1) )
    Allocate ( D(2:nyg-1), DL(2:nyg-2), DU(2:nyg-2) )
     
    ! second derivative matrix for pressure (full data in y assumed)
    Dyy = 0d0
    Do j=3,nyg-2

       a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
       b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
       c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 

       Dyy(j,j+1) = a
       Dyy(j,j-1) = c 
       Dyy(j,j  ) = b

    End Do

    ! Boundary conditions for pressure (full data in y assumed)
    j = 2
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    ! Dirichlet in V: p(1)==p(2) 
    Dyy(2,2)   = b + c 
    Dyy(2,3)   = a
    coef_bc_1  = c

    j = nyg-1
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    ! Dirichlet in V: p(nyg)==p(nyg-1) 
    Dyy(nyg-1,nyg-1) = a + b
    Dyy(nyg-1,nyg-2) = c
    coef_bc_2        = a

    Allocate ( bc_1(2:nxg-1,2:nzg-1), bc_2(2:nxg-1,2:nzg-1) )
    Allocate ( bc_1_hat(0:mx,0:mz),   bc_2_hat(0:mx,0:mz)   )

    ! some parameters for linear solver
    nr   = nym
    nrhs = 1

    !--------------------interpolation weights--------------------!
    in1 = 1
    in2 = 2 
    If ( in2==1 ) Write(*,*) 'Conservative interpolations'
    Allocate ( weight_y_0(ny), weight_y_1(ny) )
    weight_y_0 = ( yg(2:nyg) - y(1:ny) ) / ( yg(2:nyg) - yg(1:nyg-1)  )
    weight_y_1 = 1d0 - weight_y_0

    !------------------------statistics---------------------------!
    Allocate (  Umean(nyg), Vmean(nyg),  Wmean(nyg)              )
    Allocate ( U2mean(nyg), V2mean(nyg), W2mean(nyg), UVmean(nyg))
    Umean  = 0d0
    Vmean  = 0d0
    Wmean  = 0d0
    U2mean = 0d0
    V2mean = 0d0
    W2mean = 0d0
    UVmean = 0d0

    !------------------------Runge-Kutta 3-------------------------!
    If ( myid==0 ) Write(*,*) 'initializing time integration...'
    Allocate( rk_coef(3,3), rk_t(3) )

    rk_t(1)      =  8d0/15d0
    rk_t(2)      =  2d0/3d0
    rk_t(3)      =  1d0

    rk_coef      =  0d0
    rk_coef(1,1) =  8d0/15d0
    rk_coef(2,1) =  1d0/4d0
    rk_coef(2,2) =  5d0/12d0
    rk_coef(3,1) =  1d0/4d0
    rk_coef(3,2) =  0d0
    rk_coef(3,3) =  3d0/4d0

    Allocate ( Fu1 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu2 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu3 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )

    Allocate ( Fv1 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv2 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv3 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )

    Allocate ( Fw1 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw2 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw3 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )

    !-------------------compute initial mass flow-----------------!
    Call compute_mean_mass_flow_U(U,Qflow_x_0)

    !----------------------sgs model------------------------------!
    If ( myid==0 ) Write(*,*) 'initializing SGS model...'

    ! Interior points only
    If (LES_model>0) Then
       Allocate( Lij    (2:nxg  ,2:nyg-1,2:nzm+2 ,6) ) ! These need one more points in each periodic direction (for filtering)
       Allocate( Mij    (2:nxg-1,2:nyg-1,2:nzm+1, 6) )
       Allocate( Sij    (2:nxg  ,2:nyg  ,2:nzm+2 ,6) ) ! These need one more points in each direction (for filtering)
       Allocate( Rij    (2:nxg  ,2:nyg  ,2:nzm+2 ,6) ) ! These need one more points in each direction (for filtering)
       Allocate( S      (2:nxg-1,2:nyg-1,2:nzm+1)    )
       
       ! Tensor buffer
       Allocate( ten_buf  (1:nxg  ,1:nyg  ,1:nzm+2 ,6) )
       Allocate( ten_buf_2(1:nxg  ,1:nyg  ,1:nzm+2 ,9) )
       
       ! Filtered velocities
       Allocate( Uf (1:nx ,1:nyg,1:nzm+2) )
       Allocate( Vf (1:nxg,1:ny ,1:nzm+2) )
       Allocate( Wf (1:nxg,1:nyg,1:nz ) )
       
    End If

    ! Eddy viscosity
    Allocate( nu_t     (1:nxg,1:nyg,1:nzg) )
    Allocate( nu_t_hat (1:nxg,1:nyg,1:nzg) )

    ! filter size
    ! (moved to filter)

    !--------------------Wall-models---------------------!

    ! ui = alpha_i dui/dy
    Allocate( alpha_x(1:nx ,1:2,1:nzg) )
    Allocate( alpha_z(1:nxg,1:2,1:nz ) )    

    alpha_x = 0d0
    alpha_z = 0d0

    !--------------------SMARTIES------------------------!

    nagents = NUM_AGENTS/nprocs/2
    nx_nagents = (nx-2)/nagents
    nz_nagents = (nz-2) 
    Allocate( actions(2,nagents,NUM_ACTIONS) )
    Allocate( rewards(2,nagents), rewards_old(2,nagents) )
    Allocate( states(2,nagents,STATE_SIZE) )

    !-------------------------Done--------------------------------!
    Call Mpi_barrier(comm,ierr)

    ! Measure time
    time1 = MPI_WTIME()
    
  End Subroutine initialize

  
End Module initialization
