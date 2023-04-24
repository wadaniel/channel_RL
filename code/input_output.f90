!--------------------------------------!
!          Module for I/O              !
!--------------------------------------!
Module input_output

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use pressure

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------!
  !  Read input parameters from a txt file !
  !----------------------------------------!
  Subroutine read_input_parameters

    Character(200) :: dummy_line
    Real(Int64)    :: Rossby_plus, utau_

    ! processor 0 reads the data
    If ( myid==0 ) Then
       Open(1,file=filepar,access='stream',form='formatted',action='Read')

       Read(1,*) dummy_line
       Read(1,*) nx_global, ny_global, nz_global
       
       Read(1,*) dummy_line
       Read(1,*) CFL
       
       Read(1,*) dummy_line
       Read(1,*) nu   
       
       Read(1,*) dummy_line
       Read(1,*) dPdx, x_mass_cte

       Read(1,*) dummy_line
       Read(1,*) LES_model
       
       Read(1,*) dummy_line
       Read(1,*) nsteps, nsave, nstats, nmonitor
       
       Read(1,*) dummy_line
       Read(1,*) filein
       
       Read(1,*) dummy_line
       Read(1,*) fileout

       Read(1,*) dummy_line
       Read(1,*) nstep_init

       Read(1,*) dummy_line
       Read(1,*) random_init

       Read(1,*) dummy_line
       Read(1,*) Dirichlet_nu_t

       Close(1)

       utau_    = dPdx**0.5d0
       dPdx_ref = dPdx
       
    End If

    ! broadcast data to all processors
    Call Mpi_bcast ( nx_global,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast ( ny_global,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast ( nz_global,1,MPI_integer,0,comm,ierr )

    Call Mpi_bcast (      CFL,1,MPI_real8,0,comm,ierr )
    Call Mpi_bcast (       nu,1,MPI_real8,0,comm,ierr )
    Call Mpi_bcast (     dPdx,1,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( dPdx_ref,1,MPI_real8,0,comm,ierr )

    Call Mpi_bcast (   LES_model,1,MPI_integer,0,comm,ierr )

    Call Mpi_bcast (  nstep_init,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast (      nsteps,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast (       nsave,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast (      nstats,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast (    nmonitor,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast ( random_init,1,MPI_integer,0,comm,ierr )
    Call Mpi_bcast (  x_mass_cte,1,MPI_integer,0,comm,ierr )

    Call Mpi_bcast ( Dirichlet_nu_t,1,MPI_integer,0,comm,ierr )
    
  End Subroutine read_input_parameters

  !------------------------------------------------!
  !    Generates an initial condition for channel  !
  !                                                !
  ! Output: U,V,W                                  !
  !                                                !
  !------------------------------------------------!
  Subroutine init_flow
  
    Integer(Int32) :: ll, ii, jj, kk, i, j, k
    Real   (Int64) :: alpha

    If (myid==0) Then
       Write(*,*) 'Generating random initial condition'
    End If

    Do i=1,nx_global
       x_global(i) = Real(i-1,8)*0.2d0
    End Do
    x_global = 2d0*pi*x_global/x_global(nx_global-1)

    Do i=1,nz_global
       z_global(i) = Real(i-1,8)*0.1d0
    End Do
    z_global = 1d0*pi*z_global/z_global(nz_global-1)

    Do i=1,ny_global
       y_global(i) = Real(i-1,8)*0.3d0
    End Do
    y_global = 2d0*y_global/Maxval(y_global)-1d0
    
    y_global = y_global - Minval(y_global)
    y_global = y_global*2d0/Maxval(y_global)

    ! U
    Do jj=1,ny_global
       U(:,jj,:) = y_global(jj)*( 2d0-y_global(jj) )
    end Do
    Do ii=1,nx_global
       Do jj=1,nyg_global
          Do kk=1,nzg
             U(ii,jj,kk) = U(ii,jj,kk) + 0.5*(rand()-0.5)
          End Do
       End Do
    End Do

    ! V
    V = 0d0
    Do ii=1,nxg_global
       Do jj=1,ny_global
          Do kk=1,nzg
             V(ii,jj,kk) = V(ii,jj,kk) + 0.5*(rand()-0.5)
             !V(ii,jj,kk) = V(ii,jj,kk) + 0.5*dcos(y_global(jj))
          End Do
       End Do
    End Do
    ! remove mean (approx)
    !Do jj=1,ny_global
    !   V(:,jj,:) = V(:,jj,:) - Sum( V(:,jj,:) )/Real(2*nxg_global*nzg_global,8)
    !End Do

    ! W
    W = 0d0
    Do ii=1,nxg_global
       Do jj=1,ny_global
          Do kk=1,nz
             W(ii,jj,kk) = W(ii,jj,kk) + 0.5*(rand()-0.5)
             !W(ii,jj,kk) = W(ii,jj,kk) + 0.5*dsin(y_global(jj))
          End Do
       End Do
    End Do

    If ( myid==0 ) Then
       Write(*,*) 'Max U',MaxVal(U)
       Write(*,*) 'Max V',MaxVal(V)
       Write(*,*) 'Max W',MaxVal(W)
       
       Write(*,*) 'Mean U',sum(U)/Real(nx_global*nyg_global*nzg_global,8)
       Write(*,*) 'Mean V',sum(V)/Real(nxg_global*ny_global*nzg_global,8)
       Write(*,*) 'Mean W',sum(W)/Real(nxg_global*nyg_global*nz_global,8)
    End If
 
  End Subroutine init_flow

  !--------------------------------------------!
  !    Read binary snapshot: mesh, U,V and W   !
  !                                            !
  ! Input:  filein                             !
  ! Output: U,V,W,x,y,z                        !
  !                                            !
  !--------------------------------------------!
  Subroutine read_input_data2
   
    
    Character(200) :: filein2
    Integer(Int32) ::  nx_global_f,  ny_global_f,  nz_global_f, iproc, nze, nzge
    Integer(Int32) :: nxm_global_f, nym_global_f, nzm_global_f, nn(3)

    Integer(Int32) :: ii,jj,kk
    Real   (Int64) :: aa
    
    ! processor 0 Reads the all the data
    If ( myid==0 ) Then

       Call RANDOM_NUMBER(aa)
       If (aa < 1d0/6d0) Then
          filein2 = '../case_folder/channel950_AMD_RL.init'
       Elseif (aa < 2d0/6d0) Then
          filein2 = '../case_folder/channel2000_AMD_RL.init'
       Elseif (aa < 3d0/6d0) Then
          filein2 = '../case_folder/channel3000_AMD_RL.init'
       Elseif (aa < 4d0/6d0) Then
          filein2 = '../case_folder/channel4200_AMD_RL.init'
       Elseif (aa < 5d0/6d0) Then
          filein2 = '../case_folder/channel8000_AMD_RL.init'
       Else
          filein2 = '../case_folder/channel1e4_AMD_RL.init'
       End If

       Write(*,*) 'reading ',Trim(Adjustl(filein2)),'...'
       Open(1,file=filein,access='stream',form='unformatted',action='Read')       
       ! mesh
       Read(1) nx_global_f
       If ( nx_global_f/=nx_global ) Stop 'nx_f/=nx'
       Read(1) x_global
       
       Read(1) ny_global_f
       If ( ny_global_f/=ny_global ) Stop 'ny_f/=ny'
       Read(1) y_global
       
       Read(1) nz_global_f
       If ( nz_global_f/=nz_global ) Stop 'nz_f/=nz'
       Read(1) z_global
       
       Read(1) nxm_global_f
       If ( nxm_global_f/=nxm_global ) Stop 'nxm_f/=nxm'
       Read(1) xm_global
       
       Read(1) nym_global_f
       If ( nym_global_f/=nym_global ) Stop 'nym_f/=nym'
       Read(1) ym_global
       
       Read(1) nzm_global_f
       If ( nzm_global_f/=nzm_global ) Stop 'nzm_f/=nzm'
       Read(1) zm_global

    End If

    ! U
    If ( myid==0 ) Then
       ! read dummy
       Read(1) nn 
       If ( nn(1)/=nx_global .or. nn(2)/=nyg_global .or. nn(3)/=nzg_global ) Then 
          Write(*,*) 'nn',nn
          Stop 'Error! wrong size in input file (U)'
       End If
       ! read data for processor 0
       nzge = kg2_global(myid) - kg1_global(myid) + 1
       Read(1) U(:,:,1:nzg-1)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             Read(1) Uo(:,:,2:nzge-1)
             Call Mpi_send(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          Else ! especial case: U has different size for last processor
             Read(1) Uoo(:,:,2:nzge)
             Call Mpi_send(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(U,nx*nyg*nzg,Mpi_real8,0,myid,comm,istat,ierr)
    Endif

    ! V
    If ( myid==0 ) Then
       ! go to correct position. I dont know, if I dont do this it gets lost sometimes
       ! read dummy
       Read(1) nn
       If ( nn(1)/=nxg_global .or. nn(2)/=ny_global .or. nn(3)/=nzg_global ) Then 
          Write(*,*) 'nn',nn
          Stop 'Error! wrong size in input file (V)'
       End If
       ! read data for processor 0
       nzge = kg2_global(myid) - kg1_global(myid) + 1
       Read(1) V(:,:,1:nzg-1)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             Read(1) Vo(:,:,2:nzge-1) 
             Call Mpi_send(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          Else ! especial case: V has different size for last processor
             Read(1) Voo(:,:,2:nzge) 
             Call Mpi_send(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(V,nxg*ny*nzg,Mpi_real8,0,myid,comm,istat,ierr)
    Endif

    ! W
    If ( myid==0 ) Then
       ! go to correct position. I dont know, if I dont do this it gets lost sometimes
       ! read dummy
       Read(1) nn
       If ( nn(1)/=nxg_global .or. nn(2)/=nyg_global .or. nn(3)/=nz_global ) Then 
          Write(*,*) 'nn',nn
          Stop 'Error! wrong size in input file (W)'
       End If
       ! read data for processor 0
       nzge = k2_global(myid) - k1_global(myid) + 1
       Read(1) W(:,:,1:nz-1)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             Read(1) Wo(:,:,2:nze-1)
             Call Mpi_send(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,comm,ierr)
          Else ! especial case: W has different size for last processor
             Read(1) Woo(:,:,2:nze)
             Call Mpi_send(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,comm,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(W,nxg*nyg*nz,Mpi_real8,0,myid,comm,istat,ierr)
    Endif

    ! close file
    If (myid==0) Then
       Close(1)
    End If

    ! send data to all other processors
    ! mesh
    Call Mpi_bcast ( x_global,nx_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( y_global,ny_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( z_global,nz_global,MPI_real8,0,comm,ierr )

    Call Mpi_bcast ( xm_global,nxm_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( ym_global,nym_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( zm_global,nzm_global,MPI_real8,0,comm,ierr ) 

    ! interior region
    Call update_ghost_interior_planes(U,1)
    Call update_ghost_interior_planes(V,2)
    Call update_ghost_interior_planes(W,3)

    ! apply periodicity in x
    Call apply_periodic_bc_x(U,1)
    Call apply_periodic_bc_x(V,2)
    Call apply_periodic_bc_x(W,2)

    ! apply periodicity in z
    Call apply_periodic_bc_z(U,1)
    Call apply_periodic_bc_z(V,2)
    Call apply_periodic_bc_z(W,3)
    


    Do ii=1,nx
       Do jj=1,nyg
          Do kk=1,nzg
             U(ii,jj,kk) = U(ii,jj,kk) + 0.1*(rand()-0.5)
          End Do
       End Do
    End Do

    Do ii=1,nxg
       Do jj=1,ny
          Do kk=1,nzg
             V(ii,jj,kk) = V(ii,jj,kk) + 0.1*(rand()-0.5)
          End Do
       End Do
    End Do

    Do ii=1,nxg
       Do jj=1,nyg
          Do kk=1,nz
             W(ii,jj,kk) = W(ii,jj,kk) + 0.1*(rand()-0.5)
          End Do
       End Do
    End Do

    ! set solution for zero step
    Uo = U
    Vo = V
    Wo = W

  End Subroutine read_input_data2

  !--------------------------------------------!
  !    Read binary snapshot: mesh, U,V and W   !
  !                                            !
  ! Input:  filein                             !
  ! Output: U,V,W,x,y,z                        !
  !                                            !
  !--------------------------------------------!
  Subroutine read_input_data
   
    Integer(Int32) ::  nx_global_f,  ny_global_f,  nz_global_f, iproc, nze, nzge
    Integer(Int32) :: nxm_global_f, nym_global_f, nzm_global_f, nn(3)

    Integer(Int32) :: ii,jj,kk
    
    ! processor 0 Reads the all the data
    If ( myid==0 ) Then

       Write(*,*) 'reading ',Trim(Adjustl(filein)),'...'
       Open(1,file=filein,access='stream',form='unformatted',action='Read')       
       ! mesh
       Read(1) nx_global_f
       If ( nx_global_f/=nx_global ) Stop 'nx_f/=nx'
       Read(1) x_global
       
       Read(1) ny_global_f
       If ( ny_global_f/=ny_global ) Stop 'ny_f/=ny'
       Read(1) y_global
       
       Read(1) nz_global_f
       If ( nz_global_f/=nz_global ) Stop 'nz_f/=nz'
       Read(1) z_global
       
       Read(1) nxm_global_f
       If ( nxm_global_f/=nxm_global ) Stop 'nxm_f/=nxm'
       Read(1) xm_global
       
       Read(1) nym_global_f
       If ( nym_global_f/=nym_global ) Stop 'nym_f/=nym'
       Read(1) ym_global
       
       Read(1) nzm_global_f
       If ( nzm_global_f/=nzm_global ) Stop 'nzm_f/=nzm'
       Read(1) zm_global

    End If

    ! U
    If ( myid==0 ) Then
       ! read dummy
       Read(1) nn 
       If ( nn(1)/=nx_global .or. nn(2)/=nyg_global .or. nn(3)/=nzg_global ) Then 
          Write(*,*) 'nn',nn
          Stop 'Error! wrong size in input file (U)'
       End If
       ! read data for processor 0
       nzge = kg2_global(myid) - kg1_global(myid) + 1
       Read(1) U(:,:,1:nzg-1)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             Read(1) Uo(:,:,2:nzge-1)
             Call Mpi_send(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          Else ! especial case: U has different size for last processor
             Read(1) Uoo(:,:,2:nzge)
             Call Mpi_send(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(U,nx*nyg*nzg,Mpi_real8,0,myid,comm,istat,ierr)
    Endif

    ! V
    If ( myid==0 ) Then
       ! go to correct position. I dont know, if I dont do this it gets lost sometimes
       ! read dummy
       Read(1) nn
       If ( nn(1)/=nxg_global .or. nn(2)/=ny_global .or. nn(3)/=nzg_global ) Then 
          Write(*,*) 'nn',nn
          Stop 'Error! wrong size in input file (V)'
       End If
       ! read data for processor 0
       nzge = kg2_global(myid) - kg1_global(myid) + 1
       Read(1) V(:,:,1:nzg-1)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             Read(1) Vo(:,:,2:nzge-1) 
             Call Mpi_send(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          Else ! especial case: V has different size for last processor
             Read(1) Voo(:,:,2:nzge) 
             Call Mpi_send(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,comm,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(V,nxg*ny*nzg,Mpi_real8,0,myid,comm,istat,ierr)
    Endif

    ! W
    If ( myid==0 ) Then
       ! go to correct position. I dont know, if I dont do this it gets lost sometimes
       ! read dummy
       Read(1) nn
       If ( nn(1)/=nxg_global .or. nn(2)/=nyg_global .or. nn(3)/=nz_global ) Then 
          Write(*,*) 'nn',nn
          Stop 'Error! wrong size in input file (W)'
       End If
       ! read data for processor 0
       nzge = k2_global(myid) - k1_global(myid) + 1
       Read(1) W(:,:,1:nz-1)
       ! data for processor n>0    
       Do iproc = 1, nprocs-1
          nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
             Read(1) Wo(:,:,2:nze-1)
             Call Mpi_send(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,comm,ierr)
          Else ! especial case: W has different size for last processor
             Read(1) Woo(:,:,2:nze)
             Call Mpi_send(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,comm,ierr)
          End If
       Enddo       
    Else
       Call Mpi_recv(W,nxg*nyg*nz,Mpi_real8,0,myid,comm,istat,ierr)
    Endif

    ! close file
    If (myid==0) Then
       Close(1)
    End If

    ! send data to all other processors
    ! mesh
    Call Mpi_bcast ( x_global,nx_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( y_global,ny_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( z_global,nz_global,MPI_real8,0,comm,ierr )

    Call Mpi_bcast ( xm_global,nxm_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( ym_global,nym_global,MPI_real8,0,comm,ierr )
    Call Mpi_bcast ( zm_global,nzm_global,MPI_real8,0,comm,ierr ) 

    ! interior region
    Call update_ghost_interior_planes(U,1)
    Call update_ghost_interior_planes(V,2)
    Call update_ghost_interior_planes(W,3)

    ! apply periodicity in x
    Call apply_periodic_bc_x(U,1)
    Call apply_periodic_bc_x(V,2)
    Call apply_periodic_bc_x(W,2)

    ! apply periodicity in z
    Call apply_periodic_bc_z(U,1)
    Call apply_periodic_bc_z(V,2)
    Call apply_periodic_bc_z(W,3)
    


    Do ii=1,nx
       Do jj=1,nyg
          Do kk=1,nzg
             U(ii,jj,kk) = U(ii,jj,kk) + 0.1*(rand()-0.5)
          End Do
       End Do
    End Do

    Do ii=1,nxg
       Do jj=1,ny
          Do kk=1,nzg
             V(ii,jj,kk) = V(ii,jj,kk) + 0.1*(rand()-0.5)
          End Do
       End Do
    End Do

    Do ii=1,nxg
       Do jj=1,nyg
          Do kk=1,nz
             W(ii,jj,kk) = W(ii,jj,kk) + 0.1*(rand()-0.5)
          End Do
       End Do
    End Do

    ! set solution for zero step
    Uo = U
    Vo = V
    Wo = W

  End Subroutine read_input_data

  !--------------------------------------------!
  !    write binary snapshot: mesh, U,V and W  !
  !                                            !
  ! Input: U,V,W,x,y,z,xm,ym,zm                !
  ! Output: fileout                            !
  !                                            !
  !--------------------------------------------!
  Subroutine output_data

    Character(200)   :: fname
    Character(8)     :: ext
    Integer  (Int32) :: iproc, nze, nzge
    
    If ( Mod(istep,nsave)==0 ) then

       ! P
       If (pressure_computed.eqv..False.) Then
          Call compute_pressure
       End If

       ! processor 0 writes the data
       If ( myid==0 ) Then
          
          Write(ext,'(I8)') istep + nstep_init
          
          fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))
          Write(*,*) 'writting ',Trim(Adjustl(fname))
          Open(1,file=fname,access='stream',form='unformatted',action='write')
          
          ! mesh
          Write(1) Shape(x_global), x_global
          Write(1) Shape(y_global), y_global
          Write(1) Shape(z_global), z_global
          
          Write(1) Shape(xm_global), xm_global
          Write(1) Shape(ym_global), ym_global
          Write(1) Shape(zm_global), zm_global          
         
       End If

       ! U
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(U,nx*nyg*nzg,Mpi_real8,0,myid,comm,ierr)
       Else
          ! write U size
          Write(1) nx_global,nyg_global,nzg_global
          ! processor 0 writes its data
          Write(1) U(:,:,1:nzg-1) 
          ! processor 0 receives and writes rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Uo(:,:,2:nzge-1)
             Else
                Call Mpi_recv(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Uoo(:,:,2:nzge)
             End If
          End Do
       Endif

       ! V
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(V,nxg*ny*nzg,Mpi_real8,0,myid,comm,ierr)
       Else
          ! write V size
          Write(1) nxg_global,ny_global,nzg_global
          ! processor 0 writes its data
          Write(1) V(:,:,1:nzg-1)
          ! processor 0 receives and write rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Vo(:,:,2:nzge-1)
             Else
                Call Mpi_recv(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Voo(:,:,2:nzge)
             End If
          End Do
       Endif

       ! W
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(W,nxg*nyg*nz,Mpi_real8,0,myid,comm,ierr)
       Else
          ! write W size
          Write(1) nxg_global,nyg_global,nz_global
          ! processor 0 writes its data
          Write(1) W(:,:,1:nz-1)
          ! processor 0 receives and writes rest data
          Do iproc = 1, nprocs-1
             nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Wo(:,:,2:nze-1)
             Else
                Call Mpi_recv(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Woo(:,:,2:nze)
             End If
          End Do
       Endif

       ! P
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(P,nxg*nyg*nzg,Mpi_real8,0,myid,comm,ierr)
       Else
          ! write P size
          Write(1) nxg_global,nyg_global,nzg_global
          ! processor 0 writes its data
          Write(1) P(:,:,1:nzg-1)
          ! processor 0 receives and write rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Po,nxg*nyg*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Po(:,:,2:nzge-1)
             Else
                Call Mpi_recv(Poo,nxg*nyg*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Poo(:,:,2:nzge)
             End If
          End Do
       Endif

       ! nu_t, WARNING: destroying P
       P = nu_t
       If ( myid/=0 ) Then
          ! data from processor n>0    
          Call Mpi_send(P,nxg*nyg*nzg,Mpi_real8,0,myid,comm,ierr)
       Else
          ! write P size
          Write(1) nxg_global,nyg_global,nzg_global
          ! processor 0 writes its data
          Write(1) P(:,:,1:nzg-1)
          ! processor 0 receives and write rest data
          Do iproc = 1, nprocs-1
             nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
             If ( iproc<nprocs-1 ) Then
                Call Mpi_recv(Po,nxg*nyg*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Po(:,:,2:nzge-1)
             Else
                Call Mpi_recv(Poo,nxg*nyg*nzge,Mpi_real8,iproc,iproc,comm,istat,ierr)
                Write(1) Poo(:,:,2:nzge)
             End If
          End Do
       Endif
          
       ! close file
       If (myid==0) Then
          Close(1)
       End If
       
    End If
       
  End Subroutine output_data

  !----------------------------------------------!
  !   Write some basic statistics in a txt file  !
  !----------------------------------------------!
  Subroutine output_statistics

    Character(200) :: fname
    Character(8)   :: ext
    Integer(Int32) :: jj

    If ( myid==0 ) Then

       Write(ext,'(I8)') istep + nstep_init
       
       fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.stats.txt'
       Write(*,*) 'writting ',Trim(Adjustl(fname))
       Open(3,file=fname,form='formatted',action='write') 
       Write(3,'(A,4F15.8,4I5)') '%',t, Retau, utau, nu, nx_global, ny_global, nz_global, istep
       Do jj=1,nyg
          Write(3,'(8F15.8)') yg(jj), Umean(jj), Vmean(jj), Wmean(jj), U2mean(jj), V2mean(jj), &
                  W2mean(jj), UVmean(jj)
       End Do
       Close(3)

    End If

  End Subroutine output_statistics

End Module input_output
