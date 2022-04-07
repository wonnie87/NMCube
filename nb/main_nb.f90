program main_nb
!! Purpose: To find the dynamic response of a metabeam of bistable elements
!!         using Newmark-beta method
!!
!! Author: Written by Myungwon Hwang (hwang125@purdue.edu)
!!
!! Record of revisions:
!!    DATE        Programmer                      Decsription
!! ==========  ================  =================================================
!! 11/03/2020   Myungwon Hwang   Rev00: Initial working program
!! 08/10/2021   Myungwon Hwang   Syntax update to match that of RK solver
!! 11/30/2021   Myungwon Hwang   Rev01: Output data structure match that of RK solver
!! 03/23/2022   Myungwon Hwang   Rev02: Refactored the code
!!

use hdf5_helpers
use mpi
use sub_nb
        
implicit none
!include 'mpif.h' ! if 'use mpi' does not work

!! Data dictionary: Constants
!! Data dictionary: MPI-related
integer :: mpi_ierr, noProc, procID, status(MPI_STATUS_SIZE), mpi_errcode
integer, allocatable, dimension(:) :: scatter_sc, scatter_sc2, scatter_sc_m ! No. of data to be sent for each process for MPI_SCATTERV
integer, allocatable, dimension(:) :: scatter_disp, scatter_disp2, scatter_disp_m ! The starting counter for MPI_SCATTERV
real :: tWallclockStart, tWallclockEnd ! wall clock time of the start and the end of the execution
!! Data dictionary: hdf5-related
character(88) :: filename
integer(HID_T) :: file_id, dspace_m_id, dspace_t_id, dspace_u_id, dspace_f_id, mspace_t_id, mspace_u_id, mspace_f_id
integer(HID_T) :: aspace_PT_id, aspace_N_id, aspace_s_id, aspace_f_UC_id, attr_id
integer(HID_T) :: dset_m_id, dset_t_id, dset_u_id, dset_udot_id, dset_f_id, dset_uddot_id 
integer(HSIZE_T), dimension(1:1) :: dimsa_PT, dimsa_N, dimsa_s, dimsa_f_UC
integer(HSIZE_T), dimension(1:1) :: dims_m, dims_t, dimsm_t
integer(HSIZE_T), dimension(1:1) :: hslab_t_offset, hslab_t_count, hslab_t_stride, hslab_t_block
integer(HSIZE_T), dimension(1:2) :: dims_u, dimsm_u, dims_f, dimsm_f
integer(HSIZE_T), dimension(1:2) :: hslab_offset, hslab_count, hslab_stride, hslab_block
integer(HSIZE_T), dimension(1:2) :: hslab_f_offset, hslab_f_count, hslab_f_stride, hslab_f_block
!integer(HID_T) :: plist_id
integer :: hdf5_ierr
!integer(hsize_t), dimension(7) :: offset
character(len=80) :: dsetname
!! Data dictionary: input parameters
integer :: N_global ! Total No. of unit cells
integer :: prob_flag ! 1: pendula chain, 2: phi-4 lattice, 4: metabeam
real, allocatable, dimension(:) :: m_global ! global mass array
real, allocatable, dimension(:) :: b_global ! global damping array
integer, dimension(0:2) :: BC ! boundary condition flag
real, allocatable, dimension(:) :: s ! parameter array
integer, allocatable, dimension(:) :: LC
real, allocatable, dimension(:) :: LC_val
logical, allocatable, dimension(:) :: LC_mask
integer :: LC_dim2
integer :: LC_loc_dim2
integer, allocatable, dimension(:,:) :: LC_loc
real, allocatable, dimension(:,:) :: LC_val_loc
integer :: LC_flag ! Load case (1: single-point horizontal, 2: distributed vertical)
integer :: startnode ! global starting node of the distributed vertical load
integer :: endnode ! global end node of the distributed vertical load
real :: forceIn ! magnitude of input force
real :: freqIn ! input frequency
real :: vFlow ! flow speed
real :: t_UC ! time for the flow to travel one unit length; Or, time phase
integer :: N_NR ! max number of Newton-Raphson iterations
real :: tol_NR ! convergence criteria for NR method
integer :: N_inv ! max number of iterations for conjugate gradient method
real :: tol_inv ! tolerance for conjugate gradient method
real :: gamma_nb ! for linear acceleration method
real :: beta_nb ! for linear acceleration method
integer, parameter :: DBL=SELECTED_REAL_KIND(p=10,r=100)
real(kind=DBL) :: dt ! numerical time steps
real(kind=DBL) :: dtWrite ! output write frequency
real :: tStart, tEnd ! simulation start and end times
integer :: Nt ! No. of time steps
integer :: dNt_output ! Output results at every `dNt_output' steps
!!integer(kind=8) :: Nt ! replace the above if Nt is more than 2,147,483,647
real, allocatable, dimension(:) :: u ! (initial) displacement condition of all unit cells
real, allocatable, dimension(:) :: udot ! (initial) velocity condition of all unit cells
!! Data dictionary: global variables
real, allocatable, dimension(:) :: uddot ! derived acceleration of all unit cells
real, allocatable, dimension(:) :: a1_global ! a1=4*m/(dt**2)+2*c/(dt)
real, allocatable, dimension(:) :: a2_global ! a2=4*m/(dt)+c
real, allocatable, dimension(:) :: a3_global ! for linear acceleration method
real :: resSum ! Sum of residual elements for NR convergence
integer :: cnt_NR, cnt_inv ! number of NR and matrix inversion iterations for each time iteration
real :: t_NR_total = 0. ! cumulative time spent on NR iterations
real :: t_inv_total = 0. ! cumulative time spent on matrix inversion iterations
!! Data dictionary: local variables
integer :: N_loc ! No. of unit cells per process
integer :: N_fl, N_ceil ! Min and max number of N_loc across the processes
real, allocatable, dimension(:) :: u_loc ! displacements of unit cells
real, allocatable, dimension(:) :: udot_loc ! velocities of unit cells for each process (if needed)
real, allocatable, dimension(:) :: uddot_loc ! accelerations of unit cells for each process (if needed)
integer :: startnode_loc, endnode_loc ! corresponding local node numbers
real, allocatable, dimension(:) :: p_loc ! external force
real, allocatable, dimension(:) :: phat_loc ! phat(i+1)=p(i+1)+a1*u(i)+a2*udot(i)+m*uddot(i)
real, allocatable, dimension(:) :: a1_loc, a2_loc, a3_loc ! corresponding local array
!! Data dictionary: temporary computational variables
integer :: it, it2, ind
integer :: cnt1
real(kind=DBL) :: t ! current simulation time
real, allocatable, dimension(:) :: fs
real, allocatable, dimension(:) :: uOld_loc ! displacements of the previous time step
integer :: tmpI1
real :: tmpR1, tmpR2
character(len=20) :: tmpStr1
real, allocatable, dimension(:) :: f_val, f_val_loc ! store force value at each iteration
!! Data dictionary: file I/O
integer :: f_ierr, f_stat
character(len=80) :: f_msg

call MPI_INIT(mpi_ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpi_ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, procID, mpi_ierr)

if (procID == 0) then
    call h5open_f(hdf5_ierr)
    open(unit=1, file='../numMethod.inp', status='old', action='read', iostat=f_ierr, iomsg=f_msg)
    
    if (f_ierr == 0) then
        write (*,'(A)') " >> numMethod.inp is opened successfully."
        open(unit=2, file='../design.inp', status='old', action='read', iostat=f_ierr, iomsg=f_msg)
        if (f_ierr == 0) then
            write (*,'(A)') " >> design.inp is opened successfully."
            write (*,'(A)') " >>"

            120 format (T3, ' Simulation time ', T22, 'NR iterations', T37 'Inv iterations')
            121 format (T3, '=================', T22, '=============', T37 '==============')
!            122 format (T3, I8, T18, I8, T33, I8)
            122 format (T3, F14.8, T22, I8, T37, I8)
            123 format (T4, "The total execution time is ", F14.8, " seconds.")
            124 format (T4, "The total time for inverse calculation: ", F14.8, " seconds.")
            125 format (T4, "The total time for NR iterations: ", F14.8, " seconds.")
            open(unit=12, file='main_nb.stat', status='replace', action='write', iostat=f_ierr, iomsg=f_msg)
            write (12,120)
            write (12,121)

            read (1,*); read (1,*)
            read (1,*); read (1,*) gamma_nb
            read (1,*); read (1,*) beta_nb
            read (1,*); read (1,*) tol_NR
            read (1,*); read (1,*) tol_inv
            read (1,*); read (1,*) N_NR
            read (1,*); read (1,*) N_inv
            read (1,*); read (1,*) dt
            read (1,*); read (1,*) tStart
            read (1,*); read (1,*) tEnd
            read (1,*); read (1,*) dtWrite

            read (2,'(/A)') filename
            call h5fcreate_f("../outputs/"//trim(adjustl(filename))//"_NB.h5", H5F_ACC_TRUNC_F, file_id, hdf5_ierr)
            read (2,'(/I10)') N_global
            if (N_global < noProc) then
                write (*,'(A)') " >> The number of unit cells is less than the number of processes. &
                & Use NP value less than the number of unit cells."
                call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
            end if
            read (2,'(/I10)') prob_flag
            read (2,*); read (2,*) BC(0:2)

            if (prob_flag == 1) then
                s_dim = 4
            else if (prob_flag == 2) then
                s_dim = 5
            else if (prob_flag == 4) then
                s_dim = 12
            end if
            
            allocate(s(s_dim), STAT=f_stat)
            read (2,*); read (2,*) s(1:s_dim)
            read (2,*); read (2,*) LC_dim2

            allocate(f_val(LC_dim2), STAT=f_stat)
            allocate(LC(3*LC_dim2), STAT=f_stat)
            read (2,*)
            do it = 1, LC_dim2
                read (2,*) LC(3*it-2:3*it)
            end do
            
            allocate(LC_val(7*LC_dim2), STAT=f_stat)
            read (2,*)
            do it = 1, LC_dim2
                read (2,*) LC_val(7*it-6:7*it)
            end do

            read (2,'(/I10)') DoF

            allocate(m_global(DoF*N_global), STAT=f_stat)
            read (2,*)
            do it = 1, N_global
                read (2,*) m_global(DoF*(it-1)+1:DoF*it)
            end do

            allocate(b_global(DoF*N_global), STAT=f_stat)
            read (2,*)
            do it = 1, N_global
                read (2,*) b_global(DoF*(it-1)+1:DoF*it)
            end do

            allocate(u(DoF*N_global), STAT=f_stat)
            allocate(udot(DoF*N_global), STAT=f_stat)
            allocate(uddot(DoF*N_global), STAT=f_stat)
            read (2,*)
            do it = 1, N_global
                read (2,*) u(DoF*(it-1)+1:DoF*it)
            end do

            read (2,*)
            do it = 1, N_global
                read (2,*) udot(DoF*(it-1)+1:DoF*it)
            end do

            allocate(a1_global(DoF*N_global), STAT=f_stat)
            allocate(a2_global(DoF*N_global), STAT=f_stat)
            allocate(a3_global(DoF*N_global), STAT=f_stat)

        else
            write (*,'(A)') f_msg
            call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
        end if
    else
        write (*,'(A)') f_msg
        call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
    end if
end if

call MPI_BCAST(gamma_nb, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(beta_nb, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tol_NR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tol_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(N_NR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(N_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tStart, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tEnd, 1, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(dtWrite, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

call MPI_BCAST(N_global, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
!tol_NR = N_global*tol_NR
!tol_inv = N_global*tol_inv
call MPI_BCAST(prob_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(BC(0), 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(s_dim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
if (procID /= 0) then
    allocate(s(s_dim), STAT=f_stat)
    allocate(LC(3*LC_dim2), STAT=f_stat)
    allocate(LC_val(7*LC_dim2), STAT=f_stat)
end if
call MPI_BCAST(s(1), s_dim, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(LC(1), 3*LC_dim2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(LC_val(1), 7*LC_dim2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(DoF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)

allocate(fs(DoF), STAT=f_stat)
if (prob_flag == 1) then
    p_calc_fs => calc_fs_pendula
    p_calc_kThat => calc_kThat_pendula
else if (prob_flag == 2) then
    p_calc_fs => calc_fs_phi4
    p_calc_kThat => calc_kThat_phi4
else if (prob_flag == 4) then
    p_calc_fs => calc_fs_metabeam
    p_calc_kThat => calc_kThat_metabeam
end if

N_fl = N_global/noProc
N_ceil = N_fl + MOD(N_global, noProc)
if (procID == noProc-1) then
    N_loc = N_ceil
else
    N_loc = N_fl
end if
write (*,'(A, I4, A, I8, A)') " >> proc ", procID, " has ", N_loc, " unit cells."

LC_mask = ( procID*N_fl < LC(1:3*LC_dim2:3) .and. LC(1:3*LC_dim2:3) <= procID*N_fl+N_loc )
LC_loc_dim2 = count( LC_mask )
if ( LC_loc_dim2 > 0 ) then
    allocate(LC_loc(3, LC_loc_dim2), STAT=f_stat)
    allocate(LC_val_loc(7, LC_loc_dim2), STAT=f_stat)
    allocate(f_val_loc(LC_loc_dim2), STAT=f_stat)
end if
cnt1 = 1
do it = 1, LC_dim2
    if ( LC_mask(it) ) then
        LC_loc(:, cnt1) = LC(3*it-2:3*it)
        LC_loc(1, cnt1) = LC(3*it-2) - procID*N_fl
        LC_val_loc(:, cnt1) = LC_val(7*it-6:7*it)
        cnt1 = cnt1 + 1
    end if
end do

!! Set parameters for MPI_SCATTERV
if (procID == 0) then
    allocate(scatter_sc(noProc), STAT=f_stat)
    allocate(scatter_disp(noProc), STAT=f_stat)
    allocate(scatter_sc2(noProc), STAT=f_stat)
    allocate(scatter_disp2(noProc), STAT=f_stat)
    allocate(scatter_sc_m(noProc), STAT=f_stat)
    allocate(scatter_disp_m(noProc), STAT=f_stat)
    scatter_disp(1) = 0
    scatter_disp2(1) = 0
    scatter_disp_m(1) = 0
    do it = 1, noProc
        if (it == noProc) then
            scatter_sc(it) = DoF*(N_global/noProc + MOD(N_global,noProc))
!            scatter_sc2(it) = LC_loc_dim2
            scatter_sc_m(it) = 2*(N_global/noProc + MOD(N_global,noProc))
        else
            scatter_sc(it) = DoF*(N_global/noProc)
            scatter_disp(it+1) = scatter_disp(it) + scatter_sc(it)
!            scatter_sc2(it) = LC_loc_dim2
!            scatter_disp2(it+1) = scatter_disp2(it) + scatter_sc2(it)
            scatter_sc_m(it) = 2*(N_global/noProc)
            scatter_disp_m(it+1) = scatter_disp_m(it) + scatter_sc_m(it)
        end if
        if (it == 1) then
            scatter_sc2(it) = LC_loc_dim2
            scatter_disp2(it+1) = scatter_disp2(it) + scatter_sc2(it)
        else if (it == noProc) then
            call MPI_RECV(scatter_sc2(it), 1, MPI_INTEGER, it-1, 0, MPI_COMM_WORLD, status, mpi_ierr)
        else
            call MPI_RECV(scatter_sc2(it), 1, MPI_INTEGER, it-1, 0, MPI_COMM_WORLD, status, mpi_ierr)
            scatter_disp2(it+1) = scatter_disp2(it) + scatter_sc2(it)
        end if
    end do
else
    call MPI_SEND(LC_loc_dim2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, mpi_ierr)
end if

if (procID == 0) then
! STEP 1: Initial calculations
! Parallelization of this step may improve the scaling with number of processes
    do it=1, N_global
        if (it == 1) then
            call p_calc_fs(BC(1), 0, 2*DoF-1, u(DoF*(it-1)+1), s, fs)
        else if (it == N_global) then
            call p_calc_fs(BC(2), -DoF, DoF-1, u(DoF*(it-2)+1), s, fs)
        else
            call p_calc_fs(BC(0), -DoF, 2*DoF-1, u(DoF*(it-2)+1), s, fs)
        end if

        do it2=1, DoF
            ind = DoF*it-DoF+it2
! Need to be modified if nonzero external force p
            uddot(ind) = (-b_global(ind)*udot(ind)-fs(it2))/m_global(ind)
        end do
    end do

    a1_global = m_global/(beta_nb*dt**2) + b_global*gamma_nb/(beta_nb*dt)
    a2_global = m_global/(beta_nb*dt) + b_global*(gamma_nb/beta_nb-1)
    a3_global = m_global*(1/(2*beta_nb)-1) + b_global*dt*(gamma_nb/(2*beta_nb)-1)
! end STEP 1
end if

allocate(u_loc(DoF*(N_loc+2)), STAT=f_stat)
allocate(uOld_loc(DoF*(N_loc+2)), STAT=f_stat)
allocate(udot_loc(DoF*N_loc), STAT=f_stat)
allocate(uddot_loc(DoF*N_loc), STAT=f_stat)
allocate(p_loc(DoF*N_loc), STAT=f_stat)
allocate(phat_loc(DoF*N_loc), STAT=f_stat)
allocate(a1_loc(DoF*N_loc), STAT=f_stat)
allocate(a2_loc(DoF*N_loc), STAT=f_stat)
allocate(a3_loc(DoF*N_loc), STAT=f_stat)

call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
tWallclockStart = MPI_WTIME()

! Distribute global data to each process
call MPI_SCATTERV(u, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, u_loc(DoF+1), DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_SCATTERV(udot, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, udot_loc, DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_SCATTERV(uddot, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, uddot_loc, DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_SCATTERV(a1_global, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, a1_loc, DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_SCATTERV(a2_global, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, a2_loc, DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_SCATTERV(a3_global, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, a3_loc, DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

deallocate(a1_global, STAT=f_stat)
deallocate(a2_global, STAT=f_stat)
deallocate(a3_global, STAT=f_stat)

! Update ghost cells
if (noProc /= 1) then
    if (procID == 0) then 
        call MPI_SEND(u_loc(DoF*N_loc+1), DoF, MPI_DOUBLE_PRECISION, procID+1, 1, &
        & MPI_COMM_WORLD, mpi_ierr) 
        call MPI_RECV(u_loc(DoF*N_loc+DoF+1), DoF, MPI_DOUBLE_PRECISION, procID+1, 0, &
        & MPI_COMM_WORLD, status, mpi_ierr) 
    else if (procID == noProc-1) then
        call MPI_SEND(u_loc(DoF+1), DoF, MPI_DOUBLE_PRECISION, procID-1, 0, &
        & MPI_COMM_WORLD, mpi_ierr)
        call MPI_RECV(u_loc(1), DoF, MPI_DOUBLE_PRECISION, procID-1, 1, &
        & MPI_COMM_WORLD, status, mpi_ierr)
    else
        call MPI_SENDRECV(u_loc(DoF+1), DoF, MPI_DOUBLE_PRECISION, procID-1, 0, &
        & u_loc(1), DoF, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, mpi_ierr)
        call MPI_SENDRECV(u_loc(DoF*N_loc+1), DoF, MPI_DOUBLE_PRECISION, procID+1, 1, &
        & u_loc(DoF*N_loc+DoF+1), DoF, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, mpi_ierr)
    end if
end if

t = tStart
if ( procID == 0 ) then
    cnt1 = 0

    dimsa_PT = (/ 1 /)
    dimsa_N = (/ 1 /)
    dimsa_s = (/ s_dim /)
    dimsa_f_UC = (/ LC_dim2 /)
    call h5screate_simple_f(1, dimsa_PT, aspace_PT_id, hdf5_ierr)
    call h5screate_simple_f(1, dimsa_N, aspace_N_id, hdf5_ierr)
    call h5screate_simple_f(1, dimsa_s, aspace_s_id, hdf5_ierr)
    call h5screate_simple_f(1, dimsa_f_UC, aspace_f_UC_id, hdf5_ierr)
    call h5acreate_f(file_id, "ProblemType", H5T_NATIVE_INTEGER, aspace_PT_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, prob_flag, dimsa_PT, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "N", H5T_NATIVE_INTEGER, aspace_N_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_global, dimsa_N, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "Parameters", H5T_NATIVE_REAL, aspace_s_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, s, dimsa_s, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "f_UC", H5T_NATIVE_INTEGER, aspace_f_UC_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, LC(1:3*LC_dim2:3), dimsa_f_UC, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "f_DoF", H5T_NATIVE_INTEGER, aspace_f_UC_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, LC(2:3*LC_dim2:3), dimsa_f_UC, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5sclose_f(aspace_N_id, hdf5_ierr)
    call h5sclose_f(aspace_s_id, hdf5_ierr)
    call h5sclose_f(aspace_f_UC_id, hdf5_ierr)

    dims_m = (/ DoF*N_global /)
    call h5screate_simple_f(1, dims_m, dspace_m_id, hdf5_ierr)
    call h5dcreate_f(file_id, "/m", H5T_NATIVE_REAL, dspace_m_id, dset_m_id, hdf5_ierr)
    call h5dwrite_f(dset_m_id, H5T_NATIVE_DOUBLE, m_global, dims_m, hdf5_ierr)
    call h5dclose_f(dset_m_id, hdf5_ierr)
    call h5sclose_f(dspace_m_id, hdf5_ierr)

    dims_t = (/ NINT((tEnd-tStart)/dtWrite)+1 /)
    dims_u = (/ DoF*N_global, NINT((tEnd-tStart)/dtWrite)+1 /)
    dims_f = (/ NINT((tEnd-tStart)/dtWrite)+1, LC_dim2 /)
    call h5screate_simple_f(1, dims_t, dspace_t_id, hdf5_ierr)
    call h5screate_simple_f(2, dims_u, dspace_u_id, hdf5_ierr)
    call h5screate_simple_f(2, dims_f, dspace_f_id, hdf5_ierr)
    call h5dcreate_f(file_id, "/t", H5T_NATIVE_REAL, dspace_t_id, dset_t_id, hdf5_ierr)
    call h5dcreate_f(file_id, "/u", H5T_NATIVE_REAL, dspace_u_id, dset_u_id, hdf5_ierr)
    call h5dcreate_f(file_id, "/udot", H5T_NATIVE_REAL, dspace_u_id, dset_udot_id, hdf5_ierr)
    call h5dcreate_f(file_id, "/uddot", H5T_NATIVE_REAL, dspace_u_id, dset_uddot_id, hdf5_ierr)
    call h5dcreate_f(file_id, "/f", H5T_NATIVE_REAL, dspace_f_id, dset_f_id, hdf5_ierr)

    hslab_t_offset = (/ cnt1 /)
    hslab_t_count = (/ 1 /)
    hslab_t_stride = (/ 1 /)
    hslab_t_block = (/ 1 /)
    hslab_offset = (/ 0, cnt1 /)
    hslab_count = (/ DoF*N_global, 1 /)
    hslab_stride = (/ 1, 1 /)
    hslab_block = (/ 1, 1 /)
    hslab_f_offset = (/ cnt1, 0 /)
    hslab_f_count = (/ 1, LC_dim2 /)
    hslab_f_stride = (/ 1, 1 /)
    hslab_f_block = (/ 1, 1 /)
    call h5sselect_hyperslab_f(dspace_t_id, H5S_SELECT_SET_F, hslab_t_offset, hslab_t_count, &
    & hdf5_ierr, hslab_t_stride, hslab_t_block)
    call h5sselect_hyperslab_f(dspace_u_id, H5S_SELECT_SET_F, hslab_offset, hslab_count, hdf5_ierr, hslab_stride, hslab_block)
    call h5sselect_hyperslab_f(dspace_f_id, H5S_SELECT_SET_F, hslab_f_offset, hslab_f_count, &
    & hdf5_ierr, hslab_f_stride, hslab_f_block)
    dimsm_t = (/ 1 /)
    dimsm_u = (/ DoF*N_global, 1 /)
    dimsm_f = (/ 1, LC_dim2 /)
    call h5screate_simple_f(1, dimsm_t, mspace_t_id, hdf5_ierr)
    call h5screate_simple_f(2, dimsm_u, mspace_u_id, hdf5_ierr)
    call h5screate_simple_f(2, dimsm_f, mspace_f_id, hdf5_ierr)
    dims_t = (/ 1 /)
    dims_u = (/ DoF*N_global, 1 /)
    dims_f = (/ 1, LC_dim2 /)
    call h5dwrite_f(dset_t_id, H5T_NATIVE_DOUBLE, t, dims_t, hdf5_ierr, mspace_t_id, dspace_t_id)
    call h5dwrite_f(dset_u_id, H5T_NATIVE_DOUBLE, u(1:DoF*N_global), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)
    call h5dwrite_f(dset_udot_id, H5T_NATIVE_DOUBLE, udot(1:DoF*N_global), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)

    deallocate(m_global, STAT=f_stat)
    deallocate(b_global, STAT=f_stat)
end if

! Iterations for each time step dt
do
    t = t + dt
    p_loc = 0.
    phat_loc = 0.
    do it = 1, LC_loc_dim2
        call calc_load(LC_loc(1,it), LC_val_loc(1,it), t, &
            & p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) ), &
            & u_loc( DoF*LC_loc(1,it)+LC_loc(2,it) ) ) ! modify s.t. it can handle LC=13
        f_val_loc(it) = p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) )
    end do 

    do it = 1, DoF*N_loc
        phat_loc(it) = p_loc(it) + a1_loc(it)*u_loc(it+DoF) + a2_loc(it)*udot_loc(it) + &
        & a3_loc(it)*uddot_loc(it)
    end do
    uOld_loc = u_loc

    !! Check the following implementation in sub_nb.f90
    call NR_iterations(N_NR, tol_NR, N_loc, procID, noProc, prob_flag, BC, s, &
    & a1_loc, phat_loc, u_loc, N_inv, tol_inv, t, cnt_NR, cnt_inv, t_NR_total, t_inv_total)

    ind = DoF*(N_loc+1)
    uddot_loc = (u_loc(DoF+1:ind) - uOld_loc(DoF+1:ind))/(beta_nb*dt**2) - udot_loc/(beta_nb*dt) - uddot_loc*(1/(2*beta_nb)-1)
    udot_loc = (u_loc(DoF+1:ind) - uOld_loc(DoF+1:ind))*gamma_nb/(beta_nb*dt) + udot_loc*(1-gamma_nb/beta_nb) &
    & + uddot_loc*dt*(1-gamma_nb/(2*beta_nb))

    if (procID == 0) then
        write (12,122) t, cnt_NR, cnt_inv
    end if

    if ( MOD(t,dtWrite) < 0.1*dt .or. dtWrite-MOD(t,dtWrite) < 0.1*dt ) then
        call MPI_GATHERV(u_loc(DoF+1), DoF*N_loc, MPI_DOUBLE_PRECISION, u, scatter_sc, scatter_disp, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_GATHERV(uddot_loc(1), DoF*N_loc, MPI_DOUBLE_PRECISION, uddot, scatter_sc, scatter_disp, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_GATHERV(f_val_loc(1), LC_loc_dim2, MPI_DOUBLE_PRECISION, f_val, scatter_sc2, scatter_disp2, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

        if ( procID == 0 ) then
            cnt1 = cnt1 + 1
            hslab_t_offset = (/ cnt1 /)
            hslab_offset = (/ 0, cnt1 /)
            hslab_f_offset = (/ cnt1, 0 /)
            call h5sselect_hyperslab_f(dspace_t_id, H5S_SELECT_SET_F, hslab_t_offset, hslab_t_count, &
            & hdf5_ierr, hslab_t_stride, hslab_t_block)
            call h5sselect_hyperslab_f(dspace_u_id, H5S_SELECT_SET_F, hslab_offset, hslab_count, hdf5_ierr, &
            & hslab_stride, hslab_block)
            call h5sselect_hyperslab_f(dspace_f_id, H5S_SELECT_SET_F, hslab_f_offset, hslab_f_count, &
            & hdf5_ierr, hslab_f_stride, hslab_f_block)
            call h5screate_simple_f(1, dimsm_t, mspace_t_id, hdf5_ierr)
            call h5screate_simple_f(2, dimsm_u, mspace_u_id, hdf5_ierr)
            call h5screate_simple_f(2, dimsm_f, mspace_f_id, hdf5_ierr)
            call h5dwrite_f(dset_t_id, H5T_NATIVE_DOUBLE, t, dims_t, hdf5_ierr, mspace_t_id, dspace_t_id)
            call h5dwrite_f(dset_u_id, H5T_NATIVE_DOUBLE, u(1:DoF*N_global), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)
            call h5dwrite_f(dset_udot_id, H5T_NATIVE_DOUBLE, udot(1:DoF*N_global), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)
            call h5dwrite_f(dset_uddot_id, H5T_NATIVE_DOUBLE, uddot(1:DoF*N_global), dims_u, hdf5_ierr, &
            & mspace_u_id, dspace_u_id)
            call h5dwrite_f(dset_f_id, H5T_NATIVE_DOUBLE, f_val, dims_f, hdf5_ierr, mspace_f_id, dspace_f_id)
        end if
    end if

    if ( ABS(t-tEnd) < 0.1*dt ) EXIT

end do

call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
tWallclockEnd = MPI_WTIME()
if (procID == 0) then
    write (*, '(A)') " >>"
    write (*, '(A)') " >> Simulation has completed."
    write (*,'(A,F15.6,A)') " >> Wallclock time: ", tWallclockEnd-tWallclockStart, " s."
    write (*, '(A/)') " >> Results are saved as '"//trim(adjustl(filename))//"_NB.h5' in ./outputs folder."
    write (12,121)
    write (12,123) tWallclockEnd - tWallclockStart
    write (12,124) t_inv_total
    write (12,125) t_NR_total - t_inv_total
    call h5dclose_f(dset_t_id, hdf5_ierr)
    call h5dclose_f(dset_u_id, hdf5_ierr)
    call h5dclose_f(dset_udot_id, hdf5_ierr)
    call h5dclose_f(dset_uddot_id, hdf5_ierr)
    call h5dclose_f(dset_f_id, hdf5_ierr)
    call h5sclose_f(dspace_t_id, hdf5_ierr)
    call h5sclose_f(dspace_u_id, hdf5_ierr)
    call h5sclose_f(dspace_f_id, hdf5_ierr)
    call h5sclose_f(mspace_t_id, hdf5_ierr)
    call h5sclose_f(mspace_u_id, hdf5_ierr)
    call h5sclose_f(mspace_f_id, hdf5_ierr)
    call h5fclose_f(file_id, hdf5_ierr)
    call h5close_f(hdf5_ierr)
end if

close(unit=12)


call MPI_FINALIZE(mpi_ierr)

end program main_nb
