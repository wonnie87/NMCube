program main_rk
!! Purpose: To find the dynamic response of a metabeam of bistable elements
!!         using an explicit Runge-Kutta method.
!!
!! Author: Written by Myungwon Hwang (hwang125@purdue.edu)
!!
!! Record of revisions:
!!    DATE        Programmer                      Decsription
!! ==========  ================  =================================================
!! 10/29/2020   Myungwon Hwang   Rev00: Initial working program
!! 11/30/2021   Myungwon Hwang   Rev01: Echoing simulation status; Bug fixes.
!!

use hdf5_helpers
use mpi
use sub_rk
        
implicit none
!include 'mpif.h' ! if 'use mpi' does not work

!! Data dictionary: Pointers
!procedure (calc_f_metabeam), pointer :: p_calc_f => NULL()
!! Data dictionary: Constants
!! Data dictionary: MPI-related
integer :: mpi_ierr, noProc, procID, status(MPI_STATUS_SIZE), mpi_errcode
integer, allocatable, dimension(:) :: scatter_sc, scatter_sc2, scatter_sc3 ! number of data to be sent for each process for MPI_SCATTERV
integer, allocatable, dimension(:) :: scatter_disp, scatter_disp2, scatter_disp3 ! the starting counter for MPI_SCATTERV
real :: tWallclockStart, tWallclockEnd ! wall clock time of the star and the end of the execution
!! Data dictionary: HDF5-related
character(80) :: filename, filename_h5
integer(HID_T) :: file_id, dspace_m_id, dspace_t_id, dspace_u_id, dspace_f_id, mspace_t_id, mspace_u_id, mspace_f_id
integer(HID_T) :: aspace_PT_id, aspace_N_id, aspace_G_id, aspace_L_id, aspace_k_id, aspace_f_UC_id, attr_id
integer(HID_T) :: dset_m_id, dset_t_id, dset_u_id, dset_udot_id, dset_f_id, dset_uddot_id 
integer(HSIZE_T), dimension(1:1) :: dimsa_PT, dimsa_N, dimsa_G, dimsa_L, dimsa_k, dimsa_f_UC
integer(HSIZE_T), dimension(1:1) :: dims_m, dims_t, dimsm_t
integer(HSIZE_T), dimension(1:1) :: hslab_t_offset, hslab_t_count, hslab_t_stride, hslab_t_block
integer(HSIZE_T), dimension(1:2) :: dims_u, dimsm_u, dims_f, dimsm_f
integer(HSIZE_T), dimension(1:2) :: hslab_offset, hslab_count, hslab_stride, hslab_block
integer(HSIZE_T), dimension(1:2) :: hslab_f_offset, hslab_f_count, hslab_f_stride, hslab_f_block
integer :: hdf5_ierr
!! Data dictionary: file I/O
integer :: f_ierr, f_stat
character(len=80) :: f_msg
!! Data dictionary: input parameters
integer :: N_global ! Total number of unit cells
integer :: sol_flag ! num. sol. method (10: exp RK1, 20: exp RK2, 30: exp RK3, 40: exp RK4)
integer, parameter :: DBL=SELECTED_REAL_KIND(p=10,r=100)
real(kind=DBL) :: dt ! numerical time steps
real(kind=DBL) :: dtWrite ! output write frequency
real :: tStart, tEnd ! simulation start and end times
real :: c2, c3, c4, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4 ! RK coefficents
integer :: prob_flag ! 1: pendula chain, 2: phi-4 lattice, 4: metabeam
integer, dimension(0:2) :: BC ! boundary condition
real, dimension(9) :: k = 0.0 ! spring stiffness array
real, dimension(4) :: L ! length array containing unit cell geometry information
real :: G ! gravitational acceleration
real :: u_0 ! Initial displacement of each bistable element
real :: uDot_0 ! Initial velocity of each bistable element
integer, allocatable, dimension(:) :: LC
real, allocatable, dimension(:) :: LC_val
logical, allocatable, dimension(:) :: LC_mask
integer :: LC_dim2
integer :: LC_loc_dim2
integer, allocatable, dimension(:,:) :: LC_loc
real, allocatable, dimension(:,:) :: LC_val_loc
integer :: LC_flag ! load case (1: single-point horizontal, 2: distributed vertical)
integer :: DoF ! unit cell degrees of freedom
integer :: noState ! number of states per unit cell
real, allocatable, dimension(:) :: m_global ! global mass array
real, allocatable, dimension(:) :: b_global ! global damping array
real, allocatable, dimension(:) :: x ! global array of state space (including initial condition)
                                      ! (..., u1n,u1dn,v1n,v1dn,u2n,u2dn,v2n,v2dn,u3n,u3dn, ...)
real, allocatable, dimension(:) :: k1_global ! global reaction force array
!! Data dictionary: local variables
integer :: N_fl, N_ceil ! min and max number of N_loc across the processes
integer :: N_loc ! number of unit cells per process
real, allocatable, dimension(:) :: x_loc ! state space for each local process
real, allocatable, dimension(:) :: m_loc, b_loc ! local mass/damping arrays
real, allocatable, dimension(:) :: p_loc ! external force calculated at each time step

!! Data dictionary: temporary computational variables
integer :: it, it2, it3
integer :: cnt1
real(kind=DBL) :: t ! current simulation time 
real, allocatable, dimension(:) :: k1, k2, k3, k4 ! RK method
real, allocatable, dimension(:) :: xTmp ! ...
real, allocatable, dimension(:) :: f_val, f_val_loc ! store force value at each iteration
real :: ran ! random number


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
            read (1,'(/,I3)') sol_flag
            if (sol_flag == 10) then
                read (1,*); read (1,*) b1
            else if (sol_flag == 20) then
                read (1,*); read (1,*) c2
                read (1,*); read (1,*) a21
                read (1,*); read (1,*) b1, b2
            else if (sol_flag == 30) then
                read (1,*); read (1,*) c2, c3
                read (1,*); read (1,*) a21, a31, a32
                read (1,*); read (1,*) b1, b2, b3
            else if (sol_flag == 40) then
                read (1,*); read (1,*) c2, c3, c4
                read (1,*); read (1,*) a21, a31, a32, a41, a42, a43
                read (1,*); read (1,*) b1, b2, b3, b4
            end if
            read (1,*); read (1,*) dt
            read (1,*); read (1,*) tStart
            read (1,*); read (1,*) tEnd
            read (1,*); read (1,*) dtWrite

            read (2,'(/A)') filename
            if (sol_flag == 10) then
                filename_h5 = trim(adjustl(filename))//"_RK1.h5"
            else if (sol_flag == 20) then
                filename_h5 = trim(adjustl(filename))//"_RK2.h5"
            else if (sol_flag == 30) then
                filename_h5 = trim(adjustl(filename))//"_RK3.h5"
            else if (sol_flag == 40) then
                filename_h5 = trim(adjustl(filename))//"_RK4.h5"
            end if
            call h5fcreate_f("../outputs/"//trim(filename_h5), H5F_ACC_TRUNC_F, file_id, hdf5_ierr)

            read (2,'(/I10)') N_global
            if (N_global < noProc) then
                write (*,'(A)') " >> The number of unit cells is less than the number of processes. &
                & Use NP value less than the number of unit cells."
                call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
            end if
            read (2,'(/I10)') prob_flag
            read (2,*); read (2,*) BC(0:2)
            if (prob_flag == 1) then
                read (2,*); read (2,*) G
                read (2,*); read (2,*) k(1)
                read (2,*); read (2,*) L(1)
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
                noState = 2*DoF

                allocate(m_global(2*N_global), STAT=f_stat)
                read (2,*)
                do it = 1, N_global
                    read (2,*) m_global(2*it-1:2*it)
                end do

                allocate(b_global(DoF*N_global), STAT=f_stat)
                read (2,*)
                do it = 1, N_global
                    read (2,*) b_global(DoF*it)
                end do

            else if (prob_flag == 2) then
                read (2,*); read (2,*) k(1:5)
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
                noState = 2*DoF

                allocate(m_global(DoF*N_global), STAT=f_stat)
                read (2,*)
                do it = 1, N_global
                    read (2,*) m_global(DoF*it)
                end do

                allocate(b_global(DoF*N_global), STAT=f_stat)
                read (2,*)
                do it = 1, N_global
                    read (2,*) b_global(DoF*it)
                end do

            else if (prob_flag == 4) then
                if (BC(0) == 300) then
                    read (2,*); read (2,*) k(1:9)
                else
                    read (2,*); read (2,*) k(1:8)
                end if
                read (2,*); read (2,*) L(1:4)
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
                noState = 2*DoF

                allocate(m_global(DoF*N_global), STAT=f_stat)
                read (2,*)
                do it = 1, N_global
                    read (2,*) m_global(DoF*it-5:DoF*it)
                end do

                allocate(b_global(DoF*N_global), STAT=f_stat)
                read (2,*)
                do it = 1, N_global
                    read (2,*) b_global(DoF*it-5:DoF*it)
                end do
            end if

            allocate(x(noState*N_global), STAT=f_stat)
            read (2,*)
            do it = 1, N_global
                read (2,*) x(noState*(it-1)+1:noState*it:2)
            end do
            read (2,*)
            do it = 1, N_global
                read (2,*) x(noState*(it-1)+2:noState*it:2)
            end do
            allocate(k1_global(noState*N_global), STAT=f_stat)

        else
            write (*,'(A)') f_msg
            call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
        end if
    else
        write (*,'(A)') f_msg
        call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
    end if
end if

call MPI_BCAST(sol_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
if (sol_flag == 10) then
    call MPI_BCAST(b1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
else if (sol_flag == 20) then
    call MPI_BCAST(c2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a21, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
else if (sol_flag == 30) then
    call MPI_BCAST(c2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(c3, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a21, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a31, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a32, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
else if (sol_flag == 40) then
    call MPI_BCAST(c2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(c3, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(c4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a21, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a31, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a32, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a41, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a42, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(a43, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b3, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(b4, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
end if
call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tStart, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tEnd, 1, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(dtWrite, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

call MPI_BCAST(N_global, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(prob_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(BC(0), 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
if (prob_flag == 1) then
    call MPI_BCAST(G, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(k(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(L(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (procID /= 0) then
        allocate(LC(3*LC_dim2), STAT=f_stat)
        allocate(LC_val(7*LC_dim2), STAT=f_stat)
    end if
    call MPI_BCAST(LC(1), 3*LC_dim2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_val(1), 7*LC_dim2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(DoF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(noState, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
else if (prob_flag == 2) then
    call MPI_BCAST(k(1), 5, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (procID /= 0) then
        allocate(LC(3*LC_dim2), STAT=f_stat)
        allocate(LC_val(7*LC_dim2), STAT=f_stat)
    end if
    call MPI_BCAST(LC(1), 3*LC_dim2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_val(1), 7*LC_dim2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(DoF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(noState, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
else if (prob_flag == 4) then
    if (BC(0) == 300) then
        call MPI_BCAST(k(1), 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    else
        call MPI_BCAST(k(1), 8, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    end if
    call MPI_BCAST(L(1), 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (procID /= 0) then
        allocate(LC(3*LC_dim2), STAT=f_stat)
        allocate(LC_val(7*LC_dim2), STAT=f_stat)
    end if
    call MPI_BCAST(LC(1), 3*LC_dim2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_val(1), 7*LC_dim2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(DoF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(noState, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
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
    allocate(scatter_sc3(noProc), STAT=f_stat)
    allocate(scatter_disp3(noProc), STAT=f_stat)
    scatter_disp(1) = 0
    scatter_disp2(1) = 0
    scatter_disp3(1) = 0
    do it = 1, noProc
        if (it == noProc) then
            scatter_sc(it) = noState*(N_global/noProc + MOD(N_global,noProc))
            scatter_sc2(it) = DoF*(N_global/noProc + MOD(N_global,noProc))
!            scatter_sc3(it) = LC_loc_dim2
        else
            scatter_sc(it) = noState*(N_global/noProc)
            scatter_disp(it+1) = scatter_disp(it) + scatter_sc(it)
            scatter_sc2(it) = DoF*(N_global/noProc)
            scatter_disp2(it+1) = scatter_disp2(it) + scatter_sc2(it)
!            scatter_sc3(it) = LC_loc_dim2
!            scatter_disp3(it+1) = scatter_disp3(it) + scatter_sc3(it)
        end if
        if (it == 1) then
            scatter_sc3(it) = LC_loc_dim2
            scatter_disp3(it+1) = scatter_disp3(it) + scatter_sc3(it)
        else if (it == noProc) then
            call MPI_RECV(scatter_sc3(it), 1, MPI_INTEGER, it-1, 0, MPI_COMM_WORLD, status, mpi_ierr)
        else
            call MPI_RECV(scatter_sc3(it), 1, MPI_INTEGER, it-1, 0, MPI_COMM_WORLD, status, mpi_ierr)
            scatter_disp3(it+1) = scatter_disp3(it) + scatter_sc3(it)
        end if
    end do
else
    call MPI_SEND(LC_loc_dim2, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, mpi_ierr)
end if

allocate(x_loc(noState*(N_loc+2)), STAT=f_stat)
if (prob_flag == 1) then
    allocate(m_loc(2*N_loc), STAT=f_stat) 
else if (prob_flag == 2 .or. prob_flag == 4) then
    allocate(m_loc(DoF*N_loc), STAT=f_stat) 
end if
allocate(b_loc(DoF*N_loc), STAT=f_stat) 
allocate(p_loc(DoF*N_loc), STAT=f_stat) 
allocate(xTmp(noState*(N_loc+2)), STAT=f_stat)

allocate(k1(noState*(N_loc+2)), STAT=f_stat)
if (sol_flag /= 10) then
    allocate(k2(noState*(N_loc+2)), STAT=f_stat)
    if (sol_flag /= 20) then
        allocate(k3(noState*(N_loc+2)), STAT=f_stat)
        if (sol_flag /= 30) then
            allocate(k4(noState*(N_loc+2)), STAT=f_stat)
        end if
    end if
end if

call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
tWallclockStart = MPI_WTIME()

!! Distribute global data to each process
call MPI_SCATTERV(x, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, x_loc(noState+1), noState*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
if (prob_flag == 1) then
    call MPI_SCATTERV(m_global, scatter_sc, scatter_disp, MPI_DOUBLE_PRECISION, m_loc, 2*N_loc, &
    & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
else if (prob_flag == 2 .or. prob_flag == 4) then
    call MPI_SCATTERV(m_global, scatter_sc2, scatter_disp2, MPI_DOUBLE_PRECISION, m_loc, DoF*N_loc, &
    & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
end if
call MPI_SCATTERV(b_global, scatter_sc2, scatter_disp2, MPI_DOUBLE_PRECISION, b_loc, DoF*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

!! Update ghost cells
if (noProc /= 1) then
    if (procID == 0) then
        call MPI_SEND(x_loc(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
        & MPI_COMM_WORLD, mpi_ierr)
        call MPI_RECV(x_loc(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, &
        & MPI_COMM_WORLD, status, mpi_ierr)
    else if (procID == noProc-1) then
        call MPI_SEND(x_loc(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
        & MPI_COMM_WORLD, mpi_ierr)
        call MPI_RECV(x_loc(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, &
        & MPI_COMM_WORLD, status, mpi_ierr)
    else
        call MPI_SENDRECV(x_loc(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
        & x_loc(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, mpi_ierr)
        call MPI_SENDRECV(x_loc(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
        & x_loc(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, mpi_ierr)
    end if
end if

t = tStart
if (procID == 0) then
    cnt1 = 0

    dimsa_PT = (/ 1 /)
    dimsa_N = (/ 1 /)
    if (prob_flag == 1) then
        dimsa_G = (/ 1 /)
        dimsa_L = (/ 1 /)
        dimsa_k = (/ 1 /)
    else if (prob_flag == 2) then
        dimsa_k = (/ 5 /)
    else if (prob_flag == 4) then
        dimsa_L = (/ 4 /)
        dimsa_k = (/ 9 /)
    end if
    dimsa_f_UC = (/ LC_dim2 /)
    call h5screate_simple_f(1, dimsa_PT, aspace_PT_id, hdf5_ierr)
    call h5screate_simple_f(1, dimsa_N, aspace_N_id, hdf5_ierr)
    if (prob_flag == 1) then
        call h5screate_simple_f(1, dimsa_G, aspace_G_id, hdf5_ierr)
        call h5screate_simple_f(1, dimsa_L, aspace_L_id, hdf5_ierr)
    else if (prob_flag == 4) then
        call h5screate_simple_f(1, dimsa_L, aspace_L_id, hdf5_ierr)
    end if
    call h5screate_simple_f(1, dimsa_k, aspace_k_id, hdf5_ierr)
    call h5screate_simple_f(1, dimsa_f_UC, aspace_f_UC_id, hdf5_ierr)
    call h5acreate_f(file_id, "ProblemType", H5T_NATIVE_INTEGER, aspace_PT_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, prob_flag, dimsa_PT, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "N", H5T_NATIVE_INTEGER, aspace_N_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, N_global, dimsa_N, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    if (prob_flag == 1) then
        call h5acreate_f(file_id, "G", H5T_NATIVE_REAL, aspace_G_id, attr_id, hdf5_ierr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, G, dimsa_G, hdf5_ierr)
        call h5aclose_f(attr_id, hdf5_ierr)
        call h5acreate_f(file_id, "L", H5T_NATIVE_REAL, aspace_L_id, attr_id, hdf5_ierr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, L, dimsa_L, hdf5_ierr)
        call h5aclose_f(attr_id, hdf5_ierr)
    else if (prob_flag == 4) then
        call h5acreate_f(file_id, "L", H5T_NATIVE_REAL, aspace_L_id, attr_id, hdf5_ierr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, L, dimsa_L, hdf5_ierr)
        call h5aclose_f(attr_id, hdf5_ierr)
    end if
    call h5acreate_f(file_id, "k", H5T_NATIVE_REAL, aspace_k_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, k, dimsa_k, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "f_UC", H5T_NATIVE_INTEGER, aspace_f_UC_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, LC(1:3*LC_dim2:3), dimsa_f_UC, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5acreate_f(file_id, "f_DoF", H5T_NATIVE_INTEGER, aspace_f_UC_id, attr_id, hdf5_ierr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, LC(2:3*LC_dim2:3), dimsa_f_UC, hdf5_ierr)
    call h5aclose_f(attr_id, hdf5_ierr)
    call h5sclose_f(aspace_N_id, hdf5_ierr)
    if (prob_flag == 1) then
        call h5sclose_f(aspace_G_id, hdf5_ierr)
        call h5sclose_f(aspace_L_id, hdf5_ierr)
    else if (prob_flag == 4) then
        call h5sclose_f(aspace_L_id, hdf5_ierr)
    end if
    call h5sclose_f(aspace_k_id, hdf5_ierr)
    call h5sclose_f(aspace_f_UC_id, hdf5_ierr)

    if (prob_flag == 1) then
        dims_m = (/ 2*N_global /)
    else if (prob_flag == 2 .or. prob_flag ==4) then
        dims_m = (/ DoF*N_global /)
    end if
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
    call h5dwrite_f(dset_u_id, H5T_NATIVE_DOUBLE, x(1:noState*N_global:2), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)
    call h5dwrite_f(dset_udot_id, H5T_NATIVE_DOUBLE, x(2:noState*N_global:2), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)

    deallocate(m_global, STAT=f_stat)
    deallocate(b_global, STAT=f_stat)
end if

!! Iterate for each time step dt
do 
    !! Calculate k1
    p_loc = 0.
    do it = 1, LC_loc_dim2
        call calc_load(LC_loc(1,it), LC_val_loc(1,it), 0., t, dt, L, &
            & p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) ), &
            & x_loc( noState*LC_loc(1,it)+2*LC_loc(2,it)-1 ) )
        f_val_loc(it) = p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) )
    end do
            
    do it2 = 1, N_loc
        if (it2 == 1 .and. procID == 0) then
            if (prob_flag == 1) then
                call calc_f_pendula(BC(1), 0, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                & G, k, L, x_loc(2*it2+1), k1(2*it2+1))
            else if (prob_flag == 2) then
                call calc_f_phi4(BC(1), 0, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                & k, x_loc(2*it2+1), k1(2*it2+1))
            else if (prob_flag == 4) then
                call calc_f_metabeam(BC(1), 0, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                & k, L, x_loc(12*it2+1), k1(noState*it2+1)) 
                if ( LC_loc_dim2 > 0) then
                    do it3 = 1, LC_loc_dim2
                        if ( it2 == LC_loc(1,it3) ) then
                            if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                k1(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                k1(noState*it2+2*LC_loc(2,it3)) = 0.
                            end if
                        end if
                    end do
                end if
            end if
        else if (it2 == N_loc .and. procID == noProc-1) then
            if (prob_flag == 1) then
                call calc_f_pendula(BC(2), -2, 1, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                & G, k, L, x_loc(2*it2-1), k1(2*it2+1))
            else if (prob_flag == 2) then
                call calc_f_phi4(BC(2), -2, 1, m_loc(it2), b_loc(it2), p_loc(it2), &
                & k, x_loc(2*it2-1), k1(2*it2+1))
            else if (prob_flag == 4) then
                call calc_f_metabeam(BC(2), -12, 11, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                & k, L, x_loc(12*it2-11), k1(noState*it2+1))
                if ( LC_loc_dim2 > 0) then
                    do it3 = 1, LC_loc_dim2
                        if ( it2 == LC_loc(1,it3) ) then
                            if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                k1(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                k1(noState*it2+2*LC_loc(2,it3)) = 0.
                            end if
                        end if
                    end do
                end if
            end if
        else
            if (prob_flag == 1) then
                call calc_f_pendula(BC(0), -2, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                & G, k, L, x_loc(2*it2-1), k1(2*it2+1))
            else if (prob_flag == 2) then
                call calc_f_phi4(BC(0), -2, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                & k, x_loc(2*it2-1), k1(2*it2+1))
            else if (prob_flag == 4) then
                call calc_f_metabeam(BC(0), -12, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                & k, L, x_loc(12*it2-11), k1(noState*it2+1))
                if ( LC_loc_dim2 > 0) then
                    do it3 = 1, LC_loc_dim2
                        if ( it2 == LC_loc(1,it3) ) then
                            if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                k1(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                k1(noState*it2+2*LC_loc(2,it3)) = 0.
                            end if
                        end if
                    end do
                end if
            end if
        end if
    end do

    if ( MOD(t,dtWrite) < 0.1*dt .or. dtWrite-MOD(t,dtWrite) < 0.1*dt) then
        call MPI_GATHERV(k1(noState+1), noState*N_loc, MPI_DOUBLE_PRECISION, k1_global, scatter_sc, scatter_disp, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
        call MPI_GATHERV(f_val_loc(1), LC_loc_dim2, MPI_DOUBLE_PRECISION, f_val, scatter_sc3, scatter_disp3, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
        if (procID == 0) then
            call h5dwrite_f(dset_uddot_id, H5T_NATIVE_DOUBLE, k1_global(2:noState*N_global:2), dims_u, hdf5_ierr, &
            & mspace_u_id, dspace_u_id)
            call h5dwrite_f(dset_f_id, H5T_NATIVE_DOUBLE, f_val, dims_f, hdf5_ierr, mspace_f_id, dspace_f_id)
        end if
    end if
    if ( ABS(t-tEnd) < 0.1*dt ) EXIT

    if (sol_flag == 10) then
        x_loc = x_loc + dt*b1*k1
    else
        !! Update ghost cells
        if (noProc /= 1) then
            if (procID == 0) then
                call MPI_SEND(k1(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
                & MPI_COMM_WORLD, mpi_ierr)
                call MPI_RECV(k1(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, &
                & MPI_COMM_WORLD, status, mpi_ierr)
            else if (procID == noProc-1) then
                call MPI_SEND(k1(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
                & MPI_COMM_WORLD, mpi_ierr)
                call MPI_RECV(k1(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, &
                & MPI_COMM_WORLD, status, mpi_ierr)
            else
                call MPI_SENDRECV(k1(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
                & k1(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, mpi_ierr)
                call MPI_SENDRECV(k1(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
                & k1(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, mpi_ierr)
            end if
        end if

        !! Calculate k2
        xTmp = x_loc + dt*a21*k1

        p_loc = 0.
        do it = 1, LC_loc_dim2
            call calc_load(LC_loc(1,it), LC_val_loc(1,it), c2, t, dt, L, &
                & p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) ), &
                & xTmp( noState*LC_loc(1,it)+2*LC_loc(2,it)-1 ) )
        end do
            
        do it2 = 1, N_loc
            if (it2 == 1 .and. procID == 0) then
                if (prob_flag == 1) then
                    call calc_f_pendula(BC(1), 0, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                    & G, k, L, xTmp(2*it2+1), k2(2*it2+1))
                else if (prob_flag == 2) then
                    call calc_f_phi4(BC(1), 0, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                    & k, xTmp(2*it2+1), k2(2*it2+1))
                else if (prob_flag == 4) then
                    call calc_f_metabeam(BC(1), 0, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                    & k, L, xTmp(12*it2+1), k2(noState*it2+1)) 
                    if ( LC_loc_dim2 > 0) then
                        do it3 = 1, LC_loc_dim2
                            if ( it2 == LC_loc(1,it3) ) then
                                if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                    k2(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                    k2(noState*it2+2*LC_loc(2,it3)) = 0.
                                end if
                            end if
                        end do
                    end if
                end if
            else if (it2 == N_loc .and. procID == noProc-1) then
                if (prob_flag == 1) then
                    call calc_f_pendula(BC(2), -2, 1, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                    & G, k, L, xTmp(2*it2-1), k2(2*it2+1))
                else if (prob_flag == 2) then
                    call calc_f_phi4(BC(2), -2, 1, m_loc(it2), b_loc(it2), p_loc(it2), &
                    & k, xTmp(2*it2-1), k2(2*it2+1))
                else if (prob_flag == 4) then
                    call calc_f_metabeam(BC(2), -12, 11, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                    & k, L, xTmp(12*it2-11), k2(noState*it2+1))
                    if ( LC_loc_dim2 > 0) then
                        do it3 = 1, LC_loc_dim2
                            if ( it2 == LC_loc(1,it3) ) then
                                if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                    k2(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                    k2(noState*it2+2*LC_loc(2,it3)) = 0.
                                end if
                            end if
                        end do
                    end if
                end if
            else
                if (prob_flag == 1) then
                    call calc_f_pendula(BC(0), -2, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                    & G, k, L, xTmp(2*it2-1), k2(2*it2+1))
                else if (prob_flag == 2) then
                    call calc_f_phi4(BC(0), -2, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                    & k, xTmp(2*it2-1), k2(2*it2+1))
                else if (prob_flag == 4) then
                    call calc_f_metabeam(BC(0), -12, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                    & k, L, xTmp(12*it2-11), k2(noState*it2+1))
                    if ( LC_loc_dim2 > 0) then
                        do it3 = 1, LC_loc_dim2
                            if ( it2 == LC_loc(1,it3) ) then
                                if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                    k2(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                    k2(noState*it2+2*LC_loc(2,it3)) = 0.
                                end if
                            end if
                        end do
                    end if
                end if
            end if
        end do

        if (sol_flag == 20) then
            x_loc = x_loc + dt*(b1*k1 + b2*k2)
        else
            !! Update ghost cells
            if (noProc /= 1) then
                if (procID == 0) then
                    call MPI_SEND(k2(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
                    & MPI_COMM_WORLD, mpi_ierr)
                    call MPI_RECV(k2(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, &
                    & MPI_COMM_WORLD, status, mpi_ierr)
                else if (procID == noProc-1) then
                    call MPI_SEND(k2(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
                    & MPI_COMM_WORLD, mpi_ierr)
                    call MPI_RECV(k2(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, &
                    & MPI_COMM_WORLD, status, mpi_ierr)
                else
                    call MPI_SENDRECV(k2(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
                    & k2(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, mpi_ierr)
                    call MPI_SENDRECV(k2(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
                    & k2(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, mpi_ierr)
                end if
            end if

            !! Calculate k3
            xTmp = x_loc + dt*(a31*k1 + a32*k2)

            p_loc = 0.
            do it = 1, LC_loc_dim2
                call calc_load(LC_loc(1,it), LC_val_loc(1,it), c3, t, dt, L, &
                    & p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) ), &
                    & xTmp( noState*LC_loc(1,it)+2*LC_loc(2,it)-1 ) )
            end do
            
            do it2 = 1, N_loc
                if (it2 == 1 .and. procID == 0) then
                    if (prob_flag == 1) then
                        call calc_f_pendula(BC(1), 0, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                        & G, k, L, xTmp(2*it2+1), k3(2*it2+1))
                    else if (prob_flag == 2) then
                        call calc_f_phi4(BC(1), 0, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                        & k, xTmp(2*it2+1), k3(2*it2+1))
                    else if (prob_flag == 4) then
                        call calc_f_metabeam(BC(1), 0, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                        & k, L, xTmp(12*it2+1), k3(noState*it2+1)) 
                        if ( LC_loc_dim2 > 0) then
                            do it3 = 1, LC_loc_dim2
                                if ( it2 == LC_loc(1,it3) ) then
                                    if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                        k3(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                        k3(noState*it2+2*LC_loc(2,it3)) = 0.
                                    end if
                                end if
                            end do
                        end if
                    end if
                else if (it2 == N_loc .and. procID == noProc-1) then
                    if (prob_flag == 1) then
                        call calc_f_pendula(BC(2), -2, 1, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                        & G, k, L, xTmp(2*it2-1), k3(2*it2+1))
                    else if (prob_flag == 2) then
                        call calc_f_phi4(BC(2), -2, 1, m_loc(it2), b_loc(it2), p_loc(it2), &
                        & k, xTmp(2*it2-1), k3(2*it2+1))
                    else if (prob_flag == 4) then
                        call calc_f_metabeam(BC(2), -12, 11, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                        & k, L, xTmp(12*it2-11), k3(noState*it2+1))
                        if ( LC_loc_dim2 > 0) then
                            do it3 = 1, LC_loc_dim2
                                if ( it2 == LC_loc(1,it3) ) then
                                    if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                        k3(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                        k3(noState*it2+2*LC_loc(2,it3)) = 0.
                                    end if
                                end if
                            end do
                        end if
                    end if
                else
                    if (prob_flag == 1) then
                        call calc_f_pendula(BC(0), -2, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                        & G, k, L, xTmp(2*it2-1), k3(2*it2+1))
                    else if (prob_flag == 2) then
                        call calc_f_phi4(BC(0), -2, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                        & k, xTmp(2*it2-1), k3(2*it2+1))
                    else if (prob_flag == 4) then
                        call calc_f_metabeam(BC(0), -12, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                        & k, L, xTmp(12*it2-11), k3(noState*it2+1))
                        if ( LC_loc_dim2 > 0) then
                            do it3 = 1, LC_loc_dim2
                                if ( it2 == LC_loc(1,it3) ) then
                                    if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                        k3(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                        k3(noState*it2+2*LC_loc(2,it3)) = 0.
                                    end if
                                end if
                            end do
                        end if
                    end if
                end if
            end do

            if (sol_flag == 30) then
                x_loc = x_loc + dt*(b1*k1 + b2*k2 + b3*k3)
            else
                !! Update ghost cells
                if (noProc /= 1) then
                    if (procID == 0) then
                        call MPI_SEND(k3(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
                        & MPI_COMM_WORLD, mpi_ierr)
                        call MPI_RECV(k3(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, &
                        & MPI_COMM_WORLD, status, mpi_ierr)
                    else if (procID == noProc-1) then
                        call MPI_SEND(k3(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
                        & MPI_COMM_WORLD, mpi_ierr)
                        call MPI_RECV(k3(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, &
                        & MPI_COMM_WORLD, status, mpi_ierr)
                    else
                        call MPI_SENDRECV(k3(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
                        & k3(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, mpi_ierr)
                        call MPI_SENDRECV(k3(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
                        & k3(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, mpi_ierr)
                    end if
                end if

                !! Calculate k4
                xTmp = x_loc + dt*(a41*k1 + a42*k2 + a43*k3)

                p_loc = 0.
                do it = 1, LC_loc_dim2
                    call calc_load(LC_loc(1,it), LC_val_loc(1,it), c4, t, dt, L, &
                        & p_loc( DoF*(LC_loc(1,it)-1) + LC_loc(2,it) ), &
                        & xTmp( noState*LC_loc(1,it)+2*LC_loc(2,it)-1 ) )
                end do
            
                do it2 = 1, N_loc
                    if (it2 == 1 .and. procID == 0) then
                        if (prob_flag == 1) then
                            call calc_f_pendula(BC(1), 0, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                            & G, k, L, xTmp(2*it2+1), k4(2*it2+1))
                        else if (prob_flag == 2) then
                            call calc_f_phi4(BC(1), 0, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                            & k, xTmp(2*it2+1), k4(2*it2+1))
                        else if (prob_flag == 4) then
                            call calc_f_metabeam(BC(1), 0, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                            & k, L, xTmp(12*it2+1), k4(noState*it2+1)) 
                            if ( LC_loc_dim2 > 0) then
                                do it3 = 1, LC_loc_dim2
                                    if ( it2 == LC_loc(1,it3) ) then
                                        if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                            k4(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                            k4(noState*it2+2*LC_loc(2,it3)) = 0.
                                        end if
                                    end if
                                end do
                            end if
                        end if
                    else if (it2 == N_loc .and. procID == noProc-1) then
                        if (prob_flag == 1) then
                            call calc_f_pendula(BC(2), -2, 1, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                            & G, k, L, xTmp(2*it2-1), k4(2*it2+1))
                        else if (prob_flag == 2) then
                            call calc_f_phi4(BC(2), -2, 1, m_loc(it2), b_loc(it2), p_loc(it2), &
                            & k, xTmp(2*it2-1), k4(2*it2+1))
                        else if (prob_flag == 4) then
                            call calc_f_metabeam(BC(2), -12, 11, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                            & k, L, xTmp(12*it2-11), k4(noState*it2+1))
                            if ( LC_loc_dim2 > 0) then
                                do it3 = 1, LC_loc_dim2
                                    if ( it2 == LC_loc(1,it3) ) then
                                        if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                            k4(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                            k4(noState*it2+2*LC_loc(2,it3)) = 0.
                                        end if
                                    end if
                                end do
                            end if
                        end if
                    else
                        if (prob_flag == 1) then
                            call calc_f_pendula(BC(0), -2, 3, m_loc(2*it2-1), b_loc(it2), p_loc(it2), &
                            & G, k, L, xTmp(2*it2-1), k4(2*it2+1))
                        else if (prob_flag == 2) then
                            call calc_f_phi4(BC(0), -2, 3, m_loc(it2), b_loc(it2), p_loc(it2), &
                            & k, xTmp(2*it2-1), k4(2*it2+1))
                        else if (prob_flag == 4) then
                            call calc_f_metabeam(BC(0), -12, 23, m_loc(6*it2-5), b_loc(6*it2-5), p_loc(6*it2-5), &
                            & k, L, xTmp(12*it2-11), k4(noState*it2+1))
                            if ( LC_loc_dim2 > 0) then
                                do it3 = 1, LC_loc_dim2
                                    if ( it2 == LC_loc(1,it3) ) then
                                        if ( LC_loc(3,it3) == 11 .or. LC_loc(3,it3) == 12 .or. LC_loc(3,it3) == 13) then
                                            k4(noState*it2+2*LC_loc(2,it3)-1) = 0.
                                            k4(noState*it2+2*LC_loc(2,it3)) = 0.
                                        end if
                                    end if
                                end do
                            end if
                        end if
                    end if
                end do
                x_loc = x_loc + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
            end if
        end if
    end if

    !! Update ghost cells
    if (noProc /= 1) then
        if (procID == 0) then
            call MPI_SEND(x_loc(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
            & MPI_COMM_WORLD, mpi_ierr)
            call MPI_RECV(x_loc(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, &
            & MPI_COMM_WORLD, status, mpi_ierr)
        else if (procID == noProc-1) then
            call MPI_SEND(x_loc(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
            & MPI_COMM_WORLD, mpi_ierr)
            call MPI_RECV(x_loc(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, &
            & MPI_COMM_WORLD, status, mpi_ierr)
        else
            call MPI_SENDRECV(x_loc(noState+1), noState, MPI_DOUBLE_PRECISION, procID-1, 0, &
            & x_loc(1), noState, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, mpi_ierr)
            call MPI_SENDRECV(x_loc(noState*N_loc+1), noState, MPI_DOUBLE_PRECISION, procID+1, 1, &
            & x_loc(noState*N_loc+noState+1), noState, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, mpi_ierr)
        end if
    end if

    t = t + dt

    if ( MOD(t,dtWrite) < 0.1*dt .or. dtWrite-MOD(t,dtWrite) < 0.1*dt ) then
        call MPI_GATHERV(x_loc(noState+1), noState*N_loc, MPI_DOUBLE_PRECISION, x, scatter_sc, scatter_disp, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

        if (procID == 0) then
            if ( SUM(x) /= SUM(x) ) then ! Check if NaN (overflow)
                write (*,'(A)') " >>"
                write (*,'(A)') " >> WARNING: NaN values in output displacements!!!" 
                write (*,'(A)') " >> This is most likely a numerical stability issue."
                write (*,'(A)') " >> Use a smaller time step or the NB solver."

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
                call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
            end if
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
            call h5dwrite_f(dset_u_id, H5T_NATIVE_DOUBLE, x(1:noState*N_global:2), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)
            call h5dwrite_f(dset_udot_id, H5T_NATIVE_DOUBLE, x(2:noState*N_global:2), dims_u, hdf5_ierr, mspace_u_id, dspace_u_id)
        end if
    end if
end do
        
call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
tWallclockEnd = MPI_WTIME()

if (procID == 0) then
    write (*, '(A)') " >>"
    write (*, '(A)') " >> Simulation has completed."
    write (*,'(A,F15.6,A)') " >> Wallclock time: ", tWallclockEnd-tWallclockStart, " s."
    write (*, '(A/)') " >> Results are saved as '"//trim(filename_h5)//"' in ./outputs folder."

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

call MPI_FINALIZE(mpi_ierr)

end program main_rk
