program main_NB
!! Purpose: To find the dynamic response of a metabeam of bistable elements
!!         using Newmark-beta method
!!
!! Author: Written by Myungwon Hwang (hwang125@purdue.edu)
!!
!! Using Makefile:
!!   To compile: make
!!   To run: make run NP=(No of Procs) INP=(InputFile)
!!
!! Record of revisions:
!!    DATE        Programmer                      Decsription
!! ==========  ================  =================================================
!! 11/03/2020   Myungwon Hwang   Rev00: Initial working program
!! 08/10/2021   Myungwon Hwang   Rev01: Syntax update to match that of RK solver
!!

use phdf5_helpers
use mpi
use sub_NB
        
implicit none
!include 'mpif.h' ! if 'use mpi' does not work

!! Data dictionary: Constants
real, parameter :: PI=3.141592653589793
!! Data dictionary: MPI-related
integer :: mpi_ierr, noProc, procID, status(MPI_STATUS_SIZE), mpi_errcode
integer, allocatable, dimension(:) :: scatter_sc, scatter_sc2 ! No. of data to be sent for each process for MPI_SCATTERV
integer, allocatable, dimension(:) :: scatter_disp, scatter_disp2 ! The starting counter for MPI_SCATTERV
real :: tWallclockStart, tWallclockEnd ! wall clock time of the start and the end of the execution
!! Data dictionary: hdf5-related
character(88) :: filename
integer(HID_T) :: file_id, dspace_id, aspace_id, dset_id, attr_id
integer(HSIZE_T), dimension(1:1) :: dimsf, dimsa, data_dims
!integer(HID_T) :: plist_id
integer :: hdf5_ierr
!integer(hsize_t), dimension(7) :: offset
character(len=80) :: dsetname
!! Data dictionary: input parameters
integer :: N_global ! Total No. of unit cells
integer :: prob_flag ! 1: pendula chain, 2: phi-4 lattice, 4: metabeam
integer :: DoF ! unit cell degrees of freedom
real, allocatable, dimension(:) :: m_global ! global mass array
real, allocatable, dimension(:) :: b_global ! global damping array
integer, dimension(0:2) :: BC ! boundary condition flag
real, dimension(9) :: k ! spring stiffness array
real, dimension(4) :: L ! length array containing unit cell geometry information
real :: zeta ! damping coefficient
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
real :: gamma_NB ! for linear acceleration method
real :: beta_NB ! for linear acceleration method
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
real, allocatable, dimension(:) :: a1_loc, a2_loc, a3_loc, m_loc ! corresponding local array
!! Data dictionary: temporary computational variables
integer :: it, it2, ind
integer :: cnt1
real, dimension(6) :: fs
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
    open(unit=1, file='numMethod_NB.inp', status='old', action='read', iostat=f_ierr, iomsg=f_msg)
    
    if (f_ierr == 0) then
        write (*,'(A)') "numMethod_NB.inp is opened successfully."
        open(unit=2, file='design.inp', status='old', action='read', iostat=f_ierr, iomsg=f_msg)
        if (f_ierr == 0) then
            write (*,'(A)') "design.inp is opened successfully."

            120 format (T3, ' Increment ', T18, 'NR iterations', T33 'Inv iterations')
            121 format (T3, '===========', T18, '=============', T33 '==============')
            122 format (T3, I8, T18, I8, T33, I8)
            123 format (T4, "The total execution time is ", F14.8, " seconds.")
            124 format (T4, "The total time for inverse calculation: ", F14.8, " seconds.")
            125 format (T4, "The total time for NR iterations: ", F14.8, " seconds.")
            open(unit=12, file='main_NB.stat', status='replace', action='write', iostat=f_ierr, iomsg=f_msg)
            write (12,120)
            write (12,121)

            read (1,*); read (1,*) gamma_NB
            read (1,*); read (1,*) beta_NB
            read (1,*); read (1,*) tol_NR
            read (1,*); read (1,*) tol_inv
            read (1,*); read (1,*) N_NR
            read (1,*); read (1,*) N_inv
            read (1,*); read (1,*) dt
            read (1,*); read (1,*) tStart
            read (1,*); read (1,*) tEnd
            read (1,*); read (1,*) dtWrite

            read (2,'(/A)') filename
            call h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, hdf5_ierr)
            read (2,'(/I10)') N_global
            if (N_global < noProc) then
                write (*,'(A)') "The number of unit cells is less than the number of processes. &
                & Use NP value less than the number of unit cells."
                call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
            end if
            read (2,'(/I10)') prob_flag
            read (2,*); read (2,*) BC(0:2)
            if (prob_flag == 1) then
                read (2,*); read (2,*) k(1:3)
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

                allocate(a1_global(DoF*N_global))
                allocate(a2_global(DoF*N_global))
                allocate(a3_global(DoF*N_global))

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

                allocate(a1_global(DoF*N_global))
                allocate(a2_global(DoF*N_global))
                allocate(a3_global(DoF*N_global))

            else if (prob_flag == 4) then
                read (2,*); read (2,*) k(1:8)
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

                allocate(a1_global(DoF*N_global))
                allocate(a2_global(DoF*N_global))
                allocate(a3_global(DoF*N_global))
            end if
        else
            write (*,'(A)') f_msg
            call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
        end if
    else
        write (*,'(A)') f_msg
        call MPI_ABORT(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
    end if
end if

call MPI_BCAST(gamma_NB, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(beta_NB, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tol_NR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tol_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(N_NR, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(N_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tStart, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(tEnd, 1, MPI_DOUBLE_PRECISION , 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(dtWrite, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
Nt = NINT((tEnd - tStart)/dt)
dNt_output = NINT(dtWrite/dt)

call MPI_BCAST(N_global, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(prob_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
call MPI_BCAST(BC(0), 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
if (prob_flag == 1) then
    call MPI_BCAST(k(1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(L(1), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (procID /= 0) then
        allocate(LC(3*LC_dim2), STAT=f_stat)
        allocate(LC_val(7*LC_dim2), STAT=f_stat)
    end if
    call MPI_BCAST(LC(1), 3*LC_dim2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_val(1), 7*LC_dim2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(DoF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
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
else if (prob_flag == 4) then
    call MPI_BCAST(k(1), 8, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(L(1), 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_dim2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    if (procID /= 0) then
        allocate(LC(3*LC_dim2), STAT=f_stat)
        allocate(LC_val(7*LC_dim2), STAT=f_stat)
    end if
    call MPI_BCAST(LC(1), 3*LC_dim2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(LC_val(1), 7*LC_dim2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)
    call MPI_BCAST(DoF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)
end if

N_fl = N_global/noProc
N_ceil = N_fl + MOD(N_global, noProc)
if (procID == noProc-1) then
    N_loc = N_ceil
else
    N_loc = N_fl
end if
write (*,*) "proc ", procID, "has ", N_loc, "unit cells."

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
    scatter_disp(1) = 0
    scatter_disp2(1) = 0
    do it = 1, noProc
        if (it == noProc) then
            scatter_sc(it) = DoF*(N_global/noProc + MOD(N_global,noProc))
            scatter_sc2(it) = LC_loc_dim2
        else
            scatter_sc(it) = DoF*(N_global/noProc)
            scatter_disp(it+1) = scatter_disp(it) + scatter_sc(it)
            scatter_sc2(it) = LC_loc_dim2
            scatter_disp2(it+1) = scatter_disp2(it) + scatter_sc2(it)
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
            if (prob_flag == 1) then
                call calc_fs_pendula(BC(1), 0, 1, u(it), k, L, fs)
            else if (prob_flag == 2) then
                call calc_fs_phi4(BC(1), 0, 1, u(it), k, fs)
            else if (prob_flag == 4) then
                call calc_fs(BC(1), 0, 11, u(6*it-5), k, L, fs)
            end if
        else if (it == N_global) then
            if (prob_flag == 1) then
                call calc_fs_pendula(BC(2), -1, 0, u(it-1), k, L, fs)
            else if (prob_flag == 2) then
                call calc_fs_phi4(BC(2), -1, 0, u(it-1), k, fs)
            else if (prob_flag == 4) then
                call calc_fs(BC(2), -6, 5, u(6*it-11), k, L, fs)
            end if
        else
            if (prob_flag == 1) then
                call calc_fs_pendula(BC(0), -1, 1, u(it-1), k, L, fs)
            else if (prob_flag == 2) then
                call calc_fs_phi4(BC(0), -1, 1, u(it-1), k, fs)
            else if (prob_flag == 4) then
                call calc_fs(BC(0), -6, 11, u(6*it-11), k, L, fs)
            end if
        end if

        do it2=1, DoF
            ind = DoF*it-DoF+it2
! Need to be modified if nonzero external force p
            uddot(ind) = (-b_global(ind)*udot(ind)-fs(it2))/m_global(ind)
        end do
    end do

! calculation for a, b, khat arrays
    a1_global = m_global/(beta_NB*dt**2) + b_global*gamma_NB/(beta_NB*dt)
    a2_global = m_global/(beta_NB*dt) + b_global*(gamma_NB/beta_NB-1)
    a3_global = m_global*(1/(2*beta_NB)-1) + b_global*dt*(gamma_NB/(2*beta_NB)-1)
! end STEP 1
end if

deallocate(m_global, STAT=f_stat)
deallocate(b_global, STAT=f_stat)

allocate(u_loc(DoF*(N_loc+2)))
allocate(uOld_loc(DoF*(N_loc+2)))
allocate(udot_loc(DoF*N_loc))
allocate(uddot_loc(DoF*N_loc))
allocate(p_loc(DoF*N_loc))
allocate(phat_loc(DoF*N_loc))
allocate(a1_loc(DoF*N_loc))
allocate(a2_loc(DoF*N_loc))
allocate(a3_loc(DoF*N_loc))

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

if ( procID == 0 ) then
    dimsf = (/ DoF*N_global /)
    dimsa = (/ 1 /)
    data_dims = (/ DoF*N_global /)
end if

! Iterations for each time step dt
do it = 1, Nt
    p_loc = 0.
    phat_loc = 0.
    do it2 = 1, LC_loc_dim2
        if ( LC_loc(3,it2) == 1) then
            if ( (it*dt) >= LC_val_loc(1,it2) .and. (it*dt) <= LC_val_loc(2,it2) ) then
                f_val_loc(it2) = LC_val_loc(3,it2)*SIN(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))
                p_loc( DoF*(LC_loc(1,it2)-1) + LC_loc(2,it2) ) = f_val_loc(it2)
            end if
        else if (LC_loc(3,it2) == 2) then
            if ( (it*dt) >= LC_val_loc(1,it2) .and. (it*dt) <= LC_val_loc(2,it2) ) then
                f_val_loc(it2) = LC_val_loc(3,it2)*SIN(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))*&
                & SIN(PI*((it*dt)-LC_val_loc(1,it2))/(LC_val_loc(2,it2)-LC_val_loc(1,it2)))**2
                p_loc( DoF*(LC_loc(1,it2)-1) + LC_loc(2,it2) ) = f_val_loc(it2)
            end if
        else if (LC_loc(3,it2) == 3) then
        else if (LC_loc(3,it2) == 11) then !!! To be implemented
            if ( (it*dt) >= LC_val_loc(1,it2) .and. (it*dt) <= LC_val_loc(2,it2) ) then
                u_loc(DoF*LC_loc(1,it2)+LC_loc(2,it2)-1) = LC_val_loc(3,it2)*&
                & SIN(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))
                u_loc(DoF*LC_loc(1,it2)+LC_loc(2,it2)) = 2*PI*LC_val_loc(4,it2)*&
                & LC_val_loc(3,it2)*COS(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))
            end if
        else if (LC_loc(3,it2) == 12) then !!! To be implemented
            if ( (it*dt) >= LC_val_loc(1,it2) .and. (it*dt) <= LC_val_loc(2,it2) ) then
                u_loc(DoF*LC_loc(1,it2)+LC_loc(2,it2)-1) = LC_val_loc(3,it2)*&
                & SIN(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))*&
                & SIN(PI*((it*dt)-LC_val_loc(1,it2))/(LC_val_loc(2,it2)-LC_val_loc(1,it2)))**2
                u_loc(DoF*LC_loc(1,it2)+LC_loc(2,it2)) = 2*LC_val_loc(3,it2)*PI*&
                & SIN(PI*((it*dt)-LC_val_loc(1,it2))/(LC_val_loc(2,it2)-LC_val_loc(1,it2)))*&
                & ( LC_val_loc(4,it2)*COS(2*PI*LC_val_loc(4,it2)*(it*dt)+LC_val_loc(5,it2))*&
                & SIN(PI*((it*dt)-LC_val_loc(1,it2))/(LC_val_loc(2,it2)-LC_val_loc(1,it2)))+&
                & COS(PI*((it*dt)-LC_val_loc(1,it2))/(LC_val_loc(2,it2)-LC_val_loc(1,it2)))*&
                & SIN(2*PI*LC_val_loc(4,it2)*(it*dt)+LC_val_loc(5,it2))/(LC_val_loc(2,it2)-LC_val_loc(1,it2)) )
            end if
        else if (LC_loc(3,it2) == 13) then !!! To be implemented
            if ( (it*dt) >= LC_val_loc(1,it2) .and. (it*dt) <= LC_val_loc(2,it2) ) then
                u_loc(DoF*LC_loc(1,it2)+LC_loc(2,it2)-1) = L(4) - LC_val_loc(3,it2)*&
                & COS(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))
                u_loc(DoF*LC_loc(1,it2)+LC_loc(2,it2)) = 2*PI*LC_val_loc(4,it2)*&
                & LC_val_loc(3,it2)*SIN(2*PI*LC_val_loc(4,it2)*(it*dt) + LC_val_loc(5,it2))
            end if
        end if
    end do 
    do it2 = 1, DoF*N_loc
        phat_loc(it2) = p_loc(it2) + a1_loc(it2)*u_loc(it2+DoF) + a2_loc(it2)*udot_loc(it2) + &
        & a3_loc(it2)*uddot_loc(it2)
    end do
    uOld_loc = u_loc
    call NR_iterations(N_NR, tol_NR, N_loc, procID, noProc, prob_flag, DoF, BC, k, L, &
    & a1_loc, phat_loc, u_loc, N_inv, tol_inv, it, cnt_NR, cnt_inv, t_NR_total, t_inv_total)

    ind = DoF*(N_loc+1)
    uddot_loc = (u_loc(DoF+1:ind) - uOld_loc(DoF+1:ind))/(beta_NB*dt**2) - udot_loc/(beta_NB*dt) - uddot_loc*(1/(2*beta_NB)-1)
    udot_loc = (u_loc(DoF+1:ind) - uOld_loc(DoF+1:ind))*gamma_NB/(beta_NB*dt) + udot_loc*(1-gamma_NB/beta_NB) &
    & + uddot_loc*dt*(1-gamma_NB/(2*beta_NB))

    if ( MOD(it, dNt_output) == 0) then
        call MPI_GATHERV(u_loc(DoF+1), DoF*N_loc, MPI_DOUBLE_PRECISION, u, scatter_sc, scatter_disp, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_ierr)

        if ( procID == 0 ) then
            call h5screate_simple_f(1, dimsf, dspace_id, hdf5_ierr)
            call h5screate_simple_f(1, dimsa, aspace_id, hdf5_ierr)
            write (tmpStr1,*) it/dNt_output
            dsetname = '/u'//adjustl(tmpStr1)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdf5_ierr)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, u, data_dims, hdf5_ierr)
            call h5acreate_f(dset_id, "t", H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdf5_ierr)
            call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, it*dt, dimsa, hdf5_ierr)
            call h5aclose_f(attr_id, hdf5_ierr)
            call h5dclose_f(dset_id, hdf5_ierr)
            call h5sclose_f(dspace_id, hdf5_ierr)
            call h5sclose_f(aspace_id, hdf5_ierr)
        end if
    end if

    if (procID == 0) then
        write (12,122) it, cnt_NR, cnt_inv
    end if
end do

call MPI_BARRIER(MPI_COMM_WORLD, mpi_ierr)
tWallclockEnd = MPI_WTIME()
if (procID == 0) then
    write (*,*) tWallclockEnd - tWallclockStart
    do it = DoF+1, DoF*(10+1)
        write (*,*) u_loc(it)
    end do
    write (12,121)
    write (12,123) tWallclockEnd - tWallclockStart
    write (12,124) t_inv_total
    write (12,125) t_NR_total - t_inv_total
    call h5fclose_f(file_id, hdf5_ierr)
end if

close(unit=12)


call MPI_FINALIZE(mpi_ierr)

end program main_NB
