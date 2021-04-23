program main_NB
!! Purpose: To find the dynamic response of a metabeam of bistable elements
!!         using Newmark-beta method
!!
!! To compile (gnu): mpif90 -o main.exe -fdefault-real-8 main.f90
!! To compile (intel): mpif90 -o main.exe -r8 main.f90
!! To run : mpirun -np 4 main.exe
!!

use phdf5_helpers
use mpi
use sub_NB
        
implicit none
!include 'mpif.h' ! if 'use mpi' does not work

!! Data dictionary: Constants
real, parameter :: PI=3.141592653589793
!! Data dictionary: MPI-related
integer :: ierr, noProc, procID, status(MPI_STATUS_SIZE)
integer, allocatable, dimension(:) :: scatter_sendcounts ! No. of data to be sent for each process for MPI_SCATTERV
integer, allocatable, dimension(:) :: scatter_disp ! The starting counter for MPI_SCATTERV
real :: t_start, t_end ! wall clock time of the start and the end of the execution
!! Data dictionary: hdf5-related
integer(HID_T) :: file_id, dspace_id, aspace_id, dset_id, attr_id
integer(HSIZE_T), dimension(1:1) :: dimsf, dimsa, data_dims
!integer(HID_T) :: plist_id
integer :: hdf5_ierr
!integer(hsize_t), dimension(7) :: offset
character(len=80) :: dsetname
!! Data dictionary: input parameters
integer :: N_global ! Total No. of unit cells
real, allocatable, dimension(:) :: m_global ! global mass array
real, allocatable, dimension(:) :: c_global ! global damping array
integer, dimension(0:2) :: BC ! boundary condition flag
real, dimension(9) :: k ! spring stiffness array
real, dimension(4) :: L ! length array containing unit cell geometry information
real :: zeta ! damping coefficient
integer :: LC_flag ! Load case (1: single-point horizontal, 2: distributed vertical)
integer :: startnode ! global starting node of the distributed vertical load
integer :: endnode ! global end node of the distributed vertical load
real :: forceIn ! magnitude of input force
real :: freqIn ! input frequency
real :: vFlow ! flow speed
real :: t_UC ! time for the flow to travel one unit length; Or, time phase
real :: dt ! numerical time step
integer :: Nt ! No. of time steps
integer :: dNt_output ! Output results at every `dNt_output' steps
integer :: N_NR ! max number of Newton-Raphson iterations
real :: tol_NR ! convergence criteria for NR method
integer :: N_inv ! max number of iterations for conjugate gradient method
real :: tol_inv ! tolerance for conjugate gradient method
real :: gamma_NB ! for linear acceleration method
real :: beta_NB ! for linear acceleration method
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
real, dimension(6) :: fs
real, allocatable, dimension(:) :: uOld_loc ! displacements of the previous time step
integer :: tmpI1
real :: tmpR1, tmpR2
character(len=20) :: tmpStr1
!! Data dictionary: temporary file I/O
integer :: f10_ierr, f11_ierr, f12_ierr, f13_ierr
character(len=80) :: f10_msg, f11_msg, f12_msg, f13_msg

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, procID, ierr)

if (procID == 0) then
    call h5open_f(hdf5_ierr)
    call h5fcreate_f('main_NB.h5', H5F_ACC_TRUNC_F, file_id, hdf5_ierr)
end if

120 format (T3, ' Increment ', T18, 'NR iterations', T33 'Inv iterations')
121 format (T3, '===========', T18, '=============', T33 '==============')
122 format (T3, I8, T18, I8, T33, I8)
123 format (T4, "The total execution time is ", F14.8, " seconds.")
124 format (T4, "The total time for inverse calculation: ", F14.8, " seconds.")
125 format (T4, "The total time for NR iterations: ", F14.8, " seconds.")
if (procID == 0) then
    open(unit=12, file='main_NB.stat', status='replace', action='write', iostat=f12_ierr, iomsg=f12_msg)
    write (12,120)
    write (12,121)
end if

! Read in the total number of unit cells
if (procID == 0) then
    write (*,*) 'The total number of unit cell: '
    read (*,*) N_global
!    N_global = 30
end if

! Broadcast necessary global variable to all processes
call MPI_BCAST(N_global, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 

! Initialization of global variables
dt = 1.e-4
Nt = 10000
dNt_output = Nt/1000
N_NR = 200
tol_NR = 1.e-16
N_inv = 200
tol_inv = 1.e-16
gamma_NB = 0.5
beta_NB = 0.25

BC(0) = 100 ! 100 - baseline periodic; 200 - bottom fixed periodic; 300 = elastic foundation
BC(1) = 111 ! #10 - left end free; #11 - left end fixed;
BC(2) = 120 ! #20 - right end free; #21 - right end fixed;

k(1) = 1.241
k(2) = 0.6*1.793
k(3) = 0.6
k(4) = 100.
k(5) = 100.
k(6) = 100.
k(7) = 100.
k(8) = 100.
k(9) = 0. ! vertical spring connecting bottom mass to ground

L(1) = 20.
L(2) = 40.
L(3) = 20.
L(4) = 8.

LC_flag = 1
if (LC_flag == 1) then ! horizontal single-point at the center of the first unit cell
    forceIn = 0.1
    freqIn = 8.0
else if (LC_flag == 2) then ! vertical distributed load on the top surface
    startnode = 1
    endnode = N_global/2
    forceIn = -3.
    freqIn = 150.
    vFlow = 100000.
    t_UC = L(1)/vFlow
end if
    
if (procID == 0) then
    allocate(u(6*N_global))
    allocate(udot(6*N_global))
    allocate(uddot(6*N_global))
    allocate(m_global(6*N_global))
    allocate(c_global(6*N_global))
    allocate(a1_global(6*N_global))
    allocate(a2_global(6*N_global))
    allocate(a3_global(6*N_global))

    do it=1, N_global
        m_global(6*it-5) = 2.e-6
        m_global(6*it-4) = 2.e-6
        m_global(6*it-3) = 1.e-6
        m_global(6*it-2) = 1.e-6
        m_global(6*it-1) = 1.e-6
        m_global(6*it) = 1.e-6
    end do

    zeta = 0.03
    do it=1, N_global
        c_global(6*it-5) = 2*zeta*2*PI*26.281*m_global(1)
        c_global(6*it-4) = 2*zeta*2*PI*26.281*m_global(2)
        c_global(6*it-3) = 2*zeta*2*PI*26.281*m_global(3)
        c_global(6*it-2) = 2*zeta*2*PI*26.281*m_global(4)
        c_global(6*it-1) = 2*zeta*2*PI*26.281*m_global(5)
        c_global(6*it) = 2*zeta*2*PI*26.281*m_global(6)
    end do

    u = 0.
    udot = 0.

! STEP 1: Initial calculations
    do it=1, N_global
        if (it == 1) then
            call calc_fs(BC(1), 0, 11, u(6*it-5), k, L, fs)
        else if (it == N_global) then
            call calc_fs(BC(2), -6, 5, u(6*it-11), k, L, fs)
        else
            call calc_fs(BC(0), -6, 11, u(6*it-11), k, L, fs)
        end if

        do it2=1, 6
            ind = 6*it-6+it2
! Need to be modified if nonzero external force p
            uddot(ind) = (-c_global(ind)*udot(ind)-fs(it2))/m_global(ind)
        end do
    end do

! calculation for a, b, khat arrays
    a1_global = m_global/(beta_NB*dt**2) + c_global*gamma_NB/(beta_NB*dt)
    a2_global = m_global/(beta_NB*dt) + c_global*(gamma_NB/beta_NB-1)
    a3_global = m_global*(1/(2*beta_NB)-1) + c_global*dt*(gamma_NB/(2*beta_NB)-1)
! end STEP 1
end if

! Initialization of local variables
if ( MOD(N_global, noProc) == 0 ) then
    N_fl = N_global/noProc
    N_ceil = N_fl
    N_loc = N_fl
else
    N_fl = N_global/noProc
    N_ceil = N_fl + MOD(N_global, noProc)
    if ( procID == noProc-1 ) then
        N_loc = N_ceil
    else
        N_loc = N_fl
    end if
end if

write (*,*) "proc ", procID, "has ", N_loc, "unit cells."

tmpI1 = endnode*noProc/N_global
if ( procID .lt. tmpI1 ) then
    startnode_loc = 1
    endnode_loc = N_loc
else if ( procID == tmpI1 .and. MOD(endnode*noProc,N_global) /= 0 ) then
    startnode_loc = 1
    endnode_loc = MOD(endnode, N_global/noProc)
else
    startnode_loc = 0
    endnode_loc = 0
end if

! set parameters for MPI_SCATTERV
allocate(scatter_sendcounts(noProc))
allocate(scatter_disp(noProc))
scatter_disp(1) = 0
if (procID == 0) then
    do it = 1, noProc
        if (it == noProc) then
            scatter_sendcounts(it) = 6*(N_global/noProc + MOD(N_global,noProc))
        else
            scatter_sendcounts(it) = 6*(N_global/noProc)
            scatter_disp(it+1) = scatter_disp(it) + scatter_sendcounts(it)
        end if
    end do
end if

allocate(u_loc(6*(N_loc+2)))
allocate(uOld_loc(6*(N_loc+2)))
allocate(udot_loc(6*N_loc))
allocate(uddot_loc(6*N_loc))
allocate(p_loc(6*N_loc))
allocate(phat_loc(6*N_loc))
allocate(a1_loc(6*N_loc))
allocate(a2_loc(6*N_loc))
allocate(a3_loc(6*N_loc))

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
t_start = MPI_WTIME()

! Distribute global data to each process
call MPI_SCATTERV(u, scatter_sendcounts, scatter_disp, MPI_DOUBLE_PRECISION, u_loc(7), 6*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_SCATTERV(udot, scatter_sendcounts, scatter_disp, MPI_DOUBLE_PRECISION, udot_loc, 6*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_SCATTERV(uddot, scatter_sendcounts, scatter_disp, MPI_DOUBLE_PRECISION, uddot_loc, 6*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_SCATTERV(a1_global, scatter_sendcounts, scatter_disp, MPI_DOUBLE_PRECISION, a1_loc, 6*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_SCATTERV(a2_global, scatter_sendcounts, scatter_disp, MPI_DOUBLE_PRECISION, a2_loc, 6*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_SCATTERV(a3_global, scatter_sendcounts, scatter_disp, MPI_DOUBLE_PRECISION, a3_loc, 6*N_loc, &
& MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! Update ghost cells
if (noProc /= 1) then
    if (procID == 0) then 
        call MPI_SEND(u_loc(6*N_loc+1), 6, MPI_DOUBLE_PRECISION, procID+1, 1, &
        & MPI_COMM_WORLD, ierr) 
        call MPI_RECV(u_loc(6*N_loc+7), 6, MPI_DOUBLE_PRECISION, procID+1, 0, &
        & MPI_COMM_WORLD, status, ierr) 
    else if (procID == noProc-1) then
        call MPI_SEND(u_loc(7), 6, MPI_DOUBLE_PRECISION, procID-1, 0, &
        & MPI_COMM_WORLD, ierr)
        call MPI_RECV(u_loc(1), 6, MPI_DOUBLE_PRECISION, procID-1, 1, &
        & MPI_COMM_WORLD, status, ierr)
    else
        call MPI_SENDRECV(u_loc(7), 6, MPI_DOUBLE_PRECISION, procID-1, 0, &
        & u_loc(1), 6, MPI_DOUBLE_PRECISION, procID-1, 1, MPI_COMM_WORLD, status, ierr)
        call MPI_SENDRECV(u_loc(6*N_loc+1), 6, MPI_DOUBLE_PRECISION, procID+1, 1, &
        & u_loc(6*N_loc+7), 6, MPI_DOUBLE_PRECISION, procID+1, 0, MPI_COMM_WORLD, status, ierr)
    end if
end if

!offset(2) = procID*6*N_fl
if ( procID == 0 ) then
    dimsf = (/ 6*N_global /)
    dimsa = (/ 1 /)
    data_dims = (/ 6*N_global /)
end if

! Iterations for each time step dt
do it = 1, Nt
    p_loc = 0.
    phat_loc = 0.
    if (LC_FLAG == 1) then
        if (procID == 0) then
            p_loc(1) = forceIn*SIN(2*PI*freqIn*it*dt)
        end if
    else if (LC_FLAG == 2) then
        tmpR1 = it*dt
        do it2 = startnode_loc, endnode_loc
            if (startnode_loc == 0) then
                exit
            else if (procID == noProc-1) then
                tmpR2 = (N_global-N_loc+it2-1)*t_UC
                if ( tmpR1 .gt. tmpR2 ) then
                    p_loc(6*it2-2) = forceIn*SIN(2*PI*freqIn*(tmpR1-tmpR2))
                end if
            else
                tmpR2 = (N_loc*procID+it2-1)*t_UC
                if ( tmpR1 .gt. tmpR2 ) then
                    p_loc(6*it2-2) = forceIn*SIN(2*PI*freqIn*(tmpR1-tmpR2))
                end if
            end if
        end do
    end if
    do it2 = 1, 6*N_loc
        phat_loc(it2) = p_loc(it2) + a1_loc(it2)*u_loc(it2+6) + a2_loc(it2)*udot_loc(it2) + &
        & a3_loc(it2)*uddot_loc(it2)
    end do
    uOld_loc = u_loc
    call NR_iterations(N_NR, tol_NR, N_loc, procID, noProc, BC, k, L, a1_loc, phat_loc, u_loc, &
    & N_inv, tol_inv, it, cnt_NR, cnt_inv, t_NR_total, t_inv_total)

    ind = 6*(N_loc+1)
    uddot_loc = (u_loc(7:ind) - uOld_loc(7:ind))/(beta_NB*dt**2) - udot_loc/(beta_NB*dt) - uddot_loc*(1/(2*beta_NB)-1)
    udot_loc = (u_loc(7:ind) - uOld_loc(7:ind))*gamma_NB/(beta_NB*dt) + udot_loc*(1-gamma_NB/beta_NB) &
    & + uddot_loc*dt*(1-gamma_NB/(2*beta_NB))

    if ( MOD(it, dNt_output) == 0) then
        call MPI_GATHERV(u_loc(7), 6*N_loc, MPI_DOUBLE_PRECISION, u, scatter_sendcounts, scatter_disp, &
        & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

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

call MPI_BARRIER(MPI_COMM_WORLD, ierr)
t_end = MPI_WTIME()
if (procID == 0) then
    write (*,*) t_end - t_start
    do it = 7, 6*(N_loc+1)
        write (*,*) u_loc(it)
    end do
    write (12,121)
    write (12,123) t_end - t_start
    write (12,124) t_inv_total
    write (12,125) t_NR_total - t_inv_total
    call h5fclose_f(file_id, hdf5_ierr)
end if

close(unit=12)


call MPI_FINALIZE(ierr)

end program main_NB
