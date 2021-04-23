! module phdf5_helpers
!
! Convenience subroutines for working with HDF5 files.





module phdf5_helpers

  use hdf5
  use hdf5_helpers
 ! include 'mpif.h'

contains

  !....
  !   subroutine initialise_phdf5(error)
  !   subroutine create_phdf5_file(filename,comm,info,file_id,plist_id,error)
  !   subroutine open_phdf5_file_rdonly(filename,file_id,plist_id,error)
  !   subroutine open_phdf5_file_rdonly(filename,file_id,plist_id,error)
  !   subroutine finalise_phdf5_file(file_id,plist_id,error)
  !
  !   subroutine read_phdf5_scalar_integer(name,a,file_id,error)
  !   subroutine read_phdf5_scalar_double(name,a,file_id,error)
  !   subroutine read_phdf5_vector_double(name,n,a,file_id,error)
  !   subroutine read_phdf5_3darray_double(dsetname,nii,njj,nkk,px,pz,offset,rdata_sub,file_id,error)
  !   subroutine read_phdf5_4darray_double(dsetname,nii,njj,nkk,nll,px,pz,offset,rdata_sub,file_id,error)
  !
  !   subroutine write_phdf5_scalar_integer(name,a,file_id,error)
  !   subroutine write_phdf5_scalar_double(name,a,file_id,error)
  !   subroutine write_phdf5_3darray_double(dsetname,nii,njj,nkk,px,pz,offset,subdata,file_id,error)
  !   subroutine write_phdf5_4darray_double(dsetname,nii,njj,nkk,nll,px,pz,offset,subdata,file_id,error)
  !....




  ! subroutine hdf5_helpers/initialise_hdf5
  !
  ! Startup HDF5 by running h5open_f. Required before calling any
  ! HDF5 routines.

  !!!already present in serial hdf5_helpers
  subroutine initialise_phdf5(error)

    implicit none

    integer, intent(out) :: error                 ! error flag

    call h5open_f(error)

  end subroutine initialise_phdf5





  ! subroutine hdf5_helpers/create_hdf5_file
  !
  ! Create and open a new HDF5 file.generate file_id,plist_id

  subroutine create_phdf5_file(filename,comm,info,file_id,plist_id,error)

    implicit none
!    include 'mpif.h'

    character(len=*), intent(in) :: filename      ! name of HDF5 file to create
    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
    integer(hid_t), intent(out)  :: plist_id      ! Property list identifier
    integer, intent(in)          :: comm, info
    integer, intent(inout)       :: error

     !
     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
     !
     ! Create the file collectively.
     !
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)

  end subroutine create_phdf5_file




  ! subroutine hdf5_helpers/open_hdf5_file
  !
  ! Open an existing HDF5 file for both reading and writing.

!  subroutine open_hdf5_file(filename,file_id,error)

!    implicit none
!    character(len=*), intent(in) :: filename      ! name of HDF5 file to open
!    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
!    integer, intent(out) :: error                 ! error flag

!    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)

!  end subroutine open_hdf5_file



  ! subroutine hdf5_helpers/open_hdf5_file_rdonly
  !
  ! Open an existing HDF5 file for reading only.

!  subroutine open_hdf5_file_rdonly(filename,file_id,error)

!     implicit none
!    character(len=*), intent(in) :: filename      ! name of HDF5 file to open
!    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
!    integer, intent(out) :: error                 ! error flag

!    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)

!  end subroutine open_hdf5_file_rdonly



  ! subroutine hdf5_helpers/open_hdf5_file_rdonly
  !
  ! Open an existing HDF5 file for reading only.
  
!   subroutine open_hdf5_file_rdonly(filename,file_id,error)
!
!     implicit none
!     character(len=*), intent(in) :: filename      ! name of HDF5 file to open
! !     integer(hid_t),intent(in)    :: plist_id      ! property list number
!     integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
!     integer, intent(inout) :: error                 ! error flag
!
! !     call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error,access_prp = plist_id)
!     call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)
!
!   end subroutine open_hdf5_file_rdonly
  

  subroutine open_phdf5_file_rdonly(filename,file_id,plist_id,error)

    implicit none
    character(len=*), intent(in) :: filename      ! name of HDF5 file to open
    integer(hid_t),intent(in)    :: plist_id      ! property list number
    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
    integer, intent(inout) :: error                 ! error flag

    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error,access_prp = plist_id)
    !call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)

  end subroutine open_phdf5_file_rdonly
  
  
  subroutine open_phdf5_file_rdwr(filename,file_id,plist_id,error)

    implicit none
    character(len=*), intent(in) :: filename      ! name of HDF5 file to open
    integer(hid_t),intent(in)    :: plist_id      ! property list number
    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
    integer, intent(inout) :: error                 ! error flag

    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error,access_prp = plist_id)

  end subroutine open_phdf5_file_rdwr


  ! subroutine hdf5_helpers/close_hdf5_file
  !
  ! Close a HDF5 file handle.

  subroutine finalise_phdf5_file(file_id,error)

    implicit none
    integer(hid_t), intent(in) :: file_id        ! HDF5 file handle
    integer :: error

     !
     ! Close the file.
     !
     CALL h5fclose_f(file_id, error)
     !
     ! Close FORTRAN interfaces and HDF5 library.
     !
     CALL h5close_f(error)!ok

  end subroutine finalise_phdf5_file




  ! subroutine hdf5_helpers/create_hdf5_group
  !
  ! Create a subgroup in a HDF5 file.

!  subroutine create_hdf5_group(name,file_id)

!    implicit none
!    character(len=*), intent(in) :: name          ! name of group to create
!    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle

!    integer :: error
!    integer(hid_t) :: group_id

!    call h5gcreate_f(file_id,name,group_id,error)
!    call h5gclose_f(group_id,error)

!  end subroutine create_hdf5_group






  ! function hdf5_helpers/exists_in_hdf5_file
  !
  ! Returns true or false if a dataset exists or doesn't exist in a HDF5 file.

!  function exists_in_hdf5_file(name,file_id)

 ! .... !

!  end function exists_in_hdf5_file



  ! subroutine hdf5_helpers/read_hdf5_scalar_integer
  !
  ! Read a single integer (H5T_NATIVE_INTEGER) from a HDF5 file.

  subroutine read_phdf5_scalar_integer(name,a,file_id,error)

    implicit none

    character(len=*), intent(in) :: name          ! name of dataset to read
    integer, intent(out)         :: a             ! variable to be read into
    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle
    integer, intent(out)         :: error         ! error flag

    integer(hsize_t) :: vec_dim(7) = 0
    integer(hid_t) :: dset_id
    vec_dim(1) = 1
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,a,vec_dim,error)
    call h5dclose_f(dset_id,error)

  end subroutine read_phdf5_scalar_integer



  ! subroutine hdf5_helpers/read_hdf5_scalar_double
  !
  ! Read a single double precision (H5T_NATIVE_DOUBLE) value from a HDF5 file.

  subroutine read_phdf5_scalar_double(name,a,file_id,error)

    implicit none

    character(len=*), intent(in) :: name          ! name of dataset to read
    real, intent(out)            :: a             ! variable to be read into
    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle
    integer, intent(out)         :: error         ! error flag

    integer(hsize_t) :: vec_dim(7) = 0
    integer(hid_t) :: dset_id

    vec_dim(1) = 1
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,vec_dim,error)
    call h5dclose_f(dset_id,error)

  end subroutine read_phdf5_scalar_double



  ! subroutine hdf5_helpers/read_hdf5_vector_integer
  !
  ! Read a vector of integer values from a HDF5 file.

!  subroutine read_hdf5_vector_integer(name,n,a,file_id,error)

! .. !

!  end subroutine read_hdf5_vector_integer



  ! subroutine hdf5_helpers/read_hdf5_vector_double
  !
  ! Read a vector of double precision values from a HDF5 file.

  subroutine read_phdf5_vector_double(dsetname,nii,px,offset,rdata_sub,file_id,error)
        
    implicit none

    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: nii  ! subdata size
    real,    intent(out)         :: rdata_sub(nii)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: dataspace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(1) :: count  = (/1/)
    integer(hsize_t) ,  DIMENSION(1) :: stride  = (/1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size =0


    dims(1) = nii*px
    ! xzhao: Here it has been assumed the grid partition is uniform.. To be changed later.
    data_dims = dims
    block_size(1) = nii     !
     ! Open the  dataset.
     !
     call h5dopen_f(file_id,dsetname,dset_id,error)
     !
     ! Get dataset's dataspace identifier  and select subset.
     !
     CALL h5dget_space_f(dset_id, dataspace, error)
     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          offset, count, error, stride, block_size)
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     !
     ! Create memory dataspace.
     !
     CALL h5screate_simple_f(1, block_size, memspace, error)

     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rdata_sub,  dims, error,&
                    mem_space_id =memspace,file_space_id=dataspace,xfer_prp =plist_id)


     CALL h5pclose_f(plist_id, error)!ok
     CALL h5sclose_f(dataspace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
    
  end subroutine read_phdf5_vector_double
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine write_phdf5_vector_double(dsetname,nii,px,offset,subdata,file_id,error)

      implicit none

      character(len=*), intent(in) :: dsetname
      integer, intent(in)          :: px
      integer, intent(in)          :: nii  ! subdata size
      !integer,                     :: ni, nj, nk     ! large_data size
      real, intent(in)             :: subdata(nii)
      integer(hid_t), intent(in)   :: file_id
      integer(hid_t)               :: plist_id

      integer(hid_t)   :: filespace     ! Dataspace identifier in file
      integer(hid_t)   :: memspace      ! Dataspace identifier in memory

      integer :: error
      integer(hsize_t) :: data_dims(7) = 0
      integer(hid_t)   :: dset_id
      integer(hsize_t) :: dims(7) = 0
      integer(hsize_t) ,  DIMENSION(1) :: count  = (/1/)
      integer(hsize_t) ,  DIMENSION(1) :: stride  = (/1/)
      integer(hsize_t) ,  DIMENSION(7) :: offset
      integer(hsize_t) ,  DIMENSION(7) :: block_size=0

      dims(1) = nii*px
      data_dims = dims

      block_size(1) = nii

      !
      ! Create the data space for the  dataset.
      !
      CALL h5screate_simple_f(1, dims, filespace, error)!ok
      CALL h5screate_simple_f(1, block_size, memspace, error)!ok
      !
      ! Create chunked dataset.
      !
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)!ok
      CALL h5pset_chunk_f(plist_id, 1, block_size, error)!ok
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_double, filespace, &
                        dset_id, error, plist_id) !ok
      CALL h5sclose_f(filespace, error)!ok
      !
      ! Select hyperslab in the file.
      !
      CALL h5dget_space_f(dset_id, filespace, error) !ok
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                   stride, block_size) ! ok
      !
      ! Create property list for collective dataset write
      !
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      !
      ! Write the dataset collectively.
      !
      CALL h5dwrite_f(dset_id, H5T_NATIVE_double, subdata, dims, error, &
                       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
      !DEALLOCATE(data)
      !
      ! Close dataspaces.
      !
      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      !
      ! Close the property list.
      !
      CALL h5pclose_f(plist_id, error)

  !    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
  !    call h5sclose_f(dspace_id,error)
  !    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
  !    call h5dclose_f(dset_id,error)

    end subroutine write_phdf5_vector_double


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     subroutine rewrite_phdf5_vector_double(dsetname,nii,px,offset,subdata,file_id,error)

      implicit none

      character(len=*), intent(in) :: dsetname
      integer, intent(in)          :: px
      integer, intent(in)          :: nii  ! subdata size
      !integer,                     :: ni, nj, nk     ! large_data size
      real, intent(in)             :: subdata(nii)
      integer(hid_t), intent(in)   :: file_id
      integer(hid_t)               :: plist_id

      integer(hid_t)   :: dataspace     ! Dataspace identifier in file
      integer(hid_t)   :: memspace      ! Dataspace identifier in memory

      integer :: error
      integer(hsize_t) :: data_dims(7) = 0
      integer(hid_t)   :: dset_id
      integer(hsize_t) :: dims(7) = 0
      integer(hsize_t) ,  DIMENSION(1) :: count  = (/1/)
      integer(hsize_t) ,  DIMENSION(1) :: stride  = (/1/)
      integer(hsize_t) ,  DIMENSION(7) :: offset
      integer(hsize_t) ,  DIMENSION(7) :: block_size=0

      dims(1) = nii*px
      data_dims = dims

      block_size(1) = nii    
    
      call h5dopen_f(file_id,dsetname,dset_id,error)
      !
      ! Get dataset's dataspace identifier  and select subset.
      !
      CALL h5dget_space_f(dset_id, dataspace, error)
      CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
           offset, count, error, stride, block_size)
      !
      ! Create property list for collective dataset write
      !
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      !
      ! Create memory dataspace.
      !
      CALL h5screate_simple_f(1, block_size, memspace, error)
    
      !
      CALL h5dwrite_f(dset_id, H5T_NATIVE_double, subdata, dims, error, &
                       file_space_id = dataspace, mem_space_id = memspace, xfer_prp = plist_id)
      !   CALL h5pclose_f(plist_id, error)!ok
     CALL h5sclose_f(dataspace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)

    end subroutine rewrite_phdf5_vector_double


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! subroutine hdf5_helpers/read_hdf5_string
  !
  ! Read a string from a HDF5 file.write_hdf5_file

!  subroutine read_hdf5_string(name,n,a,file_id,error)

! .. !

!  end subroutine read_hdf5_string



  ! subroutine hdf5_helpers/read_hdf5_2darray_integer
  !
  ! Read a two-dimensional array of integers from a HDF5 file
!  subroutine read_hdf5_2darray_integer(name,ni,nj,a,file_id,error)

! .. !

!  end subroutine read_hdf5_2darray_integer

!  subroutine read_hdf5_3darray_integer(name,ni,nj,nk,a,file_id,error)

! .. !

!  end subroutine read_hdf5_3darray_integer


  ! subroutine hdf5_helpers/read_hdf5_2darray_double
  !
  ! Read a two-dimensional array of double precision values from a HDF5 file.

!  subroutine read_hdf5_2darray_double(name,ni,nj,a,file_id,error)

! .. !

!  end subroutine read_hdf5_2darray_double

subroutine read_phdf5_2darray_double(dsetname,nii,njj,px,py,offset,rdata_sub,file_id,error)

  implicit none

  character(len=*), intent(in) :: dsetname
  integer, intent(in)          :: px
  integer, intent(in)          :: py
  integer, intent(in)          :: nii, njj  ! subdata size
  real,    intent(out)         :: rdata_sub(nii,njj)
  integer(hid_t), intent(in)   :: file_id
  integer(hid_t)               :: plist_id

  integer(hid_t)   :: dataspace     ! Dataspace identifier in file
  integer(hid_t)   :: memspace      ! Dataspace identifier in memory

  integer :: error
  integer(hsize_t) :: data_dims(7) = 0
  integer(hid_t)   :: dset_id
  integer(hsize_t) :: dims(7) = 0
  integer(hsize_t) ,  DIMENSION(2) :: count  = (/1,1/)
  integer(hsize_t) ,  DIMENSION(2) :: stride  = (/1,1/)
  integer(hsize_t) ,  DIMENSION(7) :: offset
  integer(hsize_t) ,  DIMENSION(7) :: block_size =0


  dims(1) = nii*px
  dims(2) = njj*py
  ! xzhao: Here it has been assumed the grid partition is uniform.. To be changed later.
  data_dims = dims
  block_size(1) = nii
  block_size(2) = njj
   !
   ! Open the  dataset.
   !
   call h5dopen_f(file_id,dsetname,dset_id,error)
   !
   ! Get dataset's dataspace identifier  and select subset.
   !
   CALL h5dget_space_f(dset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
        offset, count, error, stride, block_size)
   !
   ! Create property list for collective dataset write
   !
   CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
   CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
   !
   ! Create memory dataspace.
   !
   CALL h5screate_simple_f(2, block_size, memspace, error)

   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rdata_sub,  dims, error,&
                  mem_space_id =memspace,file_space_id=dataspace,xfer_prp =plist_id)


   CALL h5pclose_f(plist_id, error)!ok
   CALL h5sclose_f(dataspace, error)
   CALL h5sclose_f(memspace, error)
   CALL h5dclose_f(dset_id, error)

end subroutine read_phdf5_2darray_double


  ! subroutine hdf5_helpers/read_hdf5_3darray_double
  !
  ! Read a three-dimensional array of double precision values from a HDF5 file.

  subroutine read_phdf5_3darray_double(dsetname,nii,njj,nkk,px,py,pz,offset,rdata_sub,file_id,error)

    implicit none

    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: py
    integer, intent(in)          :: pz
    integer, intent(in)          :: nii, njj, nkk  ! subdata size
    real,    intent(out)         :: rdata_sub(nii,njj,nkk)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: dataspace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(3) :: count  = (/1,1,1/)
    integer(hsize_t) ,  DIMENSION(3) :: stride  = (/1,1,1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size =0


    dims(1) = nii*px
    dims(2) = njj*py
    dims(3) = nkk*pz
    ! xzhao: Here it has been assumed the grid partition is uniform.. To be changed later.
    data_dims = dims
    block_size(1) = nii
    block_size(2) = njj
    block_size(3) = nkk
     !
     ! Open the  dataset.
     !
     call h5dopen_f(file_id,dsetname,dset_id,error)
     !
     ! Get dataset's dataspace identifier  and select subset.
     !
     CALL h5dget_space_f(dset_id, dataspace, error)
     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          offset, count, error, stride, block_size)
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     !
     ! Create memory dataspace.
     !
     CALL h5screate_simple_f(3, block_size, memspace, error)

     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rdata_sub,  dims, error,&
                    mem_space_id =memspace,file_space_id=dataspace,xfer_prp =plist_id)


     CALL h5pclose_f(plist_id, error)!ok
     CALL h5sclose_f(dataspace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)

  end subroutine read_phdf5_3darray_double



! OLD::: subroutine read_phdf5_4darray_double(name,ni,nj,nk,nl,a,file_id,error)
 subroutine read_phdf5_4darray_double(dsetname,nii,njj,nkk,nll,px,py,pz,offset,rdata_sub,file_id,error)

    implicit none


    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: py
    integer, intent(in)          :: pz
    integer, intent(in)          :: nii, njj, nkk, nll  ! subdata size
    real,    intent(out)          :: rdata_sub(nii,njj,nkk,nll)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: dataspace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(4) :: count  = (/1,1,1,1/)
    integer(hsize_t) ,  DIMENSION(4) :: stride  = (/1,1,1,1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size =0

    dims(1) = nii*px
    dims(2) = njj*py
    dims(3) = nkk*pz
    dims(4) = nll*1
    data_dims = dims
    block_size(1) = nii
    block_size(2) = njj
    block_size(3) = nkk
    block_size(4) = nll
     !
     ! Open the  dataset.
     !
     call h5dopen_f(file_id,dsetname,dset_id,error)
!      write(*,*) 'error1',error
     !
     ! Get dataset's dataspace identifier  and select subset.
     !
     CALL h5dget_space_f(dset_id, dataspace, error)
!      write(*,*) 'error2',error
     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
          offset, count, error, stride, block_size)
!                write(*,*) 'error3',error

     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!           write(*,*) 'error4',error
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!           write(*,*) 'error5',error
     !
     ! Create memory dataspace.
     !
     CALL h5screate_simple_f(4, block_size, memspace, error)
!           write(*,*) 'error6',error

     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rdata_sub,  dims, error,&
                   mem_space_id =memspace ,file_space_id=dataspace ,xfer_prp =plist_id)
!                         write(*,*) 'error7',error

     CALL h5pclose_f(plist_id, error)!ok
     CALL h5sclose_f(dataspace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
!           write(*,*) 'error8',error

  end subroutine read_phdf5_4darray_double





!  subroutine read_hdf5_5darray_double(name,ni,nj,nk,nl,ns,a,file_id,error)

!    implicit none
!    character(len=*), intent(in) :: name
!    integer, intent(in)          :: ni, nj, nk, nl, ns
!    real, intent(out)   :: a(ni,nj,nk, nl,ns)
!    integer(hid_t), intent(in)   :: file_id
!    integer, intent(out)         :: error

!    integer(hsize_t) :: data_dims(5) = 0
!    integer(hid_t) :: dset_id

!    data_dims(1) = ni
!    data_dims(2) = nj
!    data_dims(3) = nk
!    data_dims(4) = nl
!    data_dims(5) = ns
!    call h5dopen_f(file_id,name,dset_id,error)
!    ! call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
!    call h5dclose_f(dset_id,error)

!  end subroutine read_hdf5_5darray_double

!!****f* hdf5_helpers/write_hdf5_scalar_integer
!!
!! NAME
!!   write_hdf5_scalar_integer - write a single integer value
!!
!! SYNOPSIS
!!   call write_hdf5_scalar_integer(name,a,file_id)
!!
!! DESCRIPTION
!!   Write a single integer (H5T_NATIVE_INTEGER) value to a HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier character string of the scalar [in]
!!   a                integer scalar [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

  subroutine write_phdf5_scalar_integer(name,a,file_id,error)

     IMPLICIT NONE

    character(len=*), intent(in) :: name
    integer, intent(in)          :: a
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t)             :: data_dims(7) = 0
    integer(hid_t)               :: dspace_id, dset_id
    integer(hsize_t)             :: dims(7) = 0

    data_dims(1) = 1
    dims(1) = 1

    call h5screate_simple_f(1,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)

  end subroutine write_phdf5_scalar_integer
  
  subroutine rewrite_phdf5_scalar_integer(name,a,file_id,error)

     IMPLICIT NONE

    character(len=*), intent(in) :: name
    integer, intent(in)          :: a
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t)             :: data_dims(7) = 0
    integer(hid_t)               :: dspace_id, dset_id
    integer(hsize_t)             :: dims(7) = 0

    data_dims(1) = 1
    dims(1) = 1
    
    call h5dopen_f(file_id,name,dset_id,error)  
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)

  end subroutine rewrite_phdf5_scalar_integer
  
  


!!****f* hdf5_helpers/write_hdf5_scalar_double
!!
!! NAME
!!   write_hdf5_scalar_double - write a single double precision value
!!
!! SYNOPSIS
!!   call write_hdf5_scalar_double(name,a,file_id)
!!
!! DESCRIPTION
!!   Write a single double precision (H5T_NATIVE_DOUBLE) value to a HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier character string of the scalar [in]
!!   a                double precision scalar [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

  subroutine write_phdf5_scalar_double(name,a,file_id,error)

    implicit none

    character(len=*), intent(in) :: name
    real, intent(in)     :: a
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0
    data_dims(1) = 1
    dims(1) = 1

    call h5screate_simple_f(1,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_double,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_double,a,data_dims,error)
    call h5dclose_f(dset_id,error)

  end subroutine write_phdf5_scalar_double
  
  subroutine rewrite_phdf5_scalar_double(name,a,file_id,error)

     IMPLICIT NONE

    character(len=*), intent(in) :: name
    real, intent(in)          :: a
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t)             :: data_dims(7) = 0
    integer(hid_t)               :: dspace_id, dset_id
    integer(hsize_t)             :: dims(7) = 0

    data_dims(1) = 1
    dims(1) = 1
    
    call h5dopen_f(file_id,name,dset_id,error)  
    call h5dwrite_f(dset_id,H5T_NATIVE_double,a,data_dims,error)
    call h5dclose_f(dset_id,error)

  end subroutine rewrite_phdf5_scalar_double

!!****f* hdf5_helpers/write_hdf5_vector_integer
!!
!! NAME
!!   write_hdf5_vector_integer - write a vector of integer values
!!
!! SYNOPSIS
!!   call write_hdf5_vector_integer(name,n,a,file_id)
!!
!! DESCRIPTION
!!   Write a vector of integer (H5T_NATIVE_INTEGER) values to a HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier character string of the vector [in]
!!   n                length of the vector [in]
!!   a                integer vector of length n [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

!  subroutine write_hdf5_vector_integer(name,n,a,file_id)

! .. !

!  end subroutine write_hdf5_vector_integer

!!****f* hdf5_helpers/write_hdf5_vector_double
!!
!! NAME
!!   write_hdf5_vector_double - write a vector of double precision values
!!
!! SYNOPSIS
!!   call write_hdf5_vector_double(name,n,a,file_id)
!!
!! DESCRIPTION
!!   Write a vector of double precision (H5T_NATIVE_DOUBLE) values to a HDF5
!!   file.
!!
!! ARGUMENTS
!!   name             identifier character string of the vector [in]
!!   n                length of the vector [in]
!!   a                double precision vector of length n [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

!  subroutine write_hdf5_vector_double(name,n,a,file_id)

! .. !

!  end subroutine write_hdf5_vector_double

!!****f* hdf5_helpers/write_hdf5_string
!!
!! NAME
!!   write_hdf5_string - write a character string
!!
!! SYNOPSIS
!!   call write_hdf5_string(name,n,a,file_id)
!!
!! DESCRIPTION
!!   Write a vector of characters (a string) to a HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier character string of the vector [in]
!!   n                length of the string [in]
!!   a                string [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

  !subroutine write_hdf5_string(name,n,a,file_id)

  !..!

  !end subroutine write_hdf5_string

!!****f* hdf5_helpers/write_hdf5_2darray_integer
!!
!! NAME
!!   write_hdf5_2darray_integer - write a 2d array of integers
!!
!! SYNOPSIS
!!   call write_hdf5_2darray_integer(name,ni,nj,a,file_id)
!!
!! DESCRIPTION
!!   Write a two-dimensional array of integer values to a HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier of the 2D array [in]
!!   ni,nj            size of the 2D array [in]
!!   a                2D array of integer values [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***



!  subroutine write_hdf5_2darray_integer(name,ni,nj,a,file_id)

!..!

!  end subroutine write_hdf5_2darray_integer




!!****f* hdf5_helpers/write_hdf5_2darray_double
!!
!! NAME
!!   write_hdf5_2darray_double - write a 2d array of double precision values
!!
!! SYNOPSIS
!!   call write_hdf5_2darray_double (name,ni,nj,a,file_id)
!!
!! DESCRIPTION
!!   Write a two-dimensional array of double precision values to a HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier of the 2D array [in]
!!   ni,nj            size of the 2D array [in]
!!   a                2D array of double precision values [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

!  subroutine write_hdf5_2darray_double(name,ni,nj,a,file_id)

  !..!

!  end subroutine write_hdf5_2darray_double

!!****f* hdf5_helpers/write_hdf5_3darray_double
!!
!! NAME
!!   write_hdf5_3darray_double - write a 3d array of double precision values
!!
!! SYNOPSIS
!!   call write_hdf5_3darray_double(name,ni,nj,nk,a,file_id)
!!
!! DESCRIPTION
!!   Write a three-dimensional array of double precision values to a
!!   HDF5 file.
!!
!! ARGUMENTS
!!   name             identifier of the 3D array [in]
!!   ni,nj,nk         size of the 3D array [in]
!!   a                3D array of double precision values [in]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_phdf5_2darray_double(dsetname,nii,njj,px,py,offset,subdata,file_id,error)

    implicit none

    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: py
    integer, intent(in)          :: nii, njj  ! subdata size
    !integer,                     :: ni, nj, nk     ! large_data size
    real, intent(in)             :: subdata(nii,njj)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: filespace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(2) :: count  = (/1,1/)
    integer(hsize_t) ,  DIMENSION(2) :: stride  = (/1,1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size=0

    dims(1) = nii*px
    dims(2) = njj*py
    data_dims = dims

    block_size(1) = nii
    block_size(2) = njj
    
!     write(*,*), dims,block_size

    !
    ! Create the data space for the  dataset.
    !
    CALL h5screate_simple_f(2, dims, filespace, error)!ok
    CALL h5screate_simple_f(2, block_size, memspace, error)!ok
    !
    ! Create chunked dataset.
    !
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)!ok
    CALL h5pset_chunk_f(plist_id, 2, block_size, error)!ok
    CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_double, filespace, &
                      dset_id, error, plist_id) !ok
    CALL h5sclose_f(filespace, error)!ok
    !
    ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, error) !ok
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block_size) ! ok
    !
    ! Create property list for collective dataset write
    !
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    !
    ! Write the dataset collectively.
    !
    CALL h5dwrite_f(dset_id, H5T_NATIVE_double, subdata, dims, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    !DEALLOCATE(data)
    !
    ! Close dataspaces.
    !
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    !
    ! Close the dataset.
    !
    CALL h5dclose_f(dset_id, error)
    !
    ! Close the property list.
    !
    CALL h5pclose_f(plist_id, error)

!    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
!    call h5sclose_f(dspace_id,error)
!    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
!    call h5dclose_f(dset_id,error)

  end subroutine write_phdf5_2darray_double


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine rewrite_phdf5_2darray_double(dsetname,nii,njj,px,py,offset,subdata,file_id,error)

    implicit none

    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: py
    integer, intent(in)          :: nii, njj  ! subdata size
    !integer,                     :: ni, nj, nk     ! large_data size
    real, intent(in)             :: subdata(nii,njj)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: dataspace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(2) :: count  = (/1,1/)
    integer(hsize_t) ,  DIMENSION(2) :: stride  = (/1,1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size=0

    dims(1) = nii*px
    dims(2) = njj*py
    data_dims = dims

    block_size(1) = nii
    block_size(2) = njj
    
    
    call h5dopen_f(file_id,dsetname,dset_id,error)
    !
    ! Get dataset's dataspace identifier  and select subset.
    !
    CALL h5dget_space_f(dset_id, dataspace, error)
    CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
         offset, count, error, stride, block_size)
    !
    ! Create property list for collective dataset write
    !
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !
    ! Create memory dataspace.
    !
    CALL h5screate_simple_f(2, block_size, memspace, error)
    
    !
    CALL h5dwrite_f(dset_id, H5T_NATIVE_double, subdata, dims, error, &
                     file_space_id = dataspace, mem_space_id = memspace, xfer_prp = plist_id)
    !   CALL h5pclose_f(plist_id, error)!ok
   CALL h5sclose_f(dataspace, error)
   CALL h5sclose_f(memspace, error)
   CALL h5dclose_f(dset_id, error)

  end subroutine rewrite_phdf5_2darray_double


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_phdf5_3darray_double(dsetname,nii,njj,nkk,px,py,pz,offset,subdata,file_id,error)

    implicit none

    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: py
    integer, intent(in)          :: pz
    integer, intent(in)          :: nii, njj, nkk  ! subdata size
    !integer,                     :: ni, nj, nk     ! large_data size
    real, intent(in)             :: subdata(nii,njj,nkk)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: filespace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(3) :: count  = (/1,1,1/)
    integer(hsize_t) ,  DIMENSION(3) :: stride  = (/1,1,1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size=0

    dims(1) = nii*px
    dims(2) = njj*py
    dims(3) = nkk*pz
    data_dims = dims

    block_size(1) = nii
    block_size(2) = njj
    block_size(3) = nkk

    !
    ! Create the data space for the  dataset.
    !
    CALL h5screate_simple_f(3, dims, filespace, error)!ok
    CALL h5screate_simple_f(3, block_size, memspace, error)!ok
    !
    ! Create chunked dataset.
    !
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)!ok
    CALL h5pset_chunk_f(plist_id, 3, block_size, error)!ok
    CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_double, filespace, &
                      dset_id, error, plist_id) !ok
    CALL h5sclose_f(filespace, error)!ok
    !
    ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, error) !ok
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block_size) ! ok
    !
    ! Create property list for collective dataset write
    !
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    !
    ! Write the dataset collectively.
    !
    CALL h5dwrite_f(dset_id, H5T_NATIVE_double, subdata, dims, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
    !DEALLOCATE(data)
    !
    ! Close dataspaces.
    !
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    !
    ! Close the dataset.
    !
    CALL h5dclose_f(dset_id, error)
    !
    ! Close the property list.
    !
    CALL h5pclose_f(plist_id, error)

!    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
!    call h5sclose_f(dspace_id,error)
!    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
!    call h5dclose_f(dset_id,error)

  end subroutine write_phdf5_3darray_double


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine write_phdf5_4darray_double(dsetname,nii,njj,nkk,nll,px,py,pz,offset,subdata,file_id,error)

    implicit none

    character(len=*), intent(in) :: dsetname
    integer, intent(in)          :: px
    integer, intent(in)          :: py
    integer, intent(in)          :: pz
    integer, intent(in)          :: nii, njj, nkk, nll  ! subdata size
    real,    intent(in)          :: subdata(nii,njj,nkk,nll)
    integer(hid_t), intent(in)   :: file_id
    integer(hid_t)               :: plist_id

    integer(hid_t)   :: filespace     ! Dataspace identifier in file
    integer(hid_t)   :: memspace      ! Dataspace identifier in memory

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t)   :: dset_id
    integer(hsize_t) :: dims(7) = 0
    integer(hsize_t) ,  DIMENSION(4) :: count  = (/1,1,1,1/)
    integer(hsize_t) ,  DIMENSION(4) :: stride  = (/1,1,1,1/)
    integer(hsize_t) ,  DIMENSION(7) :: offset
    integer(hsize_t) ,  DIMENSION(7) :: block_size =0

    dims(1) = nii*px
    dims(2) = njj*py
    dims(3) = nkk*pz
    dims(4) = nll*1
    data_dims = dims
    block_size(1) = nii
    block_size(2) = njj
    block_size(3) = nkk
    block_size(4) = nll
    !
    ! Create the data space for the  dataset.
    !
    CALL h5screate_simple_f(4, dims, filespace, error)!ok
    CALL h5screate_simple_f(4, block_size, memspace, error)!ok
    !
    ! Create chunked dataset.
    !
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)!ok
    CALL h5pset_chunk_f(plist_id, 4, block_size, error)!ok
    CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_double, filespace, &
                      dset_id, error, plist_id) !ok
    CALL h5sclose_f(filespace, error)!ok
    !
    ! Select hyperslab in the file.
    !
    CALL h5dget_space_f(dset_id, filespace, error) !ok
    CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block_size) ! ok
    !
    ! Create property list for collective dataset write
    !
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    !
    ! Write the dataset collectively.
    !
    CALL h5dwrite_f(dset_id, H5T_NATIVE_double, subdata, dims, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    !DEALLOCATE(subdata)

    !
    ! Close dataspaces.
    !
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    !
    ! Close the dataset.
    !
    CALL h5dclose_f(dset_id, error)
    !
    ! Close the property list.
    !
    CALL h5pclose_f(plist_id, error)

  end subroutine write_phdf5_4darray_double







!  subroutine write_hdf5_5darray_double(name,ni,nj,nk,nl,nm,a,file_id)

  !..!

!  end subroutine write_hdf5_5darray_double

!!****f* hdf5_helpers/copy_hdf5_scalar_integer
!!
!! NAME
!!   copy_hdf5_scalar_integer - copy an integer value from one file to another
!!
!! SYNOPSIS
!!   call copy_hdf5_scalar_integer(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a single integer value from one HDF5 file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

!  subroutine copy_hdf5_scalar_integer(name,file_id,new_name,ofile_id)

  !..!

!  end subroutine copy_hdf5_scalar_integer

!!****f* hdf5_helpers/copy_hdf5_scalar_double
!!
!! NAME
!!   copy_hdf5_scalar_double - copy a double value from one file to another
!!
!! SYNOPSIS
!!   call copy_hdf5_scalar_double(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a single double precision value from one HDF5 file
!!   to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

  !subroutine copy_hdf5_scalar_double(name,file_id,new_name,ofile_id)

  !..!

  !end subroutine copy_hdf5_scalar_double

!!****f* hdf5_helpers/copy_hdf5_vector_integer
!!
!! NAME
!!   copy_hdf5_vector_integer - copy a integer vector from one file to another
!!
!! SYNOPSIS
!!   call copy_hdf5_vector_integer(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a vector of integers from one HDF5 file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

  !subroutine copy_hdf5_vector_integer(name,file_id,new_name,ofile_id)

  !..!

  !end subroutine copy_hdf5_vector_integer

!!****f* hdf5_helpers/copy_hdf5_vector_double
!!
!! NAME
!!   copy_hdf5_vector_double - copy a double vector from one file to another
!!
!! SYNOPSIS
!!   call copy_hdf5_vector_double(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a vector of double precision values from one HDF5
!!   file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

!  subroutine copy_hdf5_vector_double(name,file_id,new_name,ofile_id)

    !..!

!  end subroutine copy_hdf5_vector_double

!!****f* hdf5_helpers/copy_hdf5_string
!!
!! NAME
!!   copy_hdf5_string - copy a string from one file to another
!!
!! SYNOPSIS
!!   call copy_hdf5_string(name,file_id,new_name,ofile_id,n)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a string from one HDF5 file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!   n                string length [in]
!!
!!***

!  subroutine copy_hdf5_string(name,file_id,new_name,ofile_id,n)

   !..!

!  end subroutine copy_hdf5_string

!!****f* hdf5_helpers/copy_hdf5_2darray_integer
!!
!! NAME
!!   copy_hdf5_2darray_integer - copy a integer 2D array
!!
!! SYNOPSIS
!!   call copy_hdf5_2darray_integer(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a two-dimensional array of integer values
!!   from one HDF5 file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

!  subroutine copy_hdf5_2darray_integer(name,file_id,new_name,ofile_id)

  !..!

!  end subroutine copy_hdf5_2darray_integer

!!****f* hdf5_helpers/copy_hdf5_2darray_double
!!
!! NAME
!!   copy_hdf5_2darray_double - copy a double 2D array
!!
!! SYNOPSIS
!!   call copy_hdf5_2darray_double(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a two-dimensional array of double values
!!   from one HDF5 file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

!  subroutine copy_hdf5_2darray_double(name,file_id,new_name,ofile_id)

    !..!

!  end subroutine copy_hdf5_2darray_double

!!****f* hdf5_helpers/copy_hdf5_3darray_double
!!
!! NAME
!!   copy_hdf5_3darray_double - copy a double 3D array from one file to another
!!
!! SYNOPSIS
!!   call copy_hdf5_3darray_double(name,file_id,new_name,ofile_id)
!!
!! DESCRIPTION
!!   Copy (and/or rename) a three-dimensional array of double precision values
!!   from one HDF5 file to another.
!!
!! ARGUMENTS
!!   name             identifier of the source value [in]
!!   file_id          HDF5 identifier for the file handle of the source [in]
!!   new_name         identifier of the target value [in]
!!   ofile_id         HDF5 identifier for the file handle of the target [in]
!!
!!***

!  subroutine copy_hdf5_3darray_double(name,file_id,new_name,ofile_id)

  !..!

!  end subroutine copy_hdf5_3darray_double

!!****f* hdf5_helpers/getdims_hdf5_vector
!!
!! NAME
!!   getdims_hdf5_vector - get the dimensions (length) of a vector
!!
!! SYNOPSIS
!!   call getdims_hdf5_vector(name,n,file_id)
!!
!! DESCRIPTION
!!   Get the dimensions (length) of a vector in a HDF5 file.
!!
!! ARGUMENTS
!!   name             array identifier [in]
!!   n                dimension array [out]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

  !subroutine getdims_hdf5_vector(name,n,file_id)

  !..!

  !end subroutine getdims_hdf5_vector

!!****f* hdf5_helpers/getdims_hdf5_2darray
!!
!! NAME
!!   getdims_hdf5_2darray - get the dimensions of a 2D array
!!
!! SYNOPSIS
!!   call getdims_hdf5_2darray(name,n,file_id)
!!
!! DESCRIPTION
!!   Get the dimensions of a 2D array in a HDF5 file.
!!
!! ARGUMENTS
!!   name             array identifier [in]
!!   n                dimension array [out]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

!  subroutine getdims_hdf5_2darray(name,n,file_id)

  !..!

!  end subroutine getdims_hdf5_2darray

!!****f* hdf5_helpers/getdims_hdf5_3darray
!!
!! NAME
!!   getdims_hdf5_3darray - get the dimensions of a 3D array
!!
!! SYNOPSIS
!!   call getdims_hdf5_3darray(name,n,file_id)
!!
!! DESCRIPTION
!!   Get the dimensions of a 3D array in a HDF5 file.
!!
!! ARGUMENTS
!!   name             array identifier [in]
!!   n                dimension array [out]
!!   file_id          HDF5 identifier for the file handle [in]
!!
!!***

!  subroutine getdims_hdf5_3darray(name,n,file_id)

    !..!

!  end subroutine getdims_hdf5_3darray

end module phdf5_helpers

! #endif
