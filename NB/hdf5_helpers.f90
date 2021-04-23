! module hdf5_helpers
!
! Convenience subroutines for working with HDF5 files.

module hdf5_helpers

  use hdf5

contains

  ! subroutine hdf5_helpers/initialise_hdf5
  !
  ! Startup HDF5 by running h5open_f. Required before calling any
  ! HDF5 routines.

  subroutine initialise_hdf5(error)

    implicit none
    integer, intent(out) :: error                 ! error flag

    call h5open_f(error)

  end subroutine initialise_hdf5

  ! subroutine hdf5_helpers/finalise_hdf5
  !
  ! Finish with HDF5 by running h5close_f. 

  subroutine finalise_hdf5(error)

    implicit none
    integer, intent(out) :: error                 ! error flag

    call h5close_f(error)

  end subroutine finalise_hdf5

  ! subroutine hdf5_helpers/create_hdf5_file
  !
  ! Create and open a new HDF5 file.

  subroutine create_hdf5_file(filename,file_id)

    implicit none
    character(len=*), intent(in) :: filename      ! name of HDF5 file to create
    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle

    integer :: error
    
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,error)
    
  end subroutine create_hdf5_file

  ! subroutine hdf5_helpers/open_hdf5_file
  !
  ! Open an existing HDF5 file for both reading and writing.

  subroutine open_hdf5_file(filename,file_id,error)

    implicit none
    character(len=*), intent(in) :: filename      ! name of HDF5 file to open
    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
    integer, intent(out) :: error                 ! error flag
    
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)

  end subroutine open_hdf5_file

  ! subroutine hdf5_helpers/open_hdf5_file_rdonly
  !
  ! Open an existing HDF5 file for reading only.

  subroutine open_hdf5_file_rdonly(filename,file_id,error)

    implicit none
    character(len=*), intent(in) :: filename      ! name of HDF5 file to open
    integer(hid_t), intent(out)  :: file_id       ! HDF5 file handle
    integer, intent(out) :: error                 ! error flag
    
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,error)

  end subroutine open_hdf5_file_rdonly

  ! subroutine hdf5_helpers/close_hdf5_file
  !
  ! Close a HDF5 file handle.

  subroutine close_hdf5_file(file_id)

    implicit none
    integer(hid_t), intent(in) :: file_id         ! HDF5 file handle

    integer :: error
    
    call h5fclose_f(file_id,error)
    
  end subroutine close_hdf5_file

  ! subroutine hdf5_helpers/create_hdf5_group
  !
  ! Create a subgroup in a HDF5 file.

  subroutine create_hdf5_group(name,file_id)

    implicit none
    character(len=*), intent(in) :: name          ! name of group to create
    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle

    integer :: error
    integer(hid_t) :: group_id

    call h5gcreate_f(file_id,name,group_id,error)
    call h5gclose_f(group_id,error)

  end subroutine create_hdf5_group

  ! function hdf5_helpers/exists_in_hdf5_file
  !
  ! Returns true or false if a dataset exists or doesn't exist in a HDF5 file.

  function exists_in_hdf5_file(name,file_id)

    implicit none
    character(len=*), intent(in) :: name         ! name of dataset to inquire
    integer(hid_t), intent(in)   :: file_id      ! HDF5 file handle
    logical :: exists_in_hdf5_file

    integer :: error
    integer(hid_t) :: dataset

    ! Turn off automatic error control
    
    call h5eset_auto_f(0,error)

    ! Try to open the dataset

    call h5dopen_f(file_id,name,dataset,error)
    if (error < 0) then
       exists_in_hdf5_file = .false.
    else
       call h5dclose_f(dataset,error)
       exists_in_hdf5_file = .true.
    endif

    ! Turn automatic error control back on

    call h5eset_auto_f(1,error)
      
  end function exists_in_hdf5_file

  ! subroutine hdf5_helpers/read_hdf5_scalar_integer
  !
  ! Read a single integer (H5T_NATIVE_INTEGER) from a HDF5 file.

  subroutine read_hdf5_scalar_integer(name,a,file_id,error)

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
    
  end subroutine read_hdf5_scalar_integer

  ! subroutine hdf5_helpers/read_hdf5_scalar_double
  !
  ! Read a single double precision (H5T_NATIVE_DOUBLE) value from a HDF5 file.

  subroutine read_hdf5_scalar_double(name,a,file_id,error)

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
    
  end subroutine read_hdf5_scalar_double

  ! subroutine hdf5_helpers/read_hdf5_vector_integer
  !
  ! Read a vector of integer values from a HDF5 file.
  
  subroutine read_hdf5_vector_integer(name,n,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name          ! name of dataset to read
    integer, intent(in)          :: n             ! length of vector
    integer, intent(out)         :: a(n)          ! vector to be read into
    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle
    integer, intent(out)         :: error         ! error flag

    integer(hsize_t) :: vec_dim(7) = 0
    integer(hid_t) :: dset_id

    vec_dim(1) = n
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,a,vec_dim,error)
    call h5dclose_f(dset_id,error)

  end subroutine read_hdf5_vector_integer

  ! subroutine hdf5_helpers/read_hdf5_vector_double
  !
  ! Read a vector of double precision values from a HDF5 file.

  subroutine read_hdf5_vector_double(name,n,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name          ! name of dataset to read
    integer, intent(in)          :: n             ! length of vector
    real, intent(out)            :: a(n)          ! vector to be read into
    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle
    integer, intent(out)         :: error         ! error flag

    integer(hsize_t) :: vec_dim(7) = 0
    integer(hid_t) :: dset_id
    
    vec_dim(1) = n
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,vec_dim,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_vector_double

  ! subroutine hdf5_helpers/read_hdf5_string
  !
  ! Read a string from a HDF5 file.write_hdf5_file

  subroutine read_hdf5_string(name,n,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name          ! name of dataset to read
    integer, intent(in)          :: n             ! length of string
    character(len=n), intent(out):: a             ! string to be read into
    integer(hid_t), intent(in)   :: file_id       ! HDF5 file handle
    integer, intent(out)         :: error         ! error flag

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id, type_id
    integer(hsize_t) :: dims(7) = 0
    integer(size_t) :: typesize

    dims(1) = 1
    data_dims = dims
    
    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,error)
    typesize = n
    call h5tset_size_f(type_id,typesize,error)
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,type_id,a,data_dims,error)
    call h5dclose_f(dset_id,error)
       
  end subroutine read_hdf5_string

  ! subroutine hdf5_helpers/read_hdf5_2darray_integer
  !
  ! Read a two-dimensional array of integers from a HDF5 file
  subroutine read_hdf5_2darray_integer(name,ni,nj,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj
    integer, intent(out)         :: a(ni,nj)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out)         :: error

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id

    data_dims(1) = ni
    data_dims(2) = nj
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_2darray_integer
  
  subroutine read_hdf5_3darray_integer(name,ni,nj,nk,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk
    integer, intent(out)         :: a(ni,nj,nk)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out)         :: error

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id

    data_dims(1) = ni
    data_dims(2) = nj
    data_dims(3) = nk
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_3darray_integer

  ! subroutine hdf5_helpers/read_hdf5_2darray_double
  !
  ! Read a two-dimensional array of double precision values from a HDF5 file.

  subroutine read_hdf5_2darray_double(name,ni,nj,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj
    real, intent(out)   :: a(ni,nj)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out)         :: error

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id

    data_dims(1) = ni
    data_dims(2) = nj
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_2darray_double

  ! subroutine hdf5_helpers/read_hdf5_3darray_double
  !
  ! Read a three-dimensional array of double precision values from a HDF5 file.

  subroutine read_hdf5_3darray_double(name,ni,nj,nk,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk
    real, intent(out)   :: a(ni,nj,nk)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out)         :: error

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id

    data_dims(1) = ni
    data_dims(2) = nj
    data_dims(3) = nk
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_3darray_double
  
  subroutine read_hdf5_4darray_double(name,ni,nj,nk,nl,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk, nl
    real, intent(out)            :: a(ni,nj,nk,nl)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out)         :: error

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id

    data_dims(1) = ni
    data_dims(2) = nj
    data_dims(3) = nk
    data_dims(4) = nl
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_4darray_double
  
  subroutine read_hdf5_5darray_double(name,ni,nj,nk,nl,ns,a,file_id,error)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk, nl, ns
    real, intent(out)   :: a(ni,nj,nk, nl,ns)
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out)         :: error

    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dset_id

    data_dims(1) = ni
    data_dims(2) = nj
    data_dims(3) = nk
    data_dims(4) = nl
    data_dims(5) = ns
    call h5dopen_f(file_id,name,dset_id,error)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine read_hdf5_5darray_double

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

  subroutine write_hdf5_scalar_integer(name,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: a
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0

    data_dims(1) = 1
    dims(1) = 1

    call h5screate_simple_f(1,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_scalar_integer

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

  subroutine write_hdf5_scalar_double(name,a,file_id)

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
    
  end subroutine write_hdf5_scalar_double

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

  subroutine write_hdf5_vector_integer(name,n,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: n
    integer, intent(in)          :: a(n)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0

    dims(1) = n
    data_dims = dims
    call h5screate_simple_f(1,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_vector_integer

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

  subroutine write_hdf5_vector_double(name,n,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: n
    real, intent(in)     :: a(n)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0

    dims(1) = n
    data_dims = dims
    call h5screate_simple_f(1,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_vector_double

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

  subroutine write_hdf5_string(name,n,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: n
    character(len=*), intent(in) :: a
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id, type_id
    integer(hsize_t) :: dims(7) = 0
    integer(size_t) :: typesize
    
    dims(1) = 1
    data_dims = dims
    
    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,error)
    typesize = n
    call h5tset_size_f(type_id,typesize,error)
    call h5screate_simple_f(1,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,type_id,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,type_id,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_string

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
  subroutine write_hdf5_2darray_integer(name,ni,nj,a,file_id)

    ! Write a 2d array of integer values to a HDF5 file

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj
    integer, intent(in)          :: a(ni,nj)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0
    
    dims(1) = ni
    dims(2) = nj
    data_dims = dims
    call h5screate_simple_f(2,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_2darray_integer

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

  subroutine write_hdf5_2darray_double(name,ni,nj,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj
    real, intent(in)          :: a(ni,nj)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0
    
    dims(1) = ni
    dims(2) = nj
    data_dims = dims
    call h5screate_simple_f(2,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_2darray_double

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

  subroutine write_hdf5_3darray_double(name,ni,nj,nk,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk
    real, intent(in)     :: a(ni,nj,nk)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0
    
    dims(1) = ni
    dims(2) = nj
    dims(3) = nk
    data_dims = dims
    call h5screate_simple_f(3,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_3darray_double

  subroutine write_hdf5_3darray_integer(name,ni,nj,nk,a,file_id)

    ! Write a 2d array of integer values to a HDF5 file

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj,nk
    integer, intent(in)          :: a(ni,nj,nk)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0
    
    dims(1) = ni
    dims(2) = nj
	dims(3) = nk
    data_dims = dims
    call h5screate_simple_f(3,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_3darray_integer

  subroutine write_hdf5_4darray_double(name,ni,nj,nk,nl,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk, nl
    real, intent(in)     :: a(ni,nj,nk,nl)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0
    
    dims(1) = ni
    dims(2) = nj
    dims(3) = nk
	dims(4) = nl
    data_dims = dims
    call h5screate_simple_f(4,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)
    
  end subroutine write_hdf5_4darray_double
  
  subroutine write_hdf5_5darray_double(name,ni,nj,nk,nl,nm,a,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in)          :: ni, nj, nk,nl,nm
    real, intent(in)     :: a(ni,nj,nk,nl,nm)
    integer(hid_t), intent(in)   :: file_id

    integer :: error
    integer(hsize_t) :: data_dims(7) = 0
    integer(hid_t) :: dspace_id, dset_id
    integer(hsize_t) :: dims(7) = 0

    dims(1) = ni
    dims(2) = nj
    dims(3) = nk
	dims(4) = nl
	dims(5) = nm
    data_dims = dims
    call h5screate_simple_f(5,dims,dspace_id,error)
    call h5dcreate_f(file_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
    call h5sclose_f(dspace_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a,data_dims,error)
    call h5dclose_f(dset_id,error)

  end subroutine write_hdf5_5darray_double

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

  subroutine copy_hdf5_scalar_integer(old_name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: old_name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    integer :: error
    integer :: a

    call read_hdf5_scalar_integer(old_name,a,file_id,error)
    call write_hdf5_scalar_integer(new_name,a,ofile_id)

  end subroutine copy_hdf5_scalar_integer

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

  subroutine copy_hdf5_scalar_double(name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    real :: a
    integer :: error
    
    call read_hdf5_scalar_double(name,a,file_id,error)
    call write_hdf5_scalar_double(new_name,a,ofile_id)
    
  end subroutine copy_hdf5_scalar_double

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

  subroutine copy_hdf5_vector_integer(name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    integer :: error
    integer :: n
    integer, allocatable :: a(:)
    
    call getdims_hdf5_vector(name,n,file_id)
    allocate(a(n))
    call read_hdf5_vector_integer(name,n,a,file_id,error)
    call write_hdf5_vector_integer(new_name,n,a,ofile_id)
    deallocate(a)
    
  end subroutine copy_hdf5_vector_integer

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

  subroutine copy_hdf5_vector_double(name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    integer :: error
    integer :: n
    real, allocatable :: a(:)
    
    call getdims_hdf5_vector(name,n,file_id)
    allocate(a(n))
    call read_hdf5_vector_double(name,n,a,file_id,error)
    call write_hdf5_vector_double(new_name,n,a,ofile_id)
    deallocate(a)
    
  end subroutine copy_hdf5_vector_double

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

  subroutine copy_hdf5_string(name,file_id,new_name,ofile_id,n)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id
    integer, intent(in)          :: n

    character(len=n) :: str
    integer :: error

    call read_hdf5_string(name,n,str,file_id,error)
    call write_hdf5_string(new_name,n,str,ofile_id)

  end subroutine copy_hdf5_string

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

  subroutine copy_hdf5_2darray_integer(name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    integer :: n(2)
    integer, allocatable :: a(:,:)
    integer :: error

    call getdims_hdf5_2darray(name,n,file_id)
    allocate(a(n(1),n(2)))
    call read_hdf5_2darray_integer(name,n(1),n(2),a,file_id,error)
    call write_hdf5_2darray_integer(new_name,n(1),n(2),a,ofile_id)
    deallocate(a)

  end subroutine copy_hdf5_2darray_integer

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

  subroutine copy_hdf5_2darray_double(name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    integer :: n(2)
    real, allocatable :: a(:,:)
    integer :: error

    call getdims_hdf5_2darray(name,n,file_id)
    allocate(a(n(1),n(2)))
    call read_hdf5_2darray_double(name,n(1),n(2),a,file_id,error)
    call write_hdf5_2darray_double(new_name,n(1),n(2),a,ofile_id)
    deallocate(a)

  end subroutine copy_hdf5_2darray_double

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

  subroutine copy_hdf5_3darray_double(name,file_id,new_name,ofile_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    character(len=*), intent(in) :: new_name
    integer(hid_t), intent(in)   :: ofile_id

    integer :: n(3)
    real, allocatable :: a(:,:,:)
    integer :: error

    call getdims_hdf5_3darray(name,n,file_id)
    allocate(a(n(1),n(2),n(3)))
    call read_hdf5_3darray_double(name,n(1),n(2),n(3),a,file_id,error)
    call write_hdf5_3darray_double(new_name,n(1),n(2),n(3),a,ofile_id)
    deallocate(a)

  end subroutine copy_hdf5_3darray_double

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

  subroutine getdims_hdf5_vector(name,n,file_id)
    
    implicit none
    character(len=*), intent(in)  :: name
    integer(hid_t), intent(in)    :: file_id
    integer, intent(out) :: n

    integer :: error
    integer(hid_t) :: dataset,dataspace
    integer(hsize_t) :: dims(7), dimsmax(1)
    
    call h5dopen_f(file_id,name,dataset,error)
    call h5dget_space_f(dataset,dataspace,error)
    call h5sget_simple_extent_dims_f(dataspace,dims,dimsmax,error)
    call h5sclose_f(dataspace,error)
    call h5dclose_f(dataset,error)
    n = dims(1)
    
  end subroutine getdims_hdf5_vector

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

  subroutine getdims_hdf5_2darray(name,n,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out) :: n(2)

    integer :: error
    integer(hid_t) :: dataset,dataspace
    integer(hsize_t) :: dims(2), dimsmax(2)
    
    call h5dopen_f(file_id,name,dataset,error)
    call h5dget_space_f(dataset,dataspace,error)
    call h5sget_simple_extent_dims_f(dataspace,dims,dimsmax,error)
    call h5sclose_f(dataspace,error)
    call h5dclose_f(dataset,error)
    n = dims
    
  end subroutine getdims_hdf5_2darray
  
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

  subroutine getdims_hdf5_3darray(name,n,file_id)

    implicit none
    character(len=*), intent(in) :: name
    integer(hid_t), intent(in)   :: file_id
    integer, intent(out) :: n(3)

    integer :: error
    integer(hid_t) :: dataset,dataspace
    integer(hsize_t) :: dims(3), dimsmax(3)
    
    call h5dopen_f(file_id,name,dataset,error)
    call h5dget_space_f(dataset,dataspace,error)
    call h5sget_simple_extent_dims_f(dataspace,dims,dimsmax,error)
    call h5sclose_f(dataspace,error)
    call h5dclose_f(dataset,error)
    n = dims
    
  end subroutine getdims_hdf5_3darray

end module hdf5_helpers
