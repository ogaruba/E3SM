module prescribed_cloud

!--------------------------------------------------------------------------
! Purpose:
!
! Reads cloud-related fields, puts the them into the physics buffer for use
! by radiation
! This code was written by Bryce Harrop with significant contributions from 
! Balwinder Singh, Brian Medeiros, and Jerry Olson
!
!--------------------------------------------------------------------------

  use shr_kind_mod,     only : r8 =>shr_kind_r8, cx =>SHR_KIND_CX, cl =>SHR_KIND_CL, &
                               cs =>SHR_KIND_CS, cxx =>SHR_KIND_CXX
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog
  use input_data_utils, only : time_coordinate
  use shr_log_mod ,     only : errMsg => shr_log_errMsg

  implicit none
  private
  save

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_cloud_init
  public :: prescribed_cloud_adv
  public :: write_prescribed_cloud_restart
  public :: read_prescribed_cloud_restart
  public :: has_prescribed_cloud
  public :: prescribed_cloud_register
  public :: init_prescribed_cloud_restart
  public :: prescribed_cloud_readnl

  logical :: has_prescribed_cloud = .false.
  integer          , parameter :: nflds             = 10
  character(len=16), parameter :: cloud_name(nflds) = (/'DEI_rad'   ,'MU_rad'  ,'LAMBDAC_rad' ,'ICIWP_rad' ,'ICLWP_rad' ,'DES_rad' , &
                                                        'ICSWP_rad' ,'CLD_rad', 'CLDFSNOW_rad', 'CONCLD_rad' /)

  character(len=16)  :: fld_name(nflds)             = (/'DEI_rad'  ,'MU_rad' ,'LAMBDAC_rad','ICIWP_rad','ICLWP_rad','DES_rad', &
                                                        'ICSWP_rad','CLD_rad', 'CLDFSNOW_rad', 'CONCLD_rad'/)
  character(len=256) :: filename                    = ' '
  character(len=256) :: filelist                    = ' '
  character(len=256) :: datapath                    = ' '
  logical            :: rmv_file                    = .false.
  character(len=32)  :: specifier(nflds)            = ''

  logical, parameter :: horz_native = .true. ! they will all have the same behavior if they
                                             ! all come from the same file
  logical            :: dimnames_set = .false. 
  integer            :: number_flds
  character(len=16)  :: spc_name_list(nflds)
  character(len=cl)  :: spc_flist(nflds),spc_fnames(nflds)
  character(len=8)   :: dim1name, dim2name
  character(len=32)  :: data_type = 'CYCLICAL'
  real(r8)           :: num_file_years = 0._r8
  character(len=32)  :: air_type  = 'CYCLICAL_LIST'

  !------------------------------------------------------------------
  !DEFINITION:
  !"native grid forcing file": A forcing file which has to be on the 
  !same grid horizontally as the model is running on. For example, 
  !if the model is running on ne30 grid, forcing file has to be on
  !ne30 grid horizontally.
  !------------------------------------------------------------------
  
  type :: forc_air_native_grid
     !------------------------------------------------------------------
     !"forc_air_native_grid" is forcing from files which has to be on the
     !native grid (only in horizontal,vertical resolution may be different 
     !from the model's grid resolution)
     !That is, forcing files has to be on same grid as the grid used for 
     !the model run
     !------------------------------------------------------------------
     
     !Number of levels in the 3D forcing file
     integer                               :: lev_frc
     
     !Data structure to store two time samples from a file to do time interpolation in the next step
     !(pcols,lev_frc,begchunk:endchunk,2)
     real(r8), pointer, dimension(:,:,:,:) :: native_grid_flds_tslices
     
     !Data structure to store data after time interpolation from two time samples
     !(pcols,lev_frc,begchunk:endchunk)
     real(r8), pointer, dimension(:,:,:)   :: native_grid_flds
     
     !Data structure to keep track of time
     type(time_coordinate) :: time_coord
     
     !specie name
     character( len = cx)  :: spc_name_ngrd
     character( len = cx)  :: spc_cname_ngrd
     character( len = cx)  :: spc_fname_ngrd


     !Level bounds read from input file
     real(r8), pointer, dimension(:,:) :: lev_bnds

     !Forcing file name
     character( len = cx)  :: input_file
     
     !Units of forcing data
     character( len = cs)  :: units
     
     !logical to control first data read
     logical               :: initialized
     
     !pbuf index to store read in data in pbuf
     integer               :: pbuf_ndx = -1
  end type forc_air_native_grid
  type(forc_air_native_grid),allocatable :: native_grid_cloud(:)

  type :: forcing_air
     real(r8)              :: mw
     character(len=cl)     :: filelist
     character(len=cl)     :: filename
     real(r8), pointer     :: times(:)
     real(r8), pointer     :: levi(:)
     character(len=11)     :: species
     character(len=8)      :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld),pointer       :: fields(:)
     type(trfile)              :: file
  end type forcing_air
  
  type(forcing_air), allocatable :: forcings_air(:)

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_cloud_register()
    use ppgrid,         only: pver, pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: i,idx

    if (has_prescribed_cloud) then
       do i = 1,nflds
          call pbuf_add_field(cloud_name(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
       enddo
    endif

  endsubroutine prescribed_cloud_register


!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_cloud_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_cloud_readnl'

   character(len=256) :: prescribed_cloud_file
   character(len=256) :: prescribed_cloud_filelist
   character(len=256) :: prescribed_cloud_datapath
   character(len=32)  :: prescribed_cloud_type
   real(r8)           :: prescribed_cloud_num_file_years

   namelist /prescribed_cloud_nl/ &
      prescribed_cloud_file,      &
      prescribed_cloud_filelist,  &
      prescribed_cloud_datapath,  &
      prescribed_cloud_type,      &
      prescribed_cloud_num_file_years
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_cloud_file     = filename
   prescribed_cloud_filelist = filelist
   prescribed_cloud_datapath = datapath
   prescribed_cloud_type     = data_type
   prescribed_cloud_num_file_years= num_file_years

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_cloud_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_cloud_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_cloud_file,     len(prescribed_cloud_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_filelist, len(prescribed_cloud_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_datapath, len(prescribed_cloud_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_type,     len(prescribed_cloud_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_cloud_num_file_years, 1, mpir8, 0, mpicom)
#endif

   ! Update module variables with user settings.
   filename   = prescribed_cloud_file
   filelist   = prescribed_cloud_filelist
   datapath   = prescribed_cloud_datapath
   data_type  = prescribed_cloud_type
   num_file_years = prescribed_cloud_num_file_years

   ! Turn on prescribed cloud if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_cloud = .true.

end subroutine prescribed_cloud_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------

  subroutine init_prescribed_cloud_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_cloud', piofile, file )

  end subroutine init_prescribed_cloud_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_cloud_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_cloud_restart

!-------------------------------------------------------------------
  subroutine read_prescribed_cloud_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_cloud', piofile, file )

  end subroutine read_prescribed_cloud_restart
!================================================================================================


  subroutine prescribed_cloud_init(state, pbuf2d)
    !-------------------------------------------------------------------
    ! **** Initialize the aircraft aerosol data handling ****
    ! called by:
    !-------------------------------------------------------------------
    use tracer_data,      only: trcdata_init
    use physics_types,    only: physics_state
    use ppgrid,           only: begchunk, endchunk, pcols
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_index
    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc
    use cam_grid_support, only: cam_grid_id, cam_grid_check
    use cam_grid_support, only: cam_grid_get_dim_names
    use dyn_grid,         only: get_horiz_grid_dim_d
    use dycore,           only: dycore_is
    use cam_pio_utils,    only: cam_pio_openfile
    use pio,              only: file_desc_t, pio_nowrite, pio_closefile, pio_inq_dimid, pio_bcast_error, &
         pio_seterrorhandling, pio_noerr, pio_inquire_dimension, pio_get_att, pio_inq_varid, pio_get_var
    
    
    implicit none
    
!BEH: check to see if state and pbuf2d really need to be passed in.  I didn't use them
!     in the init routine for prescrbed surface fluxes.
    !arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    
    !local vars
    type(file_desc_t)   :: fh
    character(len=16)   :: spc_name
    character(len=16)   :: spc_cname
    character(len=16)   :: spc_fname
    character(len=cxx)  :: err_str
    
    integer :: ndx, istat, i, astat, m, n, mm, c
    integer :: grid_id
    integer :: dimlevid, var_id, errcode, dim1id, dim2id, dim1len, dim2len
    integer :: dimbndid, nbnd
    integer :: hdim1_d, hdim2_d    ! model grid size
    real(r8) :: dtime
    logical  :: fixed, cyclical
    
    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------

    if ( has_prescribed_cloud ) then
       if ( masterproc ) then
          write(iulog,*) 'now cloud is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    do i = 1,nflds
       specifier(i) = trim(cloud_name(i))
       if ( masterproc ) then
          write(iulog,*) 'A specifier:'//specifier(i)
       endif
    end do
    
    
    !------------------------------------------------------------------
    ! For forcing files which has to be on the native grid,dimensions
    ! are set in the following if condition
    !------------------------------------------------------------------
    if( horz_native ) then
       if (.not. dimnames_set) then
          grid_id = cam_grid_id('physgrid')
          if (.not. cam_grid_check(grid_id)) then          
             call endrun('no "physgrid" grid:'//errmsg(__FILE__,__LINE__))
          endif
          !dim1name and dim2name are populated here with the grid dimension the model is running on (e.g. ne30, lat, lon etc.)
          !For SE grid, dim1name = dim2name = "ncol"
          !For FV grid, dim1name = lon, dim2name = lat
          call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
          dimnames_set = .true.
       end if
       
       !--------------------------------------------------------------------------------
       ! allocate forcings type array for native grid forcing files
       !--------------------------------------------------------------------------------

       allocate( native_grid_cloud(nflds), stat=astat )

       if( astat /= 0 ) then 
          write(err_str,*) 'failed to allocate native_grid_cloud array; error = ',astat,',',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       end if
    endif
    
    if (masterproc) write(iulog,*) ' '
    
!    if( any( .not. horz_native(:) ) ) then
    if( .not. horz_native ) then
       !-----------------------------------------------------------------------
       !       allocate forcings type array
       !-----------------------------------------------------------------------

       allocate( forcings_air(nflds), stat=astat )

       if( astat/= 0 ) then
          write(err_str,*) 'failed to allocate forcings_air array; error = ',astat,',',errmsg(__FILE__, __LINE__)
          call endrun(err_str) 
       end if
    endif
    
    
    !-----------------------------------------------------------------------
    !       setup the forcings_air type array
    !-----------------------------------------------------------------------

    species_loop : do m = 1,nflds
       
       spc_name = specifier(m)

       spc_cname = trim(cloud_name(m))
       spc_fname = trim(cloud_name(m))
       if ( masterproc ) then
          write(iulog,*) 'spc_name is: '//trim(spc_name)
          write(iulog,*) 'spc_cname is: '//trim(spc_cname)
          write(iulog,*) 'spc_fname is: '//trim(spc_fname)
       endif
       
!       if( horz_native(index_map(m))) then
       if( horz_native ) then
          !-----------------------------------------------------------------------
          !       initialize variables for native grid forcing files
          !-----------------------------------------------------------------------
       
          native_grid_cloud(m)%spc_name_ngrd  = spc_name
          native_grid_cloud(m)%spc_cname_ngrd = spc_cname
          native_grid_cloud(m)%spc_fname_ngrd = spc_fname

          fixed    = .false.
          cyclical = .false.
          select case ( data_type )
          case( 'FIXED' )
             fixed = .true.
          case( 'CYCLICAL' )
             cyclical = .true.
          case( 'SERIAL' )
             ! Do nothing
          case default 
             write(iulog,*) 'prescribed_cloud: invalid data type: '//trim(data_type)//' file: '//trim(filename)
             write(iulog,*) 'prescribed_cloud: valid data types: SERIAL | CYCLICAL | FIXED '
             call endrun('prescribed_cloud: invalid data type: '//trim(data_type)//' file: '//trim(filename))
          endselect

          native_grid_cloud(m)%input_file     = trim(datapath)//'/'//trim(filename)

          native_grid_cloud(m)%initialized    = .false.
          dtime  = 0.0_r8
          call native_grid_cloud(m)%time_coord%initialize(trim(adjustl(native_grid_cloud(m)%input_file)), &
               force_time_interp=.true., delta_days=dtime, fixed=fixed, &
               cyclical=cyclical, num_file_years=num_file_years)

          !-----------------------------------------------------------------------
          !       Open file
          !-----------------------------------------------------------------------
          call cam_pio_openfile(fh, trim(adjustl(native_grid_cloud(m)%input_file)), PIO_NOWRITE) 

          !ask PIO to return the control if it experiences an error so that we can 
          !handle it explicitly in the code
          call pio_seterrorhandling(fh, pio_bcast_error)
          
          !-----------------------------------------------------------------------
          !       Sanity checks for the native grid
          !-----------------------------------------------------------------------
          
          !if forcing file is on a different grid than the model grid
          !(e.g. model is running on an FV grid and forcing netcdf file is on an SE grid), exit with an error
          if(pio_inq_dimid(fh, trim(adjustl(dim1name)), dim1id) /= pio_noerr) then
             !pio_inq_dimid function tries to find dim1name in file with id "fh"
             !if it can't find dim1name, it means there is a mismatch in model and netcdf
             !file grid
             call endrun('grid mismatch, failed to find '//dim1name//' dimension in file:'&
                  ' '//trim(adjustl(native_grid_cloud(m)%input_file))//' '&
                  ' '//errmsg(__FILE__,__LINE__))
          endif
          
          !find if the model and netcdf file has same grid resolution
          call get_horiz_grid_dim_d(hdim1_d,hdim2_d) !get model dim lengths
          if( dycore_is('SE') )  then
             if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr) then
                if(dim1len /= hdim1_d ) then !compare model grid length with file's
                   write(err_str,*)'Netcdf file grid size(',dim1len,') should be same as model grid size(',&
                        hdim1_d,'), netcdf file is:'//trim(adjustl(native_grid_cloud(m)%input_file))
                   call endrun(err_str//errmsg(__FILE__,__LINE__))
                endif
             else
                call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_cloud(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
             endif
          elseif( dycore_is('LR')) then
             if(pio_inq_dimid(fh, trim(adjustl(dim2name)), dim2id)) then !obtain lat dimension of model
                call endrun('failed while inquiring dimension'//trim(adjustl(dim2name))//' from file:'&
                     ' '//trim(adjustl(native_grid_cloud(m)%input_file))//' '//errmsg(__FILE__,__LINE__))
             endif
             if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr .and. &
                  pio_inquire_dimension(fh, dim2id, len = dim2len) ==  pio_noerr) then !compare grid and model's dims
                if(dim1len /= hdim1_d .or. dim2len /= hdim2_d)then
                   write(err_str,*)'Netcdf file grid size(',dim1len,' x ',dim2len,') should be same as model grid size(',&
                        hdim1_d,' x ',hdim2_d,'), netcdf file is:'//trim(adjustl(native_grid_cloud(m)%input_file))
                   call endrun(err_str//errmsg(__FILE__,__LINE__))
                endif
             else
                call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_cloud(m)%input_file))//'&
                     &'//errmsg(__FILE__,__LINE__))
             endif
          else
             call endrun('Only SE or LR(FV) grids are supported currently:'//errmsg(__FILE__,__LINE__))
          endif
          
          !Find the value of vertical levels in the forcing file
          if( pio_inq_dimid(fh, 'lev', dimlevid) ==  pio_noerr ) then
             if ( pio_inquire_dimension(fh, dimlevid, len =  native_grid_cloud(m)%lev_frc) /=  pio_noerr ) then
                write(err_str,*)'failed to obtain value of "lev" dimension from file:',&
                     trim(adjustl(native_grid_cloud(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'Dimension "lev" is not found in:',&
                  trim(adjustl(native_grid_cloud(m)%input_file)),',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !get units of the data in the forcing file
          if(pio_inq_varid( fh, spc_fname, var_id ) == pio_noerr ) then
             if(pio_get_att( fh, var_id, 'units', native_grid_cloud(m)%units) .ne. pio_noerr ) then
                write(err_str,*)'failed to obtain units of variable ',trim(spc_fname),' in &
                     &file:',trim(adjustl(native_grid_cloud(m)%input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'variable ',trim(spc_fname),' not found in:',trim(adjustl(native_grid_cloud(m)%input_file)), &
                  ',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !close file
          call pio_closefile(fh)
          
          !allocate arrays to store data for interpolation in time
          allocate(native_grid_cloud(m)%native_grid_flds_tslices(pcols, native_grid_cloud(m)%lev_frc, &
               begchunk:endchunk,2), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'failed to allocate native_grid_cloud(',m,')%native_grid_flds_tslices array;&
                  error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !allocate arrays to hold data before the vertical interpolation
          allocate(native_grid_cloud(m)%native_grid_flds(pcols, native_grid_cloud(m)%lev_frc,begchunk:endchunk), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'failed to allocate native_grid_cloud(',m,')%native_grid_flds array; error = ',&
                  astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
          !get pbuf index to store the field in pbuf
          native_grid_cloud(m)%pbuf_ndx = pbuf_get_index(spc_cname,errcode)
          if(errcode < 0 ) then
             write(err_str,*)'failed to get pbuf index for specie:',spc_cname,' errorcode is:',errcode,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
          
       else
          
          allocate( forcings_air(m)%sectors(1), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'aircraft_emit_init: failed to allocate forcings_air%sectors &
                  &array; error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          end if
          
          allocate( forcings_air(m)%fields(1), stat=astat )
          if( astat/= 0 ) then
             write(err_str,*) 'aircraft_emit_init: failed to allocate forcings_air%fields &
                  &array; error = ',astat,',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          end if
          
          !-----------------------------------------------------------------------
          !         default settings
          !-----------------------------------------------------------------------
          forcings_air(m)%file%stepTime        = .true.   ! Aircraft data is not to be interpolated in time
          forcings_air(m)%file%cyclical_list   = .true.   ! Aircraft data cycles over the filename list
          forcings_air(m)%file%weight_by_lat   = .true.   ! Aircraft data -  interpolated with latitude weighting
          forcings_air(m)%file%conserve_column = .true.   ! Aircraft data - vertically interpolated to conserve the total column
          forcings_air(m)%species              = spc_name
          forcings_air(m)%sectors              = spc_name ! Only one species per file for aircraft data
          forcings_air(m)%nsectors             = 1
          forcings_air(m)%filelist             = spc_flist(m)
          forcings_air(m)%filename             = spc_fnames(m)
       endif
       
    end do species_loop

    
    
    !------------------------------------------------------------------
    !       Initialize the aircraft file processing
    !------------------------------------------------------------------

    do m=1,nflds

       number_flds = 0
       if(horz_native) then
          if (associated(native_grid_cloud(m)%native_grid_flds_tslices)) &
               number_flds = 1
          !read the forcing file once to initialize variable including time cordinate
          call advance_native_grid_data( native_grid_cloud(m) )
          native_grid_cloud(m)%initialized = .true.
       else
          allocate (forcings_air(m)%file%in_pbuf(size(forcings_air(m)%sectors)))
          forcings_air(m)%file%in_pbuf(:) = .true.
          if (associated(forcings_air(m)%fields)) number_flds = size( forcings_air(m)%fields )
          
          call trcdata_init( forcings_air(m)%sectors, forcings_air(m)%filename, forcings_air(m)%filelist, datapath, &
               forcings_air(m)%fields, forcings_air(m)%file, rmv_file, 0, 0, 0, air_type)
       endif
       
       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(err_str,*) 'There are no aircraft aerosols ',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
       end if
    end do
  end subroutine prescribed_cloud_init

  subroutine prescribed_cloud_adv( state, pbuf2d)
    !-------------------------------------------------------------------
    ! **** Advance to the next aircraft data ****
    ! called by:
    !-------------------------------------------------------------------
    
    use perf_mod,       only: t_startf, t_stopf
    use tracer_data,    only: advance_trcdata
    use physics_types,  only: physics_state
    use ppgrid,         only: begchunk, endchunk
    use ppgrid,         only: pcols, pver
    use string_utils,   only: to_lower, GLC
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    
    implicit none
    
    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    
    integer             :: ind, c, ncol, i, caseid, m, pbuf_ndx
    real(r8), pointer   :: tmpptr(:,:)
    character(len = cs) :: units_spc
    
    !------------------------------------------------------------------
    ! Return if aircraft_cnt is zero (no aircraft data to process)
    !------------------------------------------------------------------

    if( .not. has_prescribed_cloud ) return

    call t_startf('All_aircraft_emit_adv')
    
    !-------------------------------------------------------------------
    !    For each field, read more data if needed and interpolate it to the current model time
    !-------------------------------------------------------------------

    do m = 1,nflds

       if (horz_native) then
          units_spc = native_grid_cloud(m)%units
          pbuf_ndx  = native_grid_cloud(m)%pbuf_ndx
          
          !read in next time slice (if needed) and interpolate in time
          !following call just reads in time slices in horizontal
          !vertical interpolation is done in the next call
          call advance_native_grid_data( native_grid_cloud(m) )

!++BEH ---->  Testing putting data into pbuf
          !Need to store the native_grid_cloud data into the physics buffer.
          !$OMP PARALLEL DO PRIVATE (C, NCOL, TMPPTR, PBUF_CHNK)
          do c = begchunk, endchunk
             ncol = state(c)%ncol
             pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
             call pbuf_get_field(pbuf_chnk, pbuf_ndx, tmpptr)
             tmpptr(:ncol,:) = native_grid_cloud(m)%native_grid_flds(:,:,c)
          enddo
          
       else          
          units_spc = forcings_air(m)%fields(i)%units
          pbuf_ndx  = forcings_air(m)%fields(i)%pbuf_ndx
          call advance_trcdata( forcings_air(m)%fields, forcings_air(m)%file, state, pbuf2d)
       endif

    enddo
    
    call t_stopf('All_aircraft_emit_adv')
  end subroutine prescribed_cloud_adv

  subroutine advance_native_grid_data( native_grid_strct )
    !-------------------------------------------------------------------
    !    This subroutine reads the data from the native grid and
    !    interpolates in time
    ! called by:
    !-------------------------------------------------------------------
    
    use ppgrid,         only: begchunk, endchunk, pcols
    use ncdio_atm,      only: infld
    use cam_pio_utils,  only: cam_pio_openfile
    use pio,            only: file_desc_t, pio_nowrite, pio_closefile

    implicit none

    !args
    type(forc_air_native_grid), intent (inout) :: native_grid_strct 
    
    !local vars
    type(file_desc_t) :: fh
    character(len=cs) :: spc_name
    character(len=cs) :: spc_cname
    character(len=cs) :: spc_fname
    
    logical  :: read_data
    integer  :: indx2_pre_adv
    logical  :: found

    !obtain name of the specie
    spc_name  = native_grid_strct%spc_name_ngrd
    spc_cname = native_grid_strct%spc_cname_ngrd
    spc_fname = native_grid_strct%spc_fname_ngrd
    
    !Decide whether to read new data or not (e.g. data may needs to be read on month boundaries )
    read_data = native_grid_strct%time_coord%read_more() .or. .not. native_grid_strct%initialized
    
    !Find time index to decide whether to read new data or recycle previously read data
    indx2_pre_adv = native_grid_strct%time_coord%indxs(2)
    
    !compute weights for time interpolation (time_coord%wghts) by advancing in time
    call native_grid_strct%time_coord%advance()
    
    if ( read_data ) then
       
       !open file
       call cam_pio_openfile(fh, trim(adjustl(native_grid_strct%input_file)), PIO_NOWRITE)
       
       ! read time-level 1
       if (native_grid_strct%initialized .and. native_grid_strct%time_coord%indxs(1) == indx2_pre_adv) then
          ! skip the read if the needed vals for time level 1 are present in time-level 2
          native_grid_strct%native_grid_flds_tslices(:,:,:,1) = native_grid_strct%native_grid_flds_tslices(:,:,:,2)
       else
          !NOTE: infld call doesn't do any interpolation in space, it just reads in the data
          call infld(trim(spc_fname), fh, dim1name, dim2name, 'lev',&
               1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
               native_grid_strct%native_grid_flds_tslices(:,:,:,1), found, &
               gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(1))
          if (.not. found) then
             call endrun(trim(spc_fname) // ' not found '//errmsg(__FILE__,__LINE__))
          endif
       endif
       
       ! read time level 2
       call infld(trim(spc_fname), fh, dim1name, dim2name, 'lev',&
            1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
            native_grid_strct%native_grid_flds_tslices(:,:,:,2), found, &
            gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(2))
       
       if (.not. found) then
          call endrun(trim(spc_fname) // ' not found '//errmsg(__FILE__,__LINE__))
       endif
       
       !close file
       call pio_closefile(fh)
    endif
    ! interpolate between time-levels
    ! If time:bounds is in the dataset, and the dataset calendar is compatible with EAM's,
    ! then the time_coordinate class will produce time_coord%wghts(2) == 0.0,
    ! generating fluxes that are piecewise constant in time.
    
    if (native_grid_strct%time_coord%wghts(2) == 0.0_r8) then
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1)
    else
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1) + &
            native_grid_strct%time_coord%wghts(2) * (native_grid_strct%native_grid_flds_tslices(:,:,:,2) - &
            native_grid_strct%native_grid_flds_tslices(:,:,:,1))
    endif
    
  end subroutine advance_native_grid_data

end module prescribed_cloud
