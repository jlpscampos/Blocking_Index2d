PROGRAM MAIN

  USE netcdf
  USE Blocked_Flow_Index_new

  IMPLICIT NONE

  INTEGER,PARAMETER :: xdim = 360         ! Setting x dimension parameter
  INTEGER,PARAMETER :: ydim = 181         ! setting y dimension parameter
  INTEGER :: NX,NY,NT                     ! 
  INTEGER :: tdim                         ! Time dimension
  INTEGER :: BISSEXTO                     ! Leap year parameter
  REAL,DIMENSION(xdim) :: lon             ! Longitude 
  REAL,DIMENSION(ydim) :: lat             ! Latitude
  
  REAL,PARAMETER :: res = 1.0             ! grid resolution
  INTEGER,PARAMETER :: pers = 5           ! persistence for blocking events

  ! Allocatable variables index, difference, potential temperature, blocking 
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: ind, dif, theta,theta4, blk, thetas4
  INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: thet

  CHARACTER(300) :: iofile, var, att                      ! Character variable
  INTEGER :: t, y, month, k, days, time1day               ! time dependence 
  REAL ::   scale, offset                                 ! Netcdf library requirements

  INTEGER :: begin, final                                 ! Start and end dates
  CHARACTER(300) :: input_dir, output_dir, conf_file      ! Input and output folders
  
  conf_file = 'config_files/config.fortran.txt'
  CALL read_conf_file( conf_file, input_dir, output_dir, begin, final )  ! Getting some variables

  DO y = begin,final
    WRITE(*,'(A,1X,I4,A)') 'Doing year :',y,' ----------------------------> '
    k = 1
    
    IF( BISSEXTO(y) == 1 )THEN
      tdim = 366
    ELSE
      tdim = 365
    END IF
    
	DO month = 1,12
      IF( month < 10 )THEN  
        WRITE(iofile,'(A,I4,A,I4,A,I1,A)') TRIM(input_dir),y,'/theta_2PVU_',y,'-0',month,'.nc'
      ELSE
        WRITE(iofile,'(A,I4,A,I4,A,I2,A)') TRIM(input_dir),y,'/theta_2PVU_',y,'-',month,'.nc'
      END IF
      
      WRITE(*,'(A,A)') 'Reading file: ',TRIM(iofile)

      SELECT CASE( month )
        CASE(1); days = 31
        CASE(2); days = 28
        CASE(3); days = 31
        CASE(4); days = 30
        CASE(5); days = 31
        CASE(6); days = 30
        CASE(7); days = 31
        CASE(8); days = 31
        CASE(9); days = 30
        CASE(10);days = 31
        CASE(11);days = 30
        CASE(12);days = 31
      END SELECT

      ! Verify for leap year
      IF( BISSEXTO(y) == 1 .and. month == 2)THEN
        days = 29
      END IF
      
      
      IF( month==1 ) ALLOCATE(theta(xdim,ydim,tdim*4))
      ALLOCATE(thet(xdim,ydim,days*4))
      
      NX = xdim; NY = ydim; NT = 4*days;
      var = 'pt';           CALL reader(iofile,var,thet)
      att = 'scale_factor'; CALL reader_att( iofile,att,var,scale )
      att = 'add_offset'  ; CALL reader_att( iofile,att,var,offset )

      theta(:,:,k:k+days-1) = thet(:,:,:) * scale + offset
      k = k+4*days;

      !write(*,*) theta(60,46,k-1)

      DEALLOCATE(thet)

      IF( month == 1 )THEN
        var = 'longitude';   CALL reader_dim( iofile,var,lon, xdim ); 
        var = 'latitude' ;   CALL reader_dim( iofile,var,lat, ydim );
      END IF
      
    END DO

    ALLOCATE(theta4(xdim,ydim,tdim))
    time1day = 4
    CALL DAILY_MEANS(theta, theta4, xdim, ydim, 4*tdim, time1day)

    NX = xdim; NY = ydim; NT = tdim;

    ! Computing the blocking index
    ALLOCATE(ind(NX,NY,NT)); ALLOCATE(dif(NX,NY,NT)); ALLOCATE(blk(NX,NY,NT)); ALLOCATE(thetas4(NX,NY,NT))

    ! Before evaluate the blockings compute smooth the data in n degrees
    CALL SMOOTH(theta4,thetas4, xdim, ydim, tdim, res)
    
    DO t = 1, tdim
      CALL PH_Index ( ind(:,:,t), dif(:,:,t), thetas4(:,:,t), xdim, ydim, res )
      CALL GLBI ( ind(:,:,t), xdim, ydim, res )
    END DO
    
    WRITE(*,*)
    WRITE(*,*) 'Wait -----> Searching for blockings'
    WRITE(*,*)

    ! Computing Blockings with pers days
    CALL BLOCKING (blk, ind, xdim, ydim, tdim, res, pers )
   
    WRITE(iofile,'(A,i4,A)') TRIM(output_dir), y,'.blocked.nc'
    CALL write_netcdf( iofile,ind, lat, lon)
    WRITE(*,'(A,A)') 'Writing file: ',TRIM(iofile)
    
    WRITE(iofile,'(A,i4,A)') TRIM(output_dir), y,'.blocking.nc'
    CALL write_netcdf( iofile,blk, lat, lon)
    WRITE(*,'(A,A)') 'Writing file: ',TRIM(iofile)
    
    DEALLOCATE(ind); DEALLOCATE(blk); DEALLOCATE(dif); DEALLOCATE(theta);
    DEALLOCATE(theta4); DEALLOCATE(thetas4);

  END DO

  

   CONTAINS

   !_____________________________________________________________!
   SUBROUTINE read_conf_file( conf_file, input_dir, output_dir, begin, final)
       ! Subroutine to read configuration file
       
       CHARACTER(300), INTENT(IN)    :: conf_file
       CHARACTER(300), INTENT(INOUT) :: input_dir, output_dir
       INTEGER       , INTENT(INOUT) :: begin, final

       OPEN(10, file=conf_file )
       
       READ(10,'(a)') input_dir
       READ(10,'(a)') output_dir
       READ(10,'(i4)') begin
       READ(10,'(i4)') final
    
       CLOSE(10)

    END SUBROUTINE

   !_____________________________________________________________!
   SUBROUTINE reader(fname,varname,output)
      ! Subroutine to read the data from netcdf files
      
      CHARACTER(100),INTENT(in)                :: fname    ! file name
      CHARACTER(100),INTENT(in)                :: varname  ! variable to be readed
      INTEGER,DIMENSION(NX,NY,NT),INTENT(inout)   :: output   ! output variable
      INTEGER :: ncid,varid                                ! netCDF id and variables id

      !Open the file N90_NOWRITE tells netcdf wants read-only access to the file
      call check( nf90_open(trim(fname), NF90_NOWRITE, ncid) )  
      call check( nf90_inq_varid(ncid,trim(varname),varid) )   !Get the varid
      call check( nf90_get_var(ncid, varid, output) )          !Read the data
      call check( nf90_close(ncid) )                           !Close the file

    END SUBROUTINE reader

    !_____________________________________________________________!
    SUBROUTINE reader_dim( fname,varname,output, sz )
      ! This subroutine is only used to return the latitude and longitude
      INTEGER,INTENT(in)                     :: sz       ! size of the variable
      CHARACTER(100),INTENT(in)              :: fname    ! file name
      CHARACTER(100),INTENT(in)              :: varname  ! variable to be readed
      REAL,DIMENSION(sz),INTENT(out)         :: output   ! output variable (lat or lon)
      INTEGER :: ncid,varid                              ! netCDF id and variables id

      call check( nf90_open(trim(fname), NF90_NOWRITE, ncid) )  
      call check( nf90_inq_varid(ncid,trim(varname),varid) )   !Get the varid
      call check( nf90_get_var(ncid, varid, output) )          !Read the data
      call check( nf90_close(ncid) )                           !Close the file

    END SUBROUTINE reader_dim

    SUBROUTINE reader_att( fname,att,varname,output )
      ! this subroutine read the att of a netcdf file
      CHARACTER(len=100),INTENT(in)              :: fname    ! file name
      CHARACTER(len=100),INTENT(in)              :: att      ! variable to be readed
      REAL,INTENT(out)                       :: output   ! output variable (lat or lon)
      CHARACTER(100),INTENT(in)              :: varname  !
      INTEGER :: ncid,varid,attid, xtype,len             ! netCDF id and variables id

      call check( nf90_open(trim(fname), NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid,trim(varname),varid) )
      
      call check( nf90_inquire_attribute(ncid,varid,trim(att),xtype,len,attid) )
      call check( nf90_get_att(ncid,varid,trim(att),output) )
      call check( nf90_close(ncid) )
    
    END SUBROUTINE reader_att
    
    !_____________________________________________________________!
    SUBROUTINE write_netcdf( file_out ,div, lat, lon)
      ! This subroutine writes the otput into a netcdf file
 
      REAL,DIMENSION(NX,NY,NT),INTENT(in) :: div
      REAL,DIMENSION(NX),INTENT(in)       :: lon
      REAL,DIMENSION(NY),INTENT(in)       :: lat
      INTEGER :: x_dimid, y_dimid,t_dimid, z_dimid, varid , ncid
      INTEGER,DIMENSION(3) :: dimids, step, start
      REAL,DIMENSION(NT) :: time
      INTEGER :: i, lat_varid, lon_varid, div_varid, t_varid, NTclear
      CHARACTER(len=100), INTENT(in) :: file_out
      
      CHARACTER(len = *), PARAMETER :: UNITS = "units"
      CHARACTER(len = *), PARAMETER :: DIV_UNITS = "do not matter"
      CHARACTER(len = *), PARAMETER :: LAT_UNITS = "degrees_north"
      CHARACTER(len = *), PARAMETER :: LON_UNITS = "degrees_east"
      CHARACTER(len = *), PARAMETER :: TIME_UNITS= "Running Steps"

      CHARACTER(len = *), PARAMETER :: AXIS = 'axis'
      CHARACTER(len = *), PARAMETER :: XAXIS = "X"
      CHARACTER(len = *), PARAMETER :: YAXIS = "Y"
      CHARACTER(len = *), PARAMETER :: TAXIS = "T"

      time(1) = 0.
      DO i = 2, NTclear
        time(i) = time(i-1)+(1./4.)
      END DO

      call check( nf90_create(TRIM(file_out), NF90_CLOBBER, ncid) )    ! Create a netcdf file

      ! Define the dimensions
      call check( nf90_def_dim( ncid,'lat',NY, y_dimid) )          ! Define the dimensions
      call check( nf90_def_dim( ncid,'lon',NX, x_dimid) )          ! Define the dimensions
      call check( nf90_def_dim( ncid,'time',NT, t_dimid) )

      dimids = (/ x_dimid,y_dimid,t_dimid/)                        ! IDs arrays
      
      ! Define the variables
      call check( nf90_def_var(ncid,'lat',NF90_DOUBLE,y_dimid,lat_varid) )
      call check( nf90_def_var(ncid,'lon',NF90_DOUBLE,x_dimid,lon_varid) )
      call check( nf90_def_var(ncid,'time',NF90_DOUBLE,t_dimid,t_varid) )
      call check( nf90_def_var(ncid,'blk',NF90_DOUBLE, dimids, div_varid ) )

      ! Assign units attributes to coordinate variables
      call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS ) )
      call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS ) )
      call check( nf90_put_att(ncid, t_varid, UNITS, TIME_UNITS) )
      call check( nf90_put_att(ncid, div_varid, UNITS, DIV_UNITS) )

      call check( nf90_put_att(ncid, lat_varid, AXIS, YAXIS) )
      call check( nf90_put_att(ncid, lon_varid, AXIS, XAXIS) )
      call check( nf90_put_att(ncid, t_varid, AXIS, TAXIS) )

      ! End define mode
      call check( nf90_enddef(ncid) )

      ! Write coordinate variable data
      call check( nf90_put_var(ncid,lat_varid,lat) )
      call check( nf90_put_var(ncid,lon_varid,lon) )
      
      
      !Writing the data in the netcdf file
      step  = (/NX,NY,NT/)
      start = (/1,1,1/)
      call check( nf90_put_var(ncid,div_varid,div,start=start,count=step) )
      
      ! close the file
      call check( nf90_close(ncid) )

      write(*,*) 'Success file created'
 
    END SUBROUTINE write_netcdf
   
    !_____________________________________________________________!
    SUBROUTINE check(status)
       integer, intent(in) :: status

       if( status /= nf90_noerr )then
         write(*,*) trim(nf90_strerror(status))
         stop "Stopped"
       end if
     END SUBROUTINE check



END PROGRAM MAIN
  

INTEGER FUNCTION BISSEXTO( year )
  IMPLICIT NONE

  INTEGER,INTENT(in) :: year

  IF( mod(year,400) == 0 )THEN
    BISSEXTO = 1 ; RETURN
  ELSE IF( mod(year,4)==0 .and. mod(year,100)/=0 )THEN
    BISSEXTO = 1 ; RETURN
  ELSE
    BISSEXTO = 0 ; RETURN
  END IF
    
END FUNCTION
