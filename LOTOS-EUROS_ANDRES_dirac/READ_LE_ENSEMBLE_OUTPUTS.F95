  !==Program to Read LE Ensemble Outputs===
  !==Andres Yarce, Santiago Lopez-Restrepo==

program read_le_ensemble_outputs
  use module_EnKF
  use module_matrix
  use netcdf
  implicit none

  integer::time,Ens,j,i
  integer::lat,lon,level
  real*8,allocatable,dimension(:,:,:,:,:)::xb(:,:,:,:,:) !Dimensions lon,lat,hist,noise,time
  real*8,allocatable,dimension(:,:,:)::Y_xb(:,:,:),Y(:,:,:) !Dimensions lon,lat,time
  character (len=100)::FILE_NAME,Path,FILE_NAME2
  integer :: ncid, varid,status,latid,lonid,levelid,timeid,lat_varid
  real*8,allocatable,dimension(:)::latitude

! ============ Read information about NC files,location and number of ensembles ============ Modify with the template of the real ensemble members
  open(61,file='./temp/Ensembles.in',status='old')
  read(61,"(A200)") Path
  read(61,*) Ens
  close(61)
  loop_Read_Ens: do j=1,Ens
     if (j < 10) then            
	    write (FILE_NAME2,'(i1)') j ! converting integer to string using a 'internal file'
	    FILE_NAME=trim(path)//"/Ens_x0"//trim(FILE_NAME2)//'.nc'
     else            
	    write (FILE_NAME2,'(i2)') j ! converting integer to string using a 'internal file'	   
	    FILE_NAME=trim(path)//"/Ens_x"//trim(FILE_NAME2)//'.nc'
     endif
      write(*,*) 'Writing Ensemble '//trim(FILE_NAME2)
   !Get for the dimensions ID 
    status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
    status = nf90_inq_dimid(ncid, "latitude", latid)
    status = nf90_inq_dimid(ncid, "longitude", lonid)     
    status = nf90_inquire(ncid, unlimiteddimid = timeid)

    !Get dimensions length
    status = nf90_inquire_dimension(ncid, latid, len = lat)
    status = nf90_inquire_dimension(ncid, lonid, len = lon)
    status = nf90_inquire_dimension(ncid, timeid, len = time)

    !Write dimensions in file for next scripts
    if (j .eq. 1) then
	allocate(xb(lon,lat,1,1,time))
    allocate(latitude(lat))
    open(61,file='./data/dimensions.dat',status='unknown')
	write(61,*)lat
	write(61,*)lon
	write(61,*)time
	close(61)
    endif

   ! Get the varid of the data variable, based on its name.
    status= nf90_inq_varid(ncid, "dc", varid)
    status= nf90_inq_varid(ncid, "latitude", lat_varid)
    !Read the data.
    status=nf90_get_var(ncid, lat_varid, latitude)
	
    ! Read the data.
    status=nf90_get_var(ncid, varid, xb)
    status= nf90_close(ncid) 
    if (j<10) then
    	open(61,file='./data/Ens_x0'//trim(FILE_NAME2)//'.dat',status='unknown')
    else
	open(61,file='./data/Ens_x'//trim(FILE_NAME2)//'.dat',status='unknown')
    endif
    write(61,*) xb(:,:,1,1,:)
    close(61)
  end do loop_Read_Ens




! ============ Read information about NC files of Ensemble outputs ============
allocate(Y_xb(lon,lat,time))
loop_Read_Y_Ens: do j=1,Ens
     if (j < 10) then            
	    write (FILE_NAME2,'(i1)') j ! converting integer to string using a 'internal file'
	    FILE_NAME=trim(path)//"/Y_Ens_x0"//trim(FILE_NAME2)//'.nc'
     else            
	    write (FILE_NAME2,'(i2)') j ! converting integer to string using a 'internal file'	   
	    FILE_NAME=trim(path)//"/Y_Ens_x"//trim(FILE_NAME2)//'.nc'
     endif
      write(*,*) 'Writing Ensemble Output '//trim(FILE_NAME2)
   !Get for the dimensions ID 
    status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)

   ! Get the varid of the data variable, based on its name.
    status= nf90_inq_varid(ncid, "no2_column", varid)

    ! Read the data.
    status=nf90_get_var(ncid, varid, Y_xb(:,:,:))
 
    status= nf90_close(ncid)
    if (j<10) then
    	open(61,file='./data/Y_Ens_x0'//trim(FILE_NAME2)//'.dat',status='unknown')
    else
	open(61,file='./data/Y_Ens_x'//trim(FILE_NAME2)//'.dat',status='unknown')
    endif
    write(61,*) Y_xb(:,:,:)
    close(61)
  end do loop_Read_Y_Ens


! ============ Read information about NC files of observations ============ Modify with the template of the real observations
    allocate(Y(lon,lat,time))
    FILE_NAME=trim(path)//'/Y.nc'
    write(*,*) 'Writing Observations'
      

   ! Get the varid of the data variable, based on its name.
    status= nf90_inq_varid(ncid, "no2_column", varid)

    ! Read the data.
    status=nf90_get_var(ncid, varid, y(:,:,:))
    status= nf90_close(ncid)
   
    open(61,file='./data/Y.dat',status='unknown')
    write(61,*) Y(:,:,:)
    close(61)
     write(*,*) 'mark here'

  

end program read_le_ensemble_outputs
