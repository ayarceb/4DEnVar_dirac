  !==Program to Read LE Ensemble Outputs===
  !==Andres Yarce, Santiago Lopez-Restrepo==

program read_le_ensemble_outputs
  use module_EnKF
  use module_matrix
  use netcdf
  implicit none

  integer::time,Ens,j,i,k,h
  integer::lat,lon,level,initial_day,dawindows,days_simulation
  real,allocatable,dimension(:,:,:,:,:)::aux(:,:,:,:,:) !Dimensions lon,lat,hist,noise,time
  real*8,allocatable,dimension(:,:)::xb(:,:),aux_y(:,:)
  real,allocatable,dimension(:,:,:)::Y_xb(:,:,:),Y(:,:,:) !Dimensions lon,lat,time
  character (len=500)::FILE_NAME,Path,FILE_NAME2,FILE_NAME3,runid
  integer :: ncid,varid,status,latid,lonid,levelid,timeid,lat_varid,n_cells_int,indi,pixelid,pixel_len,ilatid,ilonid
  integer,allocatable,dimension(:) :: nx_int(:),ny_int(:),pixel(:),ilat(:),ilon(:)
  real,allocatable,dimension(:)::latitude
  logical :: file_exists
  integer,allocatable,dimension(:) :: agrup(:)
  
  !=====Numero de Celdas====
  n_cells_int=234  ! Andres_ este se debe modificar si se ponen mas 
	allocate(nx_int(n_cells_int))
	allocate(ny_int(n_cells_int))
	allocate(agrup(n_cells_int))
	
	
!============ Read information about NC files,location and number of ensembles ============ Modify with the template of the real ensemble members
  open(61,file='./temp/Ensembles.in',status='old')
  read(61,"(A200)") Path
  read(61,*) Ens
  close(61)
!============ Read Initial time and Assimilation Windows and runid======
  open(61,file='./DATA_4DEnVAR/parameters.in',status='old')
  read(61,*) dawindows
  close(61)
  
  open(61,file='./DATA_4DEnVAR/startdate.in',status='old')
  read(61,*) initial_day
  read(61,*) days_simulation
  close(61)
  
  open(61,file='./DATA_4DEnVAR/runid.in',status='old')
  read(61,*) runid
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
    write(*,*) 'Opening file:', FILE_NAME
    status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
    if (status .NE. 0) then
        Write(*,*) 'Yarce, no se pudo abrir el archivo ojo'
    endif
    status = nf90_inq_dimid(ncid, "latitude", latid)
    status = nf90_inq_dimid(ncid, "longitude", lonid)     
    status = nf90_inquire(ncid, unlimiteddimid = timeid)
    
    !Get dimensions length
    status = nf90_inquire_dimension(ncid, latid, len = lat)
    status = nf90_inquire_dimension(ncid, lonid, len = lon)
    status = nf90_inquire_dimension(ncid, timeid, len = time)

    !Write dimensions in file for next scripts
    if (j .eq. 1) then
			allocate(aux(lon,lat,1,1,time))
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
    status=nf90_get_var(ncid, varid, aux)
! nx_int=(/54,43,16,75,59,83,89,104,55,53,54,55,55,54,43,44,44,16,17,74,75,60,60,83,89,105/)
   ! ny_int=(/36,111,90,85,32,12,42,56,36,36,37,37,35,35,112,112,111,89,89,84,84,32,33,11,41,56/)

   ! nx_int=(/54,43,16,75,59,83,89,104,55,53,54,55,55,54,43,44,44,16,17,74,75,60,60,83,89,105,56,57/)
   ! ny_int=(/36,111,90,85,32,12,42,56,36,36,37,37,35,35,112,112,111,89,89,84,84,32,33,11,41,56,38,37/)

   ! Bogota (64,103),Medellin,Cali,       Barranquilla,Cartagena,Cucuta, Ibague, Villavicencia, Bucaramanga, SantaMarta, Valledupar, Pereira, Monteria, Pasto, Manizales, Caracas, Maracaibo, Valencia, Barquisimeto, Maracay, Quito

!    nx_int=(/64,65,63,64,65,63,64,65,63,  47,48,46,47,48,46,47,48,46, 37,38,36,37,38,36,37,38,36,   56,57,55,56,57,55,56,57,55,   48,47,49,48,47,49,48,47,49,      82,83,81,82,83,81,82,83,81,        51,52,50,51,52,50,51,52,50,                          69,70,68,69,70,68,69,70,68        75,76,74,75,76,74,75,76,74         63,64,62,63,64,62,63,64,62,            73,74,72,73,74,72,73,74,72,     46,47,45,46,47,45,46,47,45                  44,45,43,44,45,43,44,45,43,               28,29,27,28,29,27,28,29,27,       48,49,47,48,49,47,48,49,47,       144,145,143,144,145,143,144,145,143,    91,92,90,91,92,90,91,92,90,        132,133,131,132,133,131,132,133,131,                   117,118,116,117,118,116,117,118,116,      136,137,135,136,137,135,136,137,135,    15,16,14,15,16,14,15,16,14,/)


!    ny_int=(/103,103,103,104,104,104,102,102,102, 120,120,120,121,121,121,119,119,119,  89,89,89,88,88,88,90,90,90,     173,173,173,174,174,174,172,172,172     166,166,166,167,167,167,165,165,165     139,139,139,138,138,138,140,140,140 100,100,100,101,101,101,99,99,99                      97,97,97,98,98,98,96,96,96           130,130,130,131,131,131,129,129,129            176,176,176,177,177,177,175,175,175             167,167,167,168,168,168,166,166,166     105,105,105,106,106,106,104,104,104       148,148,148,149,149,149,147,147,147            65,65,65,66,66,66,64,64,64         107,107,107,108,108,108,106,106,106       168,168,168,169,169,169,167,167,167     169,169,169,170,170,170,168,168,168        164,164,164,165,165,165,163,163,163       163,163,163,164,164,164,162,162,162    165,165,165,166,166,166,164,164,164     49,49,49,50,50,50,48,48,48/)
	

  nx_int=(/64,65,63,64,65,63,64,65,63,47,48,46,47,48,46,47,48,46,37,38,36,37, &
38,36,37,38,36,56,57,55,56,57,55,56,57,55,  48,47,49,48,47,49,48,47,49,82,83,81, &
82,83,81,82,83,81,51,52,50,51,52,50,51,52,50,69,70,68,69,70,68,69,70,68,75,76, &
74,75,76,74,75,76,74,63,64,62,63,64,62,63,64,62,73,74,72,73,74,72,73,74,72,46 &
,47,45,46,47,45,46,47,45,44,45,43,44,45,43,44,45,43,28,29,27,28,29,27,28,29,27 &
,48,49,47,48,49,47,48,49,47,144,145,143,144,145,143,144,145,143,91,92,90,91,92, &
90,91,92,90,132,133,131,132,133,131,132,133,131,117,118,116,117,118,116,117,118, &
116,136,137,135,136,137,135,136,137,135,15,16,14,15,16,14,15,16,14,128,129,130, &
128,129,130,128,129,130,68,69,70,68,69,70,68,69,70,125,126,127,125,126,127,125,126,127, &
125,126,127,125,126,127,125,126,127,128,129,130,128,129,130,128,129,130/)

    ny_int=(/103,103,103,104,104,104,102,102,102,120,120,120,121,121,121,119,119, &
119,89,89,89,88,88,88,90,90,90,173,173,173,174,174,174,172,172,172,166,166,166,167, &
167,167,165,165,165, 139,139,139,138,138,138,140,140,140, 100,100,100,101,101,101,99,99, &
99,97,97,97,98,98,98,96,96,96,130,130,130,131,131,131,129,129,129,176,176,176,177,177, &
177,175,175,175,167,167,167,168,168,168,166,166,166,105,105,105,106,106,106,104,104, &
104,148,148,148,149,149,149,147,147,147,65,65,65,66,66,66,64,64,64,107,107,107,108,108, &
108,106,106,106,168,168,168,169,169,169,167,167,167,169,169,169,170,170,170,168, &
168,168,164,164,164,165,165,165,163,163,163,163,163,163,164,164,164,162,162, &
162,165,165,165,166,166,166,164,164,164,49,49,49,50,50,50,48,48,48,163,163,163, &
164,164,164,162,162,162,155,155,155,156,156,156,157,157,157,163,163,163,164,164,164,162,162,162, &
166,166,166,167,167,167,165,165,165,166,166,166,167,167,167,165,165,165/)

  	!nx_int=(/64,65,63,64,65,63,64,65,63,47,48,46,47,48,46,47,48,46,37,38,36,37,38,36,37,38,36,56,57,55,56,57,55,56,57,55/)
    !ny_int=(/103,103,103,104,104,104,102,102,102,120,120,120,121,121,121,119,119,119,89,89,89,88,88,88,90,90,90,173,173,173,174,174,174,172,172,172/)
    
    
    agrup=(/1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4, &
4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,9,9,9,9,9, &
9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12, & 
12,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,16, &
16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,19,19,19, &
19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22, &
23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25, &
26,26,26,26,26,26,26,26,26/)
    
		

    if (j .eq. 1) then
		allocate(xb(maxval(agrup),time))
    endif



    do k=1,maxval(agrup)
    	do h=1,n_cells_int
  			if (agrup(h)==k) then
  				xb(k,:)=aux(nx_int(h),ny_int(h),1,1,:) 
  				exit
  			endif
  		enddo 
    enddo
    status= nf90_close(ncid) 
    if (j<10) then
    	open(61,file='./data/Ens_x0'//trim(FILE_NAME2)//'.dat',status='unknown')
    else
	  	open(61,file='./data/Ens_x'//trim(FILE_NAME2)//'.dat',status='unknown')
    endif
    write(61,*) xb(:,:)
    close(61)
  end do loop_Read_Ens
   write(*,*) 'Size Lon: ', lon, ' Size Lat: ', lat
! ============ Read information about NC files of Ensemble outputs Revisar formato de salida satelites, preguntar al nc por numero de observaciones============
	allocate(Y_xb(lon,lat,days_simulation))
	 
        write(*,*) 'Holi'
	loop_Read_Y_Ens: do j=1,Ens
		loop_days:  do i=0,days_simulation-1	   
		if (j < 10) then
			write (FILE_NAME2,'(i1)') j ! converting integer to string using a 'internal file'
			write (FILE_NAME3,'(i8)') i+initial_day        
			FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1800_xi0"//trim(FILE_NAME2)//'a.nc'
			INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
			if (file_exists .eqv. .false.) then
			  FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1700_xi0"//trim(FILE_NAME2)//'a.nc'
		          write(*,*) "No habia archivo a las 1800, se intenta con archivo a las 1700"
			  write(*,*) FILE_NAME
			  INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
			  if (file_exists .eqv. .false.) then
				write(*,*) "Tampoco hay a las 1700, pongase verga Yarce"
			  endif
			endif
		else            
			write (FILE_NAME2,'(i2)') j ! converting integer to string using a 'internal file'
			write (FILE_NAME3,'(i8)') i+initial_day 	   
			FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1800_xi"//trim(FILE_NAME2)//'a.nc'
			INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
			if (file_exists .eqv. .false.) then
			  FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1700_xi"//trim(FILE_NAME2)//'a.nc'
		          write(*,*) "No habia archivo a las 1800, se intenta con archivo a las 1700"
			  write(*,*) FILE_NAME
			  INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
			  if (file_exists .eqv. .false.) then
				write(*,*) "Tampoco hay a las 1700, pongase verga Yarce"
			  endif
			endif
		endif
		write(*,*) 'Writing Ensemble output '//trim(FILE_NAME2)//' Day '//trim(FILE_NAME3)
	  !Get for the dimensions ID 
		status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
		status = nf90_inq_dimid(ncid, "pixel", pixelid)
		!Get dimensions length
		status = nf90_inquire_dimension(ncid, pixelid, len = pixel_len) 
        write(*,*) 'Pixel Len ', pixel_len
		! Get the varid of the data variable, based on its name.
		status= nf90_inq_varid(ncid, "y", varid) 
        write(*,*) 'Status varid y',status
		allocate(aux_y(1,pixel_len))
		allocate(ilat(pixel_len))
        allocate(ilon(pixel_len))
		! Read the data.
aux_y=1
		status=nf90_get_var(ncid, varid, aux_y)
		write(*,*) FILE_NAME        
		write(*,*) 'Status Lectura ',status
        
		status= nf90_close(ncid)
		! Get ilonid and ilatid
		FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1800.nc"
        	INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		if (file_exists .eqv. .false.) then
		  FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1700.nc"
	          write(*,*) "No habia archivo a las 1800, se intenta con archivo a las 1700"
		  write(*,*) FILE_NAME
		  INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		  if (file_exists .eqv. .false.) then
			write(*,*) "Tampoco hay a las 1700, pongase verga Yarce"
		  endif
		endif
		status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
        status= nf90_inq_varid(ncid, "ilat", ilatid)
        status= nf90_inq_varid(ncid, "ilon", ilonid)
		! Read the data.
		status=nf90_get_var(ncid, ilatid, ilat)
        write(*,*) 'Status ilat ',status
        write(*,*) 'max ilat', maxval(ilat)
		status=nf90_get_var(ncid, ilonid, ilon)
		status= nf90_close(ncid)	  
		Y_xb=0.0
		do k=1,pixel_len
            if (aux_y(1,k) .ge. 0) then
			    Y_xb(ilon(k),ilat(k),i+1)=aux_y(1,k)
            endif
		enddo
		deallocate(aux_y)
		deallocate(ilat)
        deallocate(ilon)
		!=Check if exists other satelite file for the same day
		if (j < 10) then
			write (FILE_NAME2,'(i1)') j ! converting integer to string using a 'internal file'
			write (FILE_NAME3,'(i8)') i+initial_day        
			FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1900_xi0"//trim(FILE_NAME2)//'a.nc'
		else            
			write (FILE_NAME2,'(i2)') j ! converting integer to string using a 'internal file'
			write (FILE_NAME3,'(i8)') i+initial_day 	   
			FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1900_xi"//trim(FILE_NAME2)//'a.nc'
		endif
		INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		if (file_exists) then
			write(*,*) "Tambien se encontro archivo a las 1900"			
			write(*,*) FILE_NAME
			status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
			status = nf90_inq_dimid(ncid, "pixel", pixelid)
			!Get dimensions length
			status = nf90_inquire_dimension(ncid, pixelid, len = pixel_len)
			! Get the varid of the data variable, based on its name.
			status= nf90_inq_varid(ncid, "y", varid)
		    
			allocate(aux_y(1,pixel_len))
			allocate(ilat(pixel_len))
    	allocate(ilon(pixel_len))
			! Read the data.
			status=nf90_get_var(ncid, varid, aux_y)
			status= nf90_close(ncid)
			! Get ilonid and ilatid
			FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1900.nc"
			status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
		  status= nf90_inq_varid(ncid, "ilat", ilatid)
		  status= nf90_inq_varid(ncid, "ilon", ilonid)
			! Read the data.
			status=nf90_get_var(ncid, ilatid, ilat)
			status=nf90_get_var(ncid, ilonid, ilon)
			status= nf90_close(ncid)
			do k=1,pixel_len
              if (aux_y(1,k) .ge. 0) then
			      if (Y_xb(ilon(k),ilat(k),i+1) .NE. 0.0) then
					    Y_xb(ilon(k),ilat(k),i+1)=(aux_y(1,k)+ Y_xb(ilon(k),ilat(k),i))/2
			      else
					    Y_xb(ilon(k),ilat(k),i+1)=aux_y(1,k)
			      endif
               endif
			enddo
			deallocate(aux_y)
			deallocate(ilat)
    	    deallocate(ilon)
		endif
	enddo loop_days
	if (j<10) then
		open(61,file='./data/Y_Ens_x0'//trim(FILE_NAME2)//'.dat',status='unknown')
	else
		open(61,file='./data/Y_Ens_x'//trim(FILE_NAME2)//'.dat',status='unknown')
	endif
	write(61,*) Y_xb(:,:,:)
	close(61)
	write(*,*) 'max Y_xb',i+1,'=',(maxval(Y_xb))
enddo loop_Read_Y_Ens

! ============ Read information about NC files of observations ============ Modify with the template of the real observations
allocate(Y(lon,lat,days_simulation))
loop_days_y:  do i=0,days_simulation-1
	write (FILE_NAME3,'(i8)') i+initial_day 	   
	FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1800.nc"
	INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
	if (file_exists .eqv. .false.) then
	  FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1700.nc"
          write(*,*) "No habia archivo a las 1800, se intenta con archivo a las 1700"
	  write(*,*) FILE_NAME
	  INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
	  if (file_exists .eqv. .false.) then
		write(*,*) "Tampoco hay a las 1700, pongase verga Yarce"
	  endif
	endif

	write(*,*) 'Writing Observations for Day '//trim(FILE_NAME3)
	status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
	status = nf90_inq_dimid(ncid, "pixel", pixelid)
	!Get dimensions length
	status = nf90_inquire_dimension(ncid, pixelid, len = pixel_len) 
	! Get the varid of the data variable, based on its name.
	status= nf90_inq_varid(ncid, "yr", varid) 
	status= nf90_inq_varid(ncid, "ilat", ilatid)
	status= nf90_inq_varid(ncid, "ilon", ilonid)
	allocate(aux_y(1,pixel_len))
	allocate(ilat(pixel_len))
	allocate(ilon(pixel_len))
	! Read the data.
	status= nf90_get_var(ncid, varid, aux_y)
    write(*,*) 'Status Lectura ',status
   
	status= nf90_get_var(ncid, ilatid, ilat)
	status= nf90_get_var(ncid, ilonid, ilon)	
	status= nf90_close(ncid)
	Y=0.0
	do k=1,pixel_len
        if (aux_y(1,k) .ge. 0) then
		    Y(ilon(k),ilat(k),i+1)=aux_y(1,k) 
        endif
	enddo
	deallocate(aux_y)
	deallocate(ilat)
	deallocate(ilon)
	!=Check if exists other satelite file for the same day
	FILE_NAME=trim(path)//"/LE_"//trim(runid)//"_tropomi-no2_"//trim(FILE_NAME3)//"_1900.nc"
	INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
	if (file_exists) then
		status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
		status = nf90_inq_dimid(ncid, "pixel", pixelid)
		!Get dimensions length
		status = nf90_inquire_dimension(ncid, pixelid, len = pixel_len)
		! Get the varid of the data variable, based on its name.
		status= nf90_inq_varid(ncid, "yr", varid)
	  status= nf90_inq_varid(ncid, "ilat", ilatid)
		status= nf90_inq_varid(ncid, "ilon", ilonid) 
		allocate(aux_y(1,pixel_len))
		allocate(ilat(pixel_len))
  	allocate(ilon(pixel_len))
		! Read the data.
		status=nf90_get_var(ncid, varid, aux_y)
		status=nf90_get_var(ncid, ilatid, ilat)
		status=nf90_get_var(ncid, ilonid, ilon)
		status= nf90_close(ncid)
		do k=1,pixel_len
          if (aux_y(1,k) .ge. 0) then
		      if (Y(ilon(k),ilat(k),i+1) .NE. 0.0) then
				    Y(ilon(k),ilat(k),i+1)=(aux_y(1,k)+ Y(ilon(k),ilat(k),i+1))/2
			    else
				    Y(ilon(k),ilat(k),i+1)=aux_y(1,k)
			    endif
              endif
		enddo
		deallocate(aux_y)
		deallocate(ilat)
  	deallocate(ilon)
	endif
enddo loop_days_y
open(61,file='./data/Y.dat',status='unknown')
write(61,*) Y(:,:,:)
close(61)
 write(*,*) 'max Y observation = ', (maxval(Y))
deallocate(Y)



!===Write nx_int, ny_int for next scripts===
open(61,file='./data/cells_DC.dat',status='unknown')
write(61,*) n_cells_int
write(61,*) nx_int
write(61,*) ny_int
write(61,*) agrup
close(61)

  

end program read_le_ensemble_outputs
