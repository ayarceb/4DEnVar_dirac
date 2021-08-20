module module_enkf
use netcdf
contains


subroutine read_NC_Yk(Y_xk,dawindows)
	real*8,dimension(:,:),intent(inout)::Y_xk
	integer, intent(in)::dawindows
	real*8,dimension(:,:,:),allocatable::aux
	character (len=100)::FILE_NAME
	integer::lat,lon,level,time
  integer :: ncid, varid,status,latid,lonid,levelid,timeid

	FILE_NAME='/run/media/dirac/Datos/scratch/projects/4DEnVAR_inner/Xk/output/Y_xk.nc'
	status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
! Get the varid of the data variable, based on its name.
  status = nf90_inq_dimid(ncid, "latitude", latid)
  status = nf90_inq_dimid(ncid, "longitude", lonid)     
  status = nf90_inquire(ncid, unlimiteddimid = timeid)

  !Get dimensions length
  status = nf90_inquire_dimension(ncid, latid, len = lat)
  status = nf90_inquire_dimension(ncid, lonid, len = lon)
  status = nf90_inquire_dimension(ncid, timeid, len = time)
	
	status= nf90_inq_varid(ncid, "no2_column", varid)
  allocate(aux(lon,lat,time))

  ! Read the data.
  status=nf90_get_var(ncid, varid, aux(:,:,:))
  status= nf90_close(ncid)
  Y_xk(:,:)=reshape(aux(:,:,1:dawindows),(/lat*lon,dawindows/))
        
end subroutine read_NC_Yk

subroutine read_NC_Yk_Tropomi(Y_xk,dawindows,initial_day,st_ass,lon,lat)
	real*8,dimension(:,:),intent(inout)::Y_xk
	integer,intent(in)::lon,lat
	integer, intent(in)::dawindows,initial_day,st_ass
	real*8,dimension(:,:),allocatable::aux_y
	real*8,dimension(:,:,:),allocatable::Y(:,:,:)
	integer,dimension(:),allocatable::ilat(:),ilon(:)
	character (len=500)::FILE_NAME, FILE_NAME3
	integer::level,time,pixel,n_cells_int,i,initial_forecast,k
  integer :: ncid, varid,status,ilatid,ilonid,levelid,timeid,pixelid,pixel_len
  logical :: file_exists
  initial_forecast=initial_day+st_ass
  allocate(Y(lon,lat,dawindows))
  loop_days:  do i=0,dawindows-1
    write (FILE_NAME3,'(i8)') i+initial_forecast 	   
    !=====Modificar computador Andres====
		FILE_NAME="/run/media/dirac/Datos/scratch/projects/4DEnVAR_inner/Xk/output/LE_Xk_tropomi-no2_"//trim(FILE_NAME3)//"_1800.nc"
		INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		if (file_exists .eqv. .false.) then
		  FILE_NAME="/run/media/dirac/Datos/scratch/projects/4DEnVAR_inner/Xk/output/LE_Xk_tropomi-no2_"//trim(FILE_NAME3)//"_1700.nc"
                  write(*,*) "No habia archivo a las 1800, se intenta con archivo a las 1700"
		  write(*,*) FILE_NAME
		  INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		  if (file_exists .eqv. .false.) then
			write(*,*) "Tampoco hay a las 1700, pongase verga Yarce"
		  endif
		endif
		status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
	! Get the varid of the data variable, based on its name.
		status = nf90_inq_dimid(ncid, "pixel", pixelid)
		!Get dimensions length
		status = nf90_inquire_dimension(ncid, pixelid, len = pixel_len) 
		! Get the varid of the data variable, based on its name.
		status= nf90_inq_varid(ncid, "ilat", ilatid)
		status= nf90_inq_varid(ncid, "ilon", ilonid)
		allocate(aux_y(1,pixel_len))
		allocate(ilat(pixel_len))
		allocate(ilon(pixel_len))
		! Read the data.
		status= nf90_get_var(ncid, ilatid, ilat)
		status= nf90_get_var(ncid, ilonid, ilon)	
		status= nf90_close(ncid)
		FILE_NAME="/run/media/dirac/Datos/scratch/projects/4DEnVAR_inner/Xk/output/LE_Xk_tropomi-no2_"//trim(FILE_NAME3)//"_1800_xb.nc"
		write(*,*) FILE_NAME
		INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		if (file_exists .eqv. .false.) then
		  FILE_NAME="/run/media/dirac/Datos/scratch/projects/4DEnVAR_inner/Xk/output/LE_Xk_tropomi-no2_"//trim(FILE_NAME3)//"_1700_xb.nc"
                  write(*,*) "No habia archivo a las 1800, se intenta con archivo a las 1700"
		  write(*,*) FILE_NAME
		  INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		  if (file_exists .eqv. .false.) then
			write(*,*) "Tampoco hay a las 1700, pongase verga Yarce"
		  endif
		endif
		status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
		status= nf90_inq_varid(ncid, "y", varid) 
		status= nf90_get_var(ncid, varid, aux_y)
	write(*,*) 'Status read y function = ', status
		Y=0.0
		do k=1,pixel_len
            if (aux_y(1,k) .ge. 0) then
                            Y(ilon(k),ilat(k),i+1)=aux_y(1,k) 
            endif
		enddo
        write(*,*)'Max value Y 1800 = ', maxval(Y)
		deallocate(aux_y)
		deallocate(ilat)
		deallocate(ilon)
		!=Check if exists other satelite file for the same day
		FILE_NAME="/run/media/dirac/Datos/scratch/projects/4DEnVAR_inner/Xk/output/LE_Xk_tropomi-no2_"//trim(FILE_NAME3)//"_1900.nc"
		INQUIRE(FILE=FILE_NAME, EXIST=file_exists)
		if (file_exists) then
            write(*,*) FILE_NAME
			status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
			status = nf90_inq_dimid(ncid, "pixel", pixelid)
			!Get dimensions length
			status = nf90_inquire_dimension(ncid, pixelid, len = pixel_len)
			! Get the varid of the data variable, based on its name.
			status= nf90_inq_varid(ncid, "ilat", ilatid)
			status= nf90_inq_varid(ncid, "ilon", ilonid) 
			allocate(aux_y(1,pixel_len))
			allocate(ilat(pixel_len))
			allocate(ilon(pixel_len))
			! Read the data.
			status=nf90_get_var(ncid, ilatid, ilat)
			status=nf90_get_var(ncid, ilonid, ilon)
			status= nf90_close(ncid)
			FILE_NAME="/run/media/dirac/Datos/scratch/projects//4DEnVAR_inner/Xk/output/LE_Xk_tropomi-no2_"//trim(FILE_NAME3)//"_1900_xb.nc"
			status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
			status= nf90_inq_varid(ncid, "y", varid)
			status=nf90_get_var(ncid, varid, aux_y)
write(*,*) "Read Tercer"
			do k=1,pixel_len
                if (aux_y(1,k) .ge. 0) then
				    if (Y(ilon(k),ilat(k),i+1) .NE. 0.0) then
					    Y(ilon(k),ilat(k),i+1)=(aux_y(1,k)+ Y(ilon(k),ilat(k),i+1))/2
				    else
					    Y(ilon(k),ilat(k),i+1)=aux_y(1,k)
				    endif
                endif
			enddo
write(*,*)'Max value Y 1900 = ', maxval(Y)
			deallocate(aux_y)
			deallocate(ilat)
			deallocate(ilon)
		endif
		Y_xk(:,:)=reshape(Y(:,:,:),(/lat*lon,dawindows/))
write(*,*)'Max value Y_xk = ', maxval(Y_xk)
	enddo loop_days
        
end subroutine read_NC_Yk_Tropomi



subroutine create_DC_NC(Xk,n,n_cells_int,agrup,lat,lon,time,nx_int,ny_int,initial_day,st_ass)
!===Create a NC file with DC for LE forecast simulation====
   implicit none
   integer,intent(in)  :: lat,lon,time,n_cells_int,n
   real*8,dimension(:), intent(inout)::Xk(n)
   integer,dimension(:), intent(inout):: nx_int(n_cells_int),ny_int(n_cells_int),agrup(n_cells_int)
   real*8,dimension(:,:)::DC(lon,lat)
   real*8 ::DCRnd
   real*8,dimension(:)::latitude(lat),longitude(lon),time_var(time)
   real*8,dimension(:,:,:,:,:)::DC_nc(lon,lat,1,1,time)
   integer::ncid,lat_dimid,lon_dimid,time_dimid,noise_dimid,hist_dimid,name_dimid
   integer::lat_varid,lon_varid,time_varid,dc_varid,noise_varid	,status
   integer::i,initial_day,st_ass
   character (len=100)::FILE_NAME,Path,FILE_NAME3
    
	open(61,file='./temp/Ensembles.in',status='old')
	read(61,"(A200)") Path
  close(61)
  FILE_NAME=trim(path)//'/Ens_x01.nc'

  status = nf90_open(path = FILE_NAME, mode = nf90_nowrite, ncid = ncid)
  ! Get the varid of the data variable, based on its name.
  status= nf90_inq_varid(ncid, "latitude", lat_varid)
	status= nf90_inq_varid(ncid, "longitude", lon_varid)
	status= nf90_inq_varid(ncid, "time", time_varid)
  !Read the data.
  status=nf90_get_var(ncid, lat_varid, latitude)
	status=nf90_get_var(ncid, lon_varid, longitude)
	status=nf90_get_var(ncid, time_varid, time_var)
  status =nf90_close(ncid)

!===Create new DC file initial====

	write (FILE_NAME3,'(i8)') st_ass+initial_day
  !status = nf90_create('./data/DC_inner_dc_'//trim(FILE_NAME3)//'_xa.nc',NF90_CLOBBER, ncid)  -11130
  status = nf90_create('./data/DC_inner_dc_  -11130_xa.nc',NF90_CLOBBER, ncid)
!=== Define dimensions===
  status = nf90_def_dim(ncid, "latitude", lat, lat_dimid)
  status = nf90_def_dim(ncid, "longitude", lon, lon_dimid)
  status = nf90_def_dim(ncid, "noise", 1, noise_dimid)
  status = nf90_def_dim(ncid, "hist", 1, hist_dimid)
  status = nf90_def_dim(ncid, "time", time, time_dimid)
  status = nf90_def_dim(ncid, "namelen", 32, name_dimid)
!=== Define variables===
  status = nf90_def_var(ncid, "latitude", NF90_REAL8, (/ lat_dimid/), lat_varid)
  status = nf90_put_att(ncid, lat_varid, "standard_name", 'latitude')
  status = nf90_put_att(ncid, lat_varid, "units", 'degrees_north')
  status = nf90_put_att(ncid, lat_varid, "_CoordinateAxisType", 'Lat')

  status = nf90_def_var(ncid, "longitude", NF90_REAL8, (/lon_dimid/), lon_varid)
  status = nf90_put_att(ncid, lon_varid, "standard_name", 'longitude')
  status = nf90_put_att(ncid, lon_varid, "units", 'degrees_east')
  status = nf90_put_att(ncid, lon_varid, "_CoordinateAxisType", 'Lon')	

  status = nf90_def_var(ncid, "time", NF90_REAL8, (/ time_dimid/), time_varid)
  status = nf90_put_att(ncid, time_varid, "standard_name", 'time')
  status = nf90_put_att(ncid, time_varid, "long_name", 'time')
  status = nf90_put_att(ncid, time_varid, "units", 'seconds since 2019-01-01 00:00:00 UTC')
  status = nf90_put_att(ncid, time_varid, "calendar", 'standard')
  status = nf90_put_att(ncid, time_varid, "_CoordinateAxisType", 'Time')	

  status = nf90_def_var(ncid, "noise_name", NF90_CHAR, (/name_dimid, noise_dimid/), noise_varid)

  status = nf90_def_var(ncid, "dc", NF90_REAL8, (/ lon_dimid,lat_dimid,noise_dimid,hist_dimid,time_dimid/), dc_varid)
   
  status = nf90_enddef(ncid)
!=== Write Data===
  status =nf90_put_var(ncid, lat_varid, latitude)
  status =nf90_put_var(ncid, lon_varid, longitude)
  status =nf90_put_var(ncid, time_varid, time_var)
  DC=1
  CALL RANDOM_NUMBER(DCrnd)
		do i=1,n_cells_int
			if(i>1)then
				if(agrup(i).ne.agrup(i-1)) then
					CALL RANDOM_NUMBER(DCrnd)
				end if			
			end if
			
			write(*,*)'Dcnrd = ', DCrnd
			DC(nx_int(i),ny_int(i))=max(Xk(agrup(i)),0.5*DCrnd)
		        
		


		enddo
  do i=1,time
		DC_nc(:,:,1,1,i)=DC	
  enddo
	status =nf90_put_var(ncid, dc_varid, DC_nc)
	status =nf90_put_var(ncid, noise_varid, 'NOx_emis                        ')
	status =nf90_close(ncid)

end subroutine create_DC_NC
  
subroutine ensmean(Xen,Xmean)
    !===Calculate de Mean of a Matrix====
   implicit none
   integer :: nrens,j
   real*8,dimension(:,:), intent(in)  :: Xen
   real*8,dimension(size(Xen,1)),intent(out) :: Xmean
   nrens=size(xen,2)
   write(*,*) 'size xen', shape(Xen)
   Xmean(:)=Xen(:,1)
   do j=2,nrens
      Xmean(:)=Xmean(:)+Xen(:,j)
   enddo
   Xmean=(1.0/real(nrens))*Xmean

end subroutine ensmean


subroutine anomalies(Xen,Xmean,dXen)

   !===Calculate de Anomalies Matrix====
   implicit none
   integer :: nrens,j
   real*8,dimension(:,:), intent(in)  :: Xen
   real*8,dimension(:),intent(in)::Xmean
   real*8,dimension(size(Xen,1),size(Xen,2)),intent(out)::dXen
   nrens=size(Xen,2)
   do j=1,nrens
      dXen(:,j)=(Xen(:,j)-Xmean(:))
   end do
   dXen=(1.0/sqrt(real(nrens-1)))*dXen


end subroutine anomalies

subroutine Analysis_EnKF(Xben,P,H,R,sigmaobs,y,Xaen)
  use module_matrix
  implicit none
  real*8,dimension(:,:),intent(in)::Xben,P,H,R
  real*8,dimension(:),intent(in)::y
  real*8,intent(in)::sigmaobs
  real*8,dimension(size(Xben,1),size(Xben,2)),intent(out)::Xaen
  real*8,dimension(size(Xben,1))::z
  real*8,dimension(size(y,1))::v,d
  real*8,dimension(size(Xben,1),size(H,1))::Kalman,PHT
  real*8,dimension(size(H,1),size(H,1))::A
  integer ::k

  PHT=prod(P,transpose(H))
  A=prod(H,PHT)
  A=A+R
  A=inv(A)
  Kalman=prod(PHT,A)

  do k=1,size(Xben,2)
     call RANDOM_NUMBER(v(:))
     z(:)= prod_matvec(H(:,:), Xben(:,k))
      d(:)=y(:)-z(:)+  (sigmaobs)*v(:)
     Xaen(:,k)=Xben(:,k)+prod_matvec(Kalman,d)
  end do

  
end subroutine Analysis_EnKF



  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      !print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check




end module module_EnKF
