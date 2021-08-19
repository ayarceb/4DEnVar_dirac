module module_EnKF
use netcdf
contains


subroutine create_DC_NC(DC,lat,lon,time)
!===Create a NC file with DC for LE forecast simulation====
   implicit none
   integer,intent(in)  :: lat,lon,time
   real*8,dimension(:,:), intent(inout)::DC(lat,lon)
   real*8,dimension(:)::latitude(lat),longitude(lon),time_var(time)
   real*8,dimension(:,:,:,:,:)::DC_nc(lat,lon,1,1,time)
   integer::ncid,lat_dimid,lon_dimid,time_dimid,noise_dimid,hist_dimid
   integer::lat_varid,lon_varid,time_varid,dc_varid,noise_varid	,status
   character (len=100)::FILE_NAME,Path
    !DC=3 !ojo borrar
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


   status = nf90_create('./data/Putona.nc',NF90_CLOBBER, ncid)
!=== Define dimensions===
   status = nf90_def_dim(ncid, "latitude", lat, lat_dimid)
   status = nf90_def_dim(ncid, "longitude", lon, lon_dimid)
   status = nf90_def_dim(ncid, "noise", 1, noise_dimid)
   status = nf90_def_dim(ncid, "hist", 1, hist_dimid)
   status = nf90_def_dim(ncid, "time", time, time_dimid)
!=== Define variables===
   status = nf90_def_var(ncid, "latitude", NF90_INT, (/ lat_dimid/), lat_varid)
   status = nf90_def_var(ncid, "longitude", NF90_INT, (/ lon_dimid/), lon_varid)
   status = nf90_def_var(ncid, "time", NF90_INT, (/ time_dimid/), time_varid)
   status = nf90_def_var(ncid, "noise_name", NF90_INT, (/ noise_dimid/), noise_varid)
   status = nf90_def_var(ncid, "dc", NF90_INT, (/ lon_dimid,lat_dimid,noise_dimid,hist_dimid,time_dimid/), lon_varid)

   status = nf90_enddef(ncid)

!=== Write Data===
    status =nf90_put_var(ncid, lat_varid, latitude)
    status =nf90_put_var(ncid, lon_varid, longitude)
	status =nf90_put_var(ncid, time_varid, time_var)
	DC_nc=reshape(DC,(/lat,lon,1,1,time/))	
	status =nf90_put_var(ncid, dc_varid, DC_nc)

	status =nf90_close(ncid)

end subroutine create_DC_NC
  
subroutine ensmean(Xen,Xmean)
    !===Calculate de Mean of a Matrix====
   implicit none
   integer :: nrens,j
   real*8,dimension(:,:), intent(in)  :: Xen
   real*8,dimension(size(Xen,2)),intent(out) :: Xmean
   nrens=size(Xen,2)
   
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
