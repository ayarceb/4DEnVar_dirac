program prueba_nc
use netcdf
use module_EnKF

implicit none
	character (len=100)::FILE_NAME,Path,FILE_NAME2  
	integer::ncid,lat_dimid,lon_dimid,time_dimid,noise_dimid,hist_dimid
	integer::lat_varid,lon_varid,time_varid,dc_varid,noise_varid	
    real*8,dimension(:)::DC(110*124)
    integer::time,lat,lon
    lon=110
	lat=124
	time=49
    DC=2.495
	call create_DC_NC(DC,lat,lon,time)
end program prueba_nc
