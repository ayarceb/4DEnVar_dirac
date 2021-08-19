  !==Program to run the Real State Of LE===
  !==Andres Yarce, Santiago Lopez-Restrepo==

program run_lotos_euros_real
  use module_EnKF
  use module_matrix
  use netcdf
  implicit none

  integer::tsim
  integer,parameter::dx=3600,dy=1800
  real*8,dimension(:,:)::emision(dx,dy)
  real*8::pertur
  integer::j
  character (len =200 )::FILE_NAME,Path
  integer :: ncid, varid,status
  
  open(61,file='./Parameters/Model_Parameters.in',status='old')
  read(61,*)pertur
  close(61)

  open(61,file='./Parameters/DA_Parameters.in',status='old')
  read(61,*)  !Number of states
  read(61,*)  !Number of Observations
  read(61,*)  !Number of Ensembles
  read(61,*) tsim
  close(61)
  emision=1

  
  open(61,file='./temp/files.dat',status='old')
  read(61,"(A200)") Path
  loop_modifie_all_emissions_file: do j=1,20
     read(61,*) FILE_NAME
     FILE_NAME=trim(path)//"/"//trim(FILE_NAME)
    
    !write(*,*) FILE_NAME
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
    ! the file.

    status= nf90_open(path = FILE_NAME, mode = nf90_write, ncid = ncid)


    ! Get the varid of the data variable, based on its name.
    status= nf90_inq_varid(ncid, "emi_nox", varid)


    ! Read the data.
    status=nf90_get_var(ncid, varid, emision)


    ! Write the data
    status=nf90_put_var(ncid, varid, emision*pertur)
    emision=emision*pertur

    ! Close the file, freeing all resources.
    status= nf90_close(ncid)
  end do Loop_Modifie_All_Emissions_File

  close(61)




  
  
  
  

end program run_lotos_euros_real
