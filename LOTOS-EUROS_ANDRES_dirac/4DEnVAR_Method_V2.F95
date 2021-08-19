program LE_4DEnVAR_Method
  !USE modulo_distribucion_normal
  use module_EnKF
  use module_matrix
  implicit none

  integer::n,m,tsim,Nen,dawindows,inner,fileidreal=1000
  integer::st_ass,n_cells_int,days_simulation
  integer::lat,lon,time !Dimensions	
  real*8,parameter::F=8
  real*8::sigmaobs
  integer::k,i,j,l,fileid     ! i for states, j for time, k for ensembles,l for inner steps
  real*8::realaux
  integer::intaux,flaq,initial_day
  real*8,parameter :: dt=0.01
  real*8::tolerance,sparse
  real*8,allocatable,dimension(:,:,:):: xa,xb,Yb,xk
  real*8,allocatable,dimension(:,:)::dXen,Robs,inR,Y_xk
  real*8,allocatable,dimension(:,:)::increment,inno,sum1,ident,Y
  real*8,allocatable,dimension(:)::sum2,norma
  real*8,allocatable,dimension(:,:,:)::dY
  real*8,allocatable,dimension(:)::aux_vect,aux_vect2,xbmean,xamean,ymean,aux_read,aux_read_y
  integer,allocatable,dimension(:) :: nx_int(:),ny_int(:)
  character (len=500)::FILE_NAME,Path,FILE_NAME2
  

  Character(len=20)::filename

!====Reading parameters from .in file====
open(61,file='./data/dimensions.dat',status='old')
read(61,*) lat
read(61,*) lon
read(61,*) time
close(61)


open(61,file='./temp/Ensembles.in',status='old')
read(61,"(A200)")
read(61,*) Nen
close(61)
open(61,file='./data/cells_DC.dat',status='old')
read(61,*) n_cells_int
allocate(nx_int(n_cells_int))
allocate(ny_int(n_cells_int))
read(61,*) nx_int
read(61,*) ny_int
close(61)
n=n_cells_int
m=lat*lon
tsim=time

open(61,file='./DATA_4DEnVAR/startdate.in',status='old')
read(61,*) initial_day
read(61,*) days_simulation
close(61)

!====Reading parameters from .in file====
open(61,file='./DATA_4DEnVAR/parameters.in',status='old')
read(61,*) dawindows
read(61,*) st_ass
read(61,*) inner
read(61,*) tolerance
read(61,*) sparse
read(61,*) sigmaobs
close(61)

!===Matriz Initialization===
allocate(xb(n,Nen,tsim))
allocate(xa(n,Nen,tsim))
allocate(Yb(m,Nen,days_simulation))
allocate(Y(m,days_simulation))
allocate(dXen(n,Nen))
allocate(dY(m,Nen,days_simulation))
allocate(Xamean(n))
allocate(Xbmean(n))
allocate(ymean(m))
allocate(Robs(m,m))
allocate(inR(m,m))
allocate(increment(Nen,inner))
allocate(inno(m,days_simulation))
allocate(ident(Nen,Nen))
allocate(norma(inner))
allocate(Y_xk(m,dawindows))
Robs=0
inR=0
do i=1,m
    inR(i,i)=(1.0/(sigmaobs**2))
enddo


ident=0
do i=1,Nen
	ident(i,i)=1
enddo
write(*,*) "====Starting Analysis 4DEnVAR===="
!===Read Ensembles, Ensembles Outputs and Observations ===
allocate(aux_read(n*time))
allocate(aux_read_y(m*days_simulation))



loop_Read_Ens: do k=1,Nen
   if (k < 10) then            
    write (FILE_NAME2,'(i1)') k ! converting integer to string using a 'internal file'
    FILE_NAME='./data/Ens_x0'//trim(FILE_NAME2)//'.dat'
   else            
    write (FILE_NAME2,'(i2)') k ! converting integer to string using a 'internal file'	   
    FILE_NAME='./data/Ens_x'//trim(FILE_NAME2)//'.dat'
   endif
   open(61,file=FILE_NAME)
   read(61,*) aux_read(:)

   xb(:,k,:)=reshape(aux_read,(/n,tsim/))

   close(61)	
   if (k < 10) then            
    write (FILE_NAME2,'(i1)') k ! converting integer to string using a 'internal file'
    FILE_NAME='./data/Y_Ens_x0'//trim(FILE_NAME2)//'.dat'
   else            
    write (FILE_NAME2,'(i2)') k ! converting integer to string using a 'internal file'	   
    FILE_NAME='./data/Y_Ens_x'//trim(FILE_NAME2)//'.dat'
   endif
   open(61,file=FILE_NAME)
   read(61,*) aux_read_y(:)
   close(61)
   Yb(:,k,:)=reshape(aux_read_y,(/m,days_simulation/))
    write(*,*) 'Max Y ensemble', k, ' ==', maxval(Yb(:,k,:))
   	
 enddo loop_Read_Ens

 FILE_NAME='./data/Y.dat'
 open(61,file=FILE_NAME)
 read(61,*) aux_read_y(:)
 close(61)
 Y(:,:)=reshape(aux_read_y,(/m,days_simulation/))
write(*,*) 'Max Y observation', k, ' == ', maxval(Y(:,:))

!=== Ensemble initialization===

loop_time_Y:do j=1,days_simulation    
   call ensmean(Yb(:,:,j),ymean(:)) !Calculate the  Ensemble Output mean, from module_EnKF
   call anomalies(Yb(:,:,j),ymean(:),dY(:,:,j)) !Calculate anomalies matrix dy, from module_EnKF
end do loop_time_Y


call ensmean(xb(:,:,st_ass*24 +1),xbmean(:)) !Calculate the Ensemble mean, from module_EnKF
call anomalies(xb(:,:,st_ass*24 +1),xbmean(:),dXen) !Calculate anomalies matrix dX, from module_EnKF



!=====Analysis Step===
!===Inner Loop===
 CALL RANDOM_NUMBER(increment(:,1))
 allocate(xk(n,tsim,inner))
 l=1
write (FILE_NAME2,'(i1)') l
FILE_NAME='./data/increment_'//trim(FILE_NAME2)//'.dat'
open(61,file=FILE_NAME)
write(61,*) increment(:,1) 
close(61)
FILE_NAME='./data/Dxen_'//trim(FILE_NAME2)//'.dat'
open(61,file=FILE_NAME)
write(61,*) dXen 
close(61)

!=== Calculate DC_k initial====
xk(:,st_ass+1,1)=xbmean+prod_matvec(dXen,increment(:,1))
write (FILE_NAME2,'(i2)') l
FILE_NAME='./data/Xk_'//trim(FILE_NAME2)//'.dat'
open(61,file=FILE_NAME)
write(61,*) xk(:,st_ass+1,l) 
close(61)
FILE_NAME='./data/Xbmean_'//trim(FILE_NAME2)//'.dat'
open(61,file=FILE_NAME)
write(61,*) xbmean 
close(61)
call create_DC_NC(xk(:,st_ass+1,1),n_cells_int,lat,lon,time,nx_int,ny_int,initial_day,st_ass)
!=== Run LE using that DC_k Value Ojo cambiar Dirac===
CALL execute_command_line('./aux_launcher_inner')
write(*,*) '=== Merging Yxk==='
!=== Read DC_k and Column of DC_k Revisar esta parte===
!CALL execute_command_line('./aux_merge_Y_xk')
write(*,*) '=== Reading Yxk==='
Call read_NC_Yk_Tropomi(Y_xk,dawindows,initial_day,st_ass,lon,lat)
write(*,*) "Max Value Y_xk outside", maxval(Y_xk)
write(*,*) "Size Y_xk", shape(Y_xk)
write(*,*) '=== Yxk Readed==='

loop_time_analysis: do j=(st_ass+2),st_ass+1+dawindows
  write(*,*) 'j-(st_ass+1)=', j-(st_ass+1)
  inno(:,j-(st_ass+1))=Y(:,j-1)-Y_xk(:,j-(st_ass+1))
  write(*,*) 'Max value Y(:,', j-1,') =', maxval(Y(:,j-1))
  write(*,*) 'Max value Y(:,', j-(st_ass+1),') =', maxval(Y_xk(:,j-(st_ass+1)))
end do loop_time_analysis

write(*,*) "Max val inno inner", l, " = ", maxval(inno)
allocate(sum1(Nen,Nen))
allocate(sum2(Nen))
l=2
norma(1)=norm2(prod_matvec(dXen,increment(:,1)))
write(*,*) "Value of Norm iteration ", l, "= ",   norma(1) 
    FILE_NAME='./data/Norma_1.dat'
    open(61,file=FILE_NAME)
    write(61,*) norma(1)
    close(61)
norma(1)=tolerance
flaq=0
allocate(aux_vect(m))
loop_inner:do while (l .LE. inner) !.AND. (norma(l-1) .GE. tolerance))
  sum1=0
  sum2=0
  do i=st_ass+2,st_ass+1+dawindows
     sum1=sum1+prod(prod(transpose(dY(:,:,i-1)),inR),dY(:,:,i-1))
     sum2=sum2+ prod_matvec(prod(transpose(dY(:,:,i-1)),inR),inno(:,i-(st_ass+1)))
  end do
  FILE_NAME='./data/Sum1_'//trim(FILE_NAME2)//'.dat'
    open(61,file=FILE_NAME)
    write(61,*) sum1
    close(61)
    FILE_NAME='./data/Sum2_'//trim(FILE_NAME2)//'.dat'
    open(61,file=FILE_NAME)
    write(61,*) sum2
    close(61)
   

  write(*,*) '===Calculate Increment==='
  increment(:,l)=prod_matvec(inv(ident+sum1),sum2)
  write (FILE_NAME2,'(i2)') l
	FILE_NAME='./data/increment_'//trim(FILE_NAME2)//'.dat'
	open(61,file=FILE_NAME)
	write(61,*) increment(:,l) 
	close(61)
  !=== Incremental===
  xk(:,st_ass+1,l)=xk(:,st_ass+1,l-1)+0.1*prod_matvec(dXen,increment(:,l))
  write (FILE_NAME2,'(i2)') l
	FILE_NAME='./data/Xk_'//trim(FILE_NAME2)//'.dat'
	open(61,file=FILE_NAME)
	write(61,*) xk(:,st_ass+1,l) 
	close(61)
	!=== Write DC_k output ====
	Call create_DC_NC(xk(:,st_ass+1,l),n_cells_int,lat,lon,time,nx_int,ny_int,initial_day,st_ass)
	!=== Run LE using that DC_k Value Ojo cambiar Dirac===
	CALL execute_command_line('./aux_launcher_inner')
	write(*,*) '=== Merging Yxk, iteration ', l 
	!=== Read DC_k and Column of DC_k===
	!CALL execute_command_line('./aux_merge_Y_xk')
	write(*,*) '=== Reading Yxk iteration', l 
	Call read_NC_Yk_Tropomi(Y_xk,dawindows,initial_day,st_ass,lon,lat)
	write(*,*) '=== Yxk Readed iteration', l 
     loop_time_analysis_inner: do j=(st_ass+2),st_ass+1+dawindows
          write(*,*) 'j-(st_ass+1)=', j-(st_ass+1)
          inno(:,j-(st_ass+1))=Y(:,j-1)-Y_xk(:,j-(st_ass+1))
          write(*,*) 'Max value Y(:,', j-1,') =', maxval(Y(:,j-1))
          write(*,*) 'Max value Y(:,', j-(st_ass+1),') =', maxval(Y_xk(:,j-(st_ass+1)))
    end do loop_time_analysis_inner
write(*,*) "Max val inno inner", l, " = ", maxval(inno)
  
    FILE_NAME='./data/innovation_'//trim(FILE_NAME2)//'.dat'
	open(61,file=FILE_NAME)
	write(61,*) inno(:,l) 
	close(61)
	FILE_NAME='./data/Dxen_'//trim(FILE_NAME2)//'.dat'
	open(61,file=FILE_NAME)
	write(61,*) dXen 
	close(61)
  norma(l)=norm2(prod_matvec(dXen,increment(:,l)))
  FILE_NAME='./data/Norma_'//trim(FILE_NAME2)//'.dat'
    open(61,file=FILE_NAME)
    write(61,*) norma(l)
    close(61)
  write(*,*) "Value of Norm iteration ", l, "= ",   norma(l) 
  if (norma(l) .LE. norma(l-1)) then
     flaq=flaq+1
  end if     
  l=l+1
end do loop_inner
 


   
end program LE_4DEnVAR_Method
