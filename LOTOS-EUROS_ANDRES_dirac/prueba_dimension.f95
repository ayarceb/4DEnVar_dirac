program prueba_dimension
implicit none
real*8,dimension(:,:,:,:,:)::xb(2,3,1,1,4)
real*8,dimension(:,:)::xa(2*3,4)
real*8,dimension(:)::aux_xb(2*3*4)
integer::j,i,t

do j=1,2
	do i=1,3
		do t=1,4
					xb(j,i,1,1,t)=i*j*t
		enddo
	enddo
enddo

open(61,file='Prueba_dimension.in',status='unknown')
write(61,*) xb(:,:,1,1,:)
 close(61)

open(61,file='Prueba_dimension.in',status='unknown')
read(61,*) aux_xb(:)
 close(61)

xa=reshape(aux_xb,(/2*3,4/))

do t=1,4
write(*,*) "Nuevo tiempo"
		write(*,*)xa(:,t)
enddo




end program
