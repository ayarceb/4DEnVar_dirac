program prueba_lapack
  !USE modulo_distribucion_normal
  use module_EnKF
  use module_matrix
  implicit none

  real*8,dimension(:,:)::x(2,2)
  real*8,dimension(:)::a(2)
  
  x=2
  write(*,*) x  
a=0
  
  write(*,*) inv(x)
   write(*,*) prod_matvec(x,a)
   write(*,*) prod(x,x)
  call ensmean(x,a)   


end program prueba_lapack
