module MODULE_MATRIX

contains

!===Invertir una matriz ====
 function inv(A) result(Ainv)
   integer,parameter::dp=8
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv


function   prod(A,B) result(C)
  external  DGEMM
   double precision, dimension(:,:), intent(in) :: A,B
   double precision,dimension(size(A,1),size(B,2))::C
   double precision::alpha=1.0,beta=0.0
   integer::n,m,k,amax,bmax,cmax
   m=size(A,DIM=1)
   if (size(A,2) /= size(B,1)) then
      stop 'Inner matrix dimensions must agree.'
   end if
   k=size(A,DIM=2)
   n=size(B,DIM=2)
   amax=size(A,DIM=1)
   bmax=size(B,DIM=1)
   cmax=size(C,DIM=1)

   call DGEMM('n','n',m,n,k,alpha,A,amax,B,bmax,beta,C,cmax)
 end function prod



 function   prod_matvec(A,B) result(C)
  external  DGEMM
  double precision, dimension(:,:), intent(in) :: A
 double precision, dimension(:),intent(in):: B
   double precision,dimension(size(A,1))::C
   double precision::alpha=1.0,beta=0.0
   integer::n,m,k,amax,bmax,cmax
   m=size(A,DIM=1)
   if (size(A,2) /= size(B,1)) then
      stop 'Inner matrix dimensions must agree.'
   end if
   k=size(A,DIM=2)
   n=1
   amax=size(A,DIM=1)
   bmax=size(B,DIM=1)
   cmax=size(C,DIM=1)

   call DGEMM('n','n',m,n,k,alpha,A,amax,B,bmax,beta,C,cmax)
end function prod_matvec



subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r8mat_print_some








subroutine print_matrix ( m, n, a, title )
!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end subroutine print_matrix


subroutine randperm(N,m,p,pout)
    integer, intent(in)                :: N
    integer, dimension(:), intent(inout) :: p
    integer, dimension(:), intent(out) :: pout
    integer :: j, k
    real :: u
    p = 0
    do j=1,N
      call random_number(u)
      k = floor(j*u) + 1
      p(j) = p(k)
      p(k) = j
   end do
   pout=p(1:m)
   
end subroutine randperm





end module MODULE_MATRIX
