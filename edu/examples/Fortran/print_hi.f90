
module memories

  implicit none

  integer saved_array(10)

end module memories


subroutine print_hi() bind(C)
  implicit none
  write(*,*) "Hello from Fortran."
end subroutine print_hi

subroutine print_double(i, x, y) bind(C)
  implicit none
  integer :: i
  real :: x
  real :: y
  y = 2.0 * x
  write(*,*) "Integer in Fortran: ", i
  write(*,*) "Real in Fortran: ", x
end subroutine print_double

subroutine pass_arrays(array, twotimes) bind(C)

  use memories
  
  implicit none
  integer :: array(10)
  integer :: twotimes(10)
  integer :: i
  do i = 1, 10
     write(*,*) "Value ", i, ": ", array(i)
     twotimes(i) = array(i) * 2
     saved_array(i) = twotimes(i)
  enddo
end subroutine pass_arrays

subroutine get_array(array_back) bind(C)

  use memories
  
  implicit none

  integer :: array_back(10)
  integer :: i

  do i = 1, 10
     array_back(i) = saved_array(i)
  enddo
  
end subroutine get_array
