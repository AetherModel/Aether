
module memories

  implicit none

  integer saved_array(10)

end module memories

subroutine convert_int_to_string(inArray, outString, length)
  implicit none
  integer, intent(in) :: length
  integer, dimension(length), intent(in) :: inArray
  character(length), intent(out) :: outString

  integer :: i
  do i = 1, length
    if (inArray(i) > 1) then
      write(*,*) 'for : ', i, inArray(i)
      outString(i:i) = achar(inArray(i))
    else
      outString(i:i) = achar(0)
    endif
  enddo

end subroutine convert_int_to_string

subroutine test_passing_string(intArray) bind(C, name = 'test_passing_string')
  implicit none
  integer, parameter :: iLength_ = 100
  integer :: intArray(iLength_)
  character(iLength_) :: charArray
  integer :: i

  call convert_int_to_string(intArray, charArray, iLength_)

!  do i = 1, iLength_
!    if (intArray(i) > 1) then
!      write(*,*) 'for : ', i, intArray(i)
!      charArray(i:i) = achar(intArray(i))
!    else
!      charArray(i:i) = ' '
!    endif
!  enddo

  write(*,*) "string : ", charArray
end subroutine test_passing_string

subroutine print_hi() bind(C, name = 'print_hi')
  implicit none
  write(*,*) "Hello from Fortran."
end subroutine print_hi

subroutine print_double(i, x, y) bind(C, name = 'print_double')
  implicit none
  integer :: i
  real :: x
  real :: y
  y = 2.0 * x
  write(*,*) "Integer in Fortran: ", i
  write(*,*) "Real in Fortran: ", x
end subroutine print_double

subroutine pass_arrays(array, twotimes) bind(C, name = 'pass_arrays')

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

subroutine get_array(array_back) bind(C, name = 'get_array')

  use memories
  
  implicit none

  integer :: array_back(10)
  integer :: i

  do i = 1, 10
     array_back(i) = saved_array(i)
  enddo
  
end subroutine get_array
