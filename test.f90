program test
implicit none

real, dimension (10) :: test_array
integer :: i

do i =1,10
test_array(i) = 0.1
write (*,*) test_array(i)
end do
end program test

