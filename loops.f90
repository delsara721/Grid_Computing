
program test 
implicit none

integer, dimension (20) :: test_array
integer :: i
do i=1,20
test_array(i) = i
end do
call print_array(test_array, 20)
contains 

subroutine print_array(test_array,big)

integer, intent (in):: big
integer, intent (in):: test_array(big)
integer :: i

do i =1, big
write (*,*) test_array(i)
end do 
end subroutine 
end program test 

