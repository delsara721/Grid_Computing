
program test 
implicit none

integer, dimension (2,20) :: test_array
integer :: i,j, add
add = 0
do i=1,20
add = add+i
do j = 1,2
if (j==1) then 
test_array(i,j) = i

else 
test_array(i,j) = add
end if 
end do
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

