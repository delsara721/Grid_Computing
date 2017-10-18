program main
implicit none 
real, dimension (4) :: position_array
real, dimension (4) :: time_a
real, dimension (4, 6):: T
integer :: ROWS =4, COLS = 6, ROWS_X =4
real :: MIN_POSITION = -1.5, MAX_POSITION = 1.5
real :: MAX_TIME = 0.3, MIN_TIME = 0.0, dx = 1, dt = 0.1
real :: constant=1.0

call fillTime()
call fillPosition()
call base()
call fill()
call print_array()

contains 


subroutine fillPosition()
integer :: counter 
real :: i
i = MIN_POSITION
do counter =1, ROWS_X
position_array(counter) = i
i = i+dx
end do
!call printPosition()
end subroutine fillPosition
subroutine printPosition()
integer :: j
do j = 1, ROWS_X
write (*,*) position_array(j)
end do
end subroutine printPosition

subroutine fillTime()
integer :: counter
real :: i 
i = MIN_TIME
do counter=1, ROWS
time_a(counter) = i
i = i +dt
end do
!call printTime()
end subroutine fillTime
subroutine printTime()
integer :: j
do j = 1, ROWS
write (*,*) time_a(j)
end do
end subroutine printTime


subroutine base()
real :: result_1
integer :: i
do i=2,ROWS_X+1
result_1 = cos(position_array(i-1))
T(1,i) = result_1
end do
T(1,1)=T(1,ROWS_X+1)
T(1,COLS)=T(1,2)
end subroutine base

subroutine fill() 
real :: r,k
integer :: i
integer :: j
r = (constant*dt)/(dx*dx)
k = 1-(2*r)

do i=2,ROWS
do j=2, ROWS_X+1
T(i,j) = (r*T(i-1,j-1))+(k*T(i-1,j))+(r*T(i-1,j+1))
end do
T(i,1) = T(i,COLS-1)
T(i,COLS) = T(i,2)
end do
end subroutine fill 



subroutine print_array()
integer :: i,j
do i =1, ROWS
write (*,*) (T(i,j),j=1,COLS)
end do
end subroutine print_array

end program main

