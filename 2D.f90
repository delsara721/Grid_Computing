program main
include 'mpif.h'


real, dimension(0:3) :: X
real, dimension (0:3) :: Y
real, dimension (0:3,0:3) :: ul_new
real, dimension (0:3,0:3) :: ul_old
real, dimension (0:3,0:3) :: dl_new
real, dimension (0:3,0:3) :: dl_old
real, dimension (0:3,0:3) :: ur_new
real, dimension (0:3,0:3) :: ur_old
real, dimension (0:3,0:3) :: dr_new
real, dimension (0:3,0:3) :: dr_old
real:: r,s
integer :: world_rank,i,j,k,period, ierror, status(MPI_STATUS_SIZE)
integer :: ROWS_X = 3, ROWS_Y = 3
real:: c = 1.0, dx=1.0, dy=1.0, dt=0.1
real:: MAX_TIME = 0.3
r=(c*dt)/(dx*dx)
s=(c*dt)/(dy*dy)
period = MAX_TIME/dt
world_rank = 0
X=(/-1.5,-0.5,0.5,1.5/)
Y=(/-1.5,-0.5,0.5,1.5/)

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,world_size,ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,world_rank, ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (world_rank ==0) then 
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
ul_old(i,j) = exp(-1*((X(i-1)*X(i-1))+(Y(j-1)*Y(j-1))))

call MPI_Recv(ul_old(0,j), 1, MPI_REAL, 2, 0, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(ul_old(ROWS_X,j), 1, MPI_REAL, 2, 1, MPI_COMM_WORLD, status,ierror)
call MPI_Send(ul_old(ROWS_X-1,j), 1, MPI_REAL, 2, 4, MPI_COMM_WORLD,ierror)
call MPI_Send(ul_old(1,j), 1, MPI_REAL, 2, 5, MPI_COMM_WORLD, ierror)
	enddo	
call MPI_Recv(ul_old(i,0), 1, MPI_REAL, 1, 8, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(ul_old(i,ROWS_Y), 1, MPI_REAL, 1, 9, MPI_COMM_WORLD, status, ierror)
call MPI_Send(ul_old(i,ROWS_Y-1), 1, MPI_REAL, 1, 12, MPI_COMM_WORLD, ierror)
call MPI_Send(ul_old(i,1), 1, MPI_REAL, 1, 13, MPI_COMM_WORLD,ierror) 
enddo
call printT(ul_old,ROWS_X,ROWS_Y)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (world_rank ==1) then
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
ur_old(i,j) = exp(-1*((X(i-1)*X(i-1))+(Y(j+1)*Y(j+1))))

call MPI_Recv(ur_old(0,j), 1, MPI_REAL, 3, 2, MPI_COMM_WORLD,status, ierror)
call MPI_Recv(ur_old(ROWS_X,j), 1, MPI_REAL, 3, 3, MPI_COMM_WORLD, status, ierror)
call MPI_Send(ur_old(ROWS_X-1,j), 1, MPI_REAL, 3, 6, MPI_COMM_WORLD, ierror)
call MPI_Send(ur_old(1,j), 1, MPI_REAL, 3, 7, MPI_COMM_WORLD , ierror)
enddo

call MPI_Send(ur_old(i,ROWS_Y-1), 1, MPI_REAL, 0, 8, MPI_COMM_WORLD, ierror)
call MPI_Send(ur_old(i,1), 1, MPI_REAL, 0, 9, MPI_COMM_WORLD, ierror)
call MPI_Recv(ur_old(i,0), 1, MPI_REAL, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, status, ierror)
call MPI_Recv(ur_old(i,ROWS_Y), 1, MPI_REAL, 0, 13, MPI_COMM_WORLD,MPI_STATUS_IGNORE ,status,ierror)
enddo
call printT(ur_old,ROWS_X,ROWS_Y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else if (world_rank ==2) then
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
dl_old(i,j) = exp(-1*((X(i+1)*X(i+1))+(Y(j-1)*Y(j-1))));

call MPI_Send(dl_old(ROWS_X-1,j), 1, MPI_REAL, 0,0, MPI_COMM_WORLD, ierror)
call MPI_Send(dl_old(1,j), 1, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierror)
call MPI_Recv(dl_old(0,j), 1, MPI_REAL, 0, 4, MPI_COMM_WORLD, status , ierror)
call MPI_Recv(dl_old(ROWS_X,j), 1, MPI_REAL, 0, 5, MPI_COMM_WORLD, status, ierror)
enddo

call MPI_Recv(dl_old(i,0), 1, MPI_REAL, 3, 10, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(dl_old(i,ROWS_Y), 1, MPI_REAL, 3, 11, MPI_COMM_WORLD, status, ierror)
call MPI_Send(dl_old(i,ROWS_Y-1), 1, MPI_REAL, 3, 14, MPI_COMM_WORLD, ierror)
call MPI_Send(dl_old(i,1), 1, MPI_REAL, 3, 15, MPI_COMM_WORLD, ierror)
enddo
call printT(dl_old,ROWS_X,ROWS_Y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else if (world_rank==3) then
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
dr_old(i,j) = exp(-1*((X(i+1)*X(i+1))+(Y(j+1)*Y(j+1))))

call MPI_Send(dr_old(ROWS_X-1,j), 1, MPI_REAL, 1, 2, MPI_COMM_WORLD, ierror )
call MPI_Send(dr_old(1,j), 1, MPI_REAL, 1, 3, MPI_COMM_WORLD, ierror)
call MPI_Recv(dr_old(0,j), 1, MPI_REAL, 1, 6, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(dr_old(ROWS_X,j), 1, MPI_REAL, 1, 7, MPI_COMM_WORLD, status, ierror)
enddo
 
call MPI_Send(dr_old(i,ROWS_Y-1), 1, MPI_REAL, 2, 10, MPI_COMM_WORLD, ierror)
call MPI_Send(dr_old(i,1), 1, MPI_REAL, 2, 11, MPI_COMM_WORLD, ierror)
call MPI_Recv(dr_old(i,0), 1, MPI_REAL, 2, 14, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(dr_old(i,ROWS_Y), 1, MPI_REAL, 2, 15, MPI_COMM_WORLD, status, ierror)
enddo
call printT(dr_old,ROWS_X,ROWS_Y)
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=1,period
if (world_rank ==0) then
call fillul()
call printT(ul_new,ROWS_X,ROWS_Y)
else if (world_rank ==1) then 
call fillur()
call printT(ur_new,ROWS_X,ROWS_Y)
else if (world_rank ==2) then 
call filldl()
call printT(dl_new,ROWS_X,ROWS_Y)
else if (world_rank ==3) then
call filldr()
call printT(dr_new,ROWS_X,ROWS_Y)
end if
enddo

call MPI_FINALIZE(ierror)
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine printT(array,xd,yd)
integer :: i,j
integer:: xd, yd
real, dimension (0:xd,0:yd) :: array

do i=1, xd-1
	do j=1, yd-1
		write(*,*) array(i,j)
	enddo
enddo
write (*,*) 
end subroutine printT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fillul()
integer :: i,j,f,g
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
ul_new(i,j) = r*ul_old(i-1,j)+r*ul_old(i+1,j)+(1-2*r-2*s) *ul_old(i,j)+s*ul_old(i,j-1)+s*ul_old(i,j+1) 

call MPI_Recv(ul_new(0,j), 1, MPI_REAL, 2, 0, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(ul_new(ROWS_X,j), 1, MPI_REAL, 2, 1, MPI_COMM_WORLD, status,ierror)
call MPI_Send(ul_new(ROWS_X-1,j), 1, MPI_REAL, 2, 4, MPI_COMM_WORLD,ierror)
call MPI_Send(ul_new(1,j), 1, MPI_REAL, 2, 5, MPI_COMM_WORLD, ierror)
	enddo	
call MPI_Recv(ul_new(i,0), 1, MPI_REAL, 1, 8, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(ul_new(i,ROWS_Y), 1, MPI_REAL, 1, 9, MPI_COMM_WORLD, status, ierror)
call MPI_Send(ul_new(i,ROWS_Y-1), 1, MPI_REAL, 1, 12, MPI_COMM_WORLD, ierror)
call MPI_Send(ul_new(i,1), 1, MPI_REAL, 1, 13, MPI_COMM_WORLD,ierror) 
enddo

!copy new to old
do f =0, 3
do g =0, 3
	ul_old(f,g) = ul_new(f,g) 
enddo
enddo
end subroutine fillul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fillur()
integer :: i,j,f,g
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
ur_new(i,j) = r*ur_old(i-1,j)+r*ur_old(i+1,j)+(1-2*r-2*s) *ur_old(i,j)+s*ur_old(i,j-1)+s*ur_old(i,j+1) 


call MPI_Recv(ur_new(0,j), 1, MPI_REAL, 3, 2, MPI_COMM_WORLD,status, ierror)
call MPI_Recv(ur_new(ROWS_X,j), 1, MPI_REAL, 3, 3, MPI_COMM_WORLD, status, ierror)
call MPI_Send(ur_new(ROWS_X-1,j), 1, MPI_REAL, 3, 6, MPI_COMM_WORLD, ierror)
call MPI_Send(ur_new(1,j), 1, MPI_REAL, 3, 7, MPI_COMM_WORLD , ierror)
enddo

call MPI_Send(ur_new(i,ROWS_Y-1), 1, MPI_REAL, 0, 8, MPI_COMM_WORLD, ierror)
call MPI_Send(ur_new(i,1), 1, MPI_REAL, 0, 9, MPI_COMM_WORLD, ierror)
call MPI_Recv(ur_new(i,0), 1, MPI_REAL, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE, status, ierror)
call MPI_Recv(ur_new(i,ROWS_Y), 1, MPI_REAL, 0, 13, MPI_COMM_WORLD,MPI_STATUS_IGNORE ,status,ierror)
enddo

!copy new to old
do f =0, 3
do g =0, 3 
	ur_old(f,g) = ur_new(f,g) 
enddo
enddo
end subroutine fillur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filldl()
integer :: i,j,f,g
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
dl_new(i,j) = r*dl_old(i-1,j)+r*dl_old(i+1,j)+(1-2*r-2*s) *dl_old(i,j)+s*dl_old(i,j-1)+s*dl_old(i,j+1) 


call MPI_Send(dl_new(ROWS_X-1,j), 1, MPI_REAL, 0,0, MPI_COMM_WORLD, ierror)
call MPI_Send(dl_new(1,j), 1, MPI_REAL, 0, 1, MPI_COMM_WORLD, ierror)
call MPI_Recv(dl_new(0,j), 1, MPI_REAL, 0, 4, MPI_COMM_WORLD, status , ierror)
call MPI_Recv(dl_new(ROWS_X,j), 1, MPI_REAL, 0, 5, MPI_COMM_WORLD, status, ierror)
enddo

call MPI_Recv(dl_new(i,0), 1, MPI_REAL, 3, 10, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(dl_new(i,ROWS_Y), 1, MPI_REAL, 3, 11, MPI_COMM_WORLD, status, ierror)
call MPI_Send(dl_new(i,ROWS_Y-1), 1, MPI_REAL, 3, 14, MPI_COMM_WORLD, ierror)
call MPI_Send(dl_new(i,1), 1, MPI_REAL, 3, 15, MPI_COMM_WORLD, ierror)
enddo

!copy new to old
do f =0, 3
do g =0, 3

	dl_old(f,g) = dl_new(f,g) 
 
enddo
enddo

end subroutine filldl 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine filldr()
integer :: i,j,f,g
do i=1, ROWS_X-1
	do j =1, ROWS_Y-1
dr_new(i,j) = r*dr_old(i-1,j)+r*dr_old(i+1,j)+(1-2*r-2*s) *dr_old(i,j)+s*dr_old(i,j-1)+s*dr_old(i,j+1)

call MPI_Send(dr_new(ROWS_X-1,j), 1, MPI_REAL, 1, 2, MPI_COMM_WORLD, ierror )
call MPI_Send(dr_new(1,j), 1, MPI_REAL, 1, 3, MPI_COMM_WORLD, ierror)
call MPI_Recv(dr_new(0,j), 1, MPI_REAL, 1, 6, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(dr_new(ROWS_X,j), 1, MPI_REAL, 1, 7, MPI_COMM_WORLD, status, ierror)
enddo
 
call MPI_Send(dr_new(i,ROWS_Y-1), 1, MPI_REAL, 2, 10, MPI_COMM_WORLD, ierror)
call MPI_Send(dr_new(i,1), 1, MPI_REAL, 2, 11, MPI_COMM_WORLD, ierror)
call MPI_Recv(dr_new(i,0), 1, MPI_REAL, 2, 14, MPI_COMM_WORLD, status, ierror)
call MPI_Recv(dr_new(i,ROWS_Y), 1, MPI_REAL, 2, 15, MPI_COMM_WORLD, status, ierror)
enddo
!copy new to old
do f =0, 3
do g =0, 3

	dr_old(f,g) = dr_new(f,g) 
enddo
enddo
end subroutine filldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program main

