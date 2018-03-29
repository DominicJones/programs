module hash_m
implicit none

integer,parameter,private::nullkey = -1

type::hash_data_t
  integer::key
  integer::val = 0
end type
contains



subroutine initialise_hash(table,n)
  type(hash_data_t),dimension(:),allocatable::table
  integer::n,i

  call finalise_hash(table)
  allocate(table(n))
  do i=1,size(table)
    table(i)%key = nullkey
  end do
end subroutine



subroutine finalise_hash(table)
  type(hash_data_t),dimension(:),allocatable::table

  if(allocated(table))deallocate(table)
end subroutine



function hash_index(key,n)
  integer::key,n
  integer::hash_index

  hash_index = mod(key,n)
end function



subroutine insert_hash_value(key,val,table)
  integer::key
  integer::val
  type(hash_data_t),dimension(:)::table
  integer::indx

  indx = hash_index(key,size(table))
  do
    if(table(indx)%key == nullkey)exit
    indx=indx+1
    indx=mod(indx,size(table))
  end do
  table(indx)%key = key
  table(indx)%val = val
end subroutine



subroutine get_hash_value(key,val,table)
  integer::key
  integer::val
  type(hash_data_t),dimension(:)::table
  integer::indx

  indx = hash_index(key,size(table))
  do
    if(table(indx)%key == key)exit
    indx=indx+1
    indx=mod(indx,size(table))
  end do
  val = table(indx)%val
end subroutine



subroutine delete_hash_value(key,table)
  integer::key
  type(hash_data_t),dimension(:)::table
  integer::indx

  indx = hash_index(key,size(table))
  do
    if(table(indx)%key == key)exit
    indx=indx+1
    indx=mod(indx,size(table))
  end do
  table(indx)%key = nullkey
  table(indx)%val = 0
end subroutine
end module


! program main
!   use hash_m
!   type(hash_data_t),dimension(:),allocatable::table
!   integer,dimension(7)::elem_indx=(/23,18,41,76,54,12,98/)
!   integer::key,val
! 
!   ! initialise table larger than required size
!   call initialise_hash(table,size(elem_indx)*15/10)
! 
!   ! populate table
!   do i=1,7
!     key = elem_indx(i)
!     val = i
!     call insert_hash_value(key,val,table)
!   end do
! 
!   ! print table
!   do i=1,size(table)
!     print*,i,table(i)%key,table(i)%val
!   end do
! 
!   ! find value
!   key = 54
!   call get_hash_value(key,val,table)
!   print*,key,val
! 
!   ! delete entry
!   call delete_hash_value(key,table)
! 
!   ! print table
!   do i=1,size(table)
!     print*,i,table(i)%key,table(i)%val
!   end do
! 
!   call finalise_hash(table)
! end program

