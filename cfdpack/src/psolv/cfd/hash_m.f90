module hash_m
implicit none

integer,parameter,private::nullkey = -1

type::hash_t
  integer::key = nullkey
  integer::val = 0
end type
contains


subroutine hash_initialise(hash,n)
  type(hash_t) ,dimension(:) ,allocatable::hash
  integer::n
  integer::m,i

  call hash_finalise(hash)

  m=(n*3)/2
  allocate(hash(m))
end subroutine



subroutine hash_finalise(hash)
  type(hash_t) ,dimension(:) ,allocatable::hash

  if(allocated(hash))deallocate(hash)
end subroutine



function hash_index(key,n) result(indx)
  integer::key
  integer::n
  integer::indx

  indx = mod(key,n)
  indx = max(indx,1)
end function


subroutine hash_uniquekey_insert(hash,key,val)
  type(hash_t) ,dimension(:)::hash
  integer::key
  integer::val
  integer::val_dupl

  call hash_insert(hash,key,val)
  call hash_get_value(hash,key,val_dupl)

  select case(abs(val-val_dupl))
  case(1:)
    call hash_delete(hash,key,val)
  end select
end subroutine



subroutine hash_insert(hash,key,val)
  type(hash_t) ,dimension(:)::hash
  integer::key
  integer::val
  integer::indx,n

  n = size(hash)
  indx = hash_index(key,n)
  do
    select case(hash(indx)%key-nullkey)
    case(0)
      exit
    end select
    indx=indx+1
    indx=hash_index(indx,n)
  end do
  hash(indx)%key = key
  hash(indx)%val = val
end subroutine


subroutine hash_get_value(hash,key,val)
  type(hash_t) ,dimension(:)::hash
  integer::key
  integer::val
  integer::indx,n

  n = size(hash)
  indx = hash_index(key,n)
  do
    select case(hash(indx)%key-key)
    case(0)
      exit
    end select
    indx=indx+1
    indx=hash_index(indx,n)
  end do
  val = hash(indx)%val
end subroutine


subroutine hash_delete_key(hash,key)
  type(hash_t) ,dimension(:)::hash
  integer::key
  integer::indx,n

  n = size(hash)
  indx = hash_index(key,n)
  do
    select case(hash(indx)%key-key)
    case(0)
      exit
    end select
    indx=indx+1
    indx=hash_index(indx,n)
  end do
  hash(indx)%key = nullkey
  hash(indx)%val = 0
end subroutine



subroutine hash_delete(hash,key,val)
  type(hash_t) ,dimension(:)::hash
  integer::key
  integer::val
  integer::indx,n

  n = size(hash)
  indx = hash_index(key,n)
  do
    select case(hash(indx)%key-key)
    case(0)
      select case(hash(indx)%val-val)
      case(0)
        exit
      end select
    end select
    indx=indx+1
    indx=hash_index(indx,n)
  end do
  hash(indx)%key = nullkey
  hash(indx)%val = 0
end subroutine
end module


! program main
!   use hash_m
!   type(hash_t) ,dimension(:) ,allocatable::hash
! !   integer,dimension(7)::elem_indx=(/23,18,41,76,54,12,98/)
!   integer,dimension(7)::elem_indx=(/23,18,23,76,54,23,98/)
!   integer::key,val
!
!   ! initialise table larger than required size
!   call hash_initialise(hash,size(elem_indx))
!
!   ! populate table
!   print*
!   do i=1,size(elem_indx)
!     key = elem_indx(i)
!     val = i
!     call hash_uniquekey_insert(hash,key,val)
!   end do
!
!   ! print table
!   print*
!   do i=1,size(hash)
!     print*,i,hash(i)%key,hash(i)%val
!   end do
!
!   ! find value
!   key = 54
!   call hash_get_value(hash,key,val)
!   print*
!   print*,key,val
!
!   ! delete entry
!   call hash_delete_key(hash,key)
!
!   ! print table
!   print*
!   do i=1,size(hash)
!     print*,i,hash(i)%key,hash(i)%val
!   end do
!
!   call hash_finalise(hash)
! end program
