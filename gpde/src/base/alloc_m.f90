module alloc_m
use const_m
implicit none

interface alloc
  module procedure alloc_1dr
  module procedure alloc_2dr
  module procedure alloc_3dr
  module procedure alloc_1di
  module procedure alloc_2di
  module procedure alloc_3di
end interface

contains


subroutine alloc_1dr(arr,ref,i)
  real(rp),dimension(:),allocatable::arr
  real(rp),dimension(:),intent(in),optional::ref
  integer,intent(in),optional::i
  if(allocated(arr))deallocate(arr)
  if(present(ref))then
    allocate(arr(size(ref,1)))
    arr=0
  else if(present(i))then
    allocate(arr(i))
    arr=0
  end if
end subroutine


subroutine alloc_2dr(arr,ref,i,j)
  real(rp),dimension(:,:),allocatable::arr
  real(rp),dimension(:,:),intent(in),optional::ref
  integer,intent(in),optional::i,j
  if(allocated(arr))deallocate(arr)
  if(present(ref))then
    allocate(arr(size(ref,1),size(ref,2)))
    arr=0
  else if(present(i).and.present(j))then
    allocate(arr(i,j))
    arr=0
  end if
end subroutine


subroutine alloc_3dr(arr,ref,i,j,k)
  real(rp),dimension(:,:,:),allocatable::arr
  real(rp),dimension(:,:,:),intent(in),optional::ref
  integer,intent(in),optional::i,j,k
  if(allocated(arr))deallocate(arr)
  if(present(ref))then
    allocate(arr(size(ref,1),size(ref,2),size(ref,3)))
    arr=0
  else if(present(i).and.present(j).and.present(k))then
    allocate(arr(i,j,k))
    arr=0
  end if
end subroutine


subroutine alloc_1di(arr,ref,i)
  integer,dimension(:),allocatable::arr
  integer,dimension(:),intent(in),optional::ref
  integer,intent(in),optional::i
  if(allocated(arr))deallocate(arr)
  if(present(ref))then
    allocate(arr(size(ref,1)))
    arr=0
  else if(present(i))then
    allocate(arr(i))
    arr=0
  end if
end subroutine


subroutine alloc_2di(arr,ref,i,j)
  integer,dimension(:,:),allocatable::arr
  integer,dimension(:,:),intent(in),optional::ref
  integer,intent(in),optional::i,j
  if(allocated(arr))deallocate(arr)
  if(present(ref))then
    allocate(arr(size(ref,1),size(ref,2)))
    arr=0
  else if(present(i).and.present(j))then
    allocate(arr(i,j))
    arr=0
  end if
end subroutine


subroutine alloc_3di(arr,ref,i,j,k)
  integer,dimension(:,:,:),allocatable::arr
  integer,dimension(:,:,:),intent(in),optional::ref
  integer,intent(in),optional::i,j,k
  if(allocated(arr))deallocate(arr)
  if(present(ref))then
    allocate(arr(size(ref,1),size(ref,2),size(ref,3)))
    arr=0
  else if(present(i).and.present(j).and.present(k))then
    allocate(arr(i,j,k))
    arr=0
  end if
end subroutine

end module
