module list_m
implicit none

type::lst_t
  integer::n_lst=0
  integer,dimension(:),allocatable::lst
end type

contains



subroutine alloc_lst(lst,n_elm,n_lst)
  type(lst_t),dimension(:),allocatable::lst
  integer,intent(in)::n_elm
  integer,intent(in),optional::n_lst
  integer::l
  call dealloc_lst(lst)
  allocate(lst(n_elm))
  do l=1,n_elm
    lst(l)%n_lst=0
  end do
  if(present(n_lst))then
    do l=1,n_elm
      allocate(lst(l)%lst(n_lst))
      lst(l)%lst=0
    end do
  end if
end subroutine



subroutine dealloc_lst(lst)
  type(lst_t),dimension(:),allocatable::lst
  integer::l
  if(allocated(lst))then
    do l=1,size(lst)
      if(allocated(lst(l)%lst))then
        deallocate(lst(l)%lst)
      end if
    end do
    deallocate(lst)
  end if
end subroutine



function shortest_list(lst1,lst2) result(i)
  type(lst_t),intent(in)::lst1
  type(lst_t),dimension(:),intent(inout)::lst2
  integer::i
  integer::n,i1,l2,nl
  n=huge(n); i=1
  do i1=1,size(lst1%lst)
    l2=lst1%lst(i1)
    nl=size(lst2(l2)%lst)
    if(nl<n)then
      n=nl; i=i1
    end if
  end do
end function



subroutine invert_list(lst1,lst2)
  use alloc_m
  type(lst_t),dimension(:),intent(in)::lst1
  type(lst_t),dimension(:),allocatable::lst2
  integer::l1,l2,i1,i2,n1,n2,n_lst1,n_lst2

  n_lst1=size(lst1)
  n_lst2=0
  do l1=1,n_lst1
    n_lst2=max(n_lst2,maxval(lst1(l1)%lst))
  end do
  call alloc_lst(lst2,n_lst2)

  do l1=1,n_lst1
    n1=size(lst1(l1)%lst)
    do i1=1,n1
      l2=lst1(l1)%lst(i1)
      lst2(l2)%n_lst=lst2(l2)%n_lst+1
    end do
  end do

  do l2=1,n_lst2
    n2=lst2(l2)%n_lst
    allocate(lst2(l2)%lst(n2))
    lst2(l2)%n_lst=0
  end do

  do l1=1,n_lst1
    n1=size(lst1(l1)%lst)
    do i1=1,n1
      l2=lst1(l1)%lst(i1)
      lst2(l2)%n_lst=lst2(l2)%n_lst+1
      lst2(l2)%lst(lst2(l2)%n_lst)=l1
    end do
  end do
end subroutine



subroutine match_lists(lst1,lst2,n_mtch,cl,loc1)
  integer,dimension(:),intent(in)::lst1
  integer,dimension(:),intent(in)::lst2
  integer,intent(out)::n_mtch
  integer,intent(out),optional::cl
  integer,dimension(:),optional::loc1

  integer::l1,l2,n1,n2,n_max
  integer,dimension(size(lst1))::lloc1

  lloc1=0
  n1=size(lst1); n2=size(lst2); n_max=min(n1,n2)
  n_mtch=0

  lp1: do l1=1,n1
    lp2: do l2=1,n2
      if(n_mtch==n_max)then
        exit lp1
      end if
      if(lst1(l1)==lst2(l2))then
        lloc1(l1)=l2
        n_mtch=n_mtch+1
        exit lp2
      end if
    end do lp2
  end do lp1

  if(present(cl))then
    cl=0
    do l1=1,n1
      if(lloc1(l1)==0)then
        cl=l1; exit
      end if
    end do
  end if
  if(present(loc1))loc1=lloc1
end subroutine



subroutine list_to_matrix(arc,ja,ia)
  use sort_m

  type(lst_t),dimension(:),intent(in)::arc
  integer,dimension(:),allocatable::ja,ia

  type(lst_t),dimension(:),allocatable::node
  integer::irow,icol,nrow,ncol,row(7),ic,jc,jl,j,il

  call invert_list(arc,node)

  allocate(ia( size(node)+1 ))
  ic=1; ia(ic)=1

  allocate(ja( size(arc)*2+size(node) ))
  jc=0


  do irow=1,size(node)

    ncol=1; row(ncol)=irow
    do jl=1,size(node(irow)%lst)

      j=node(irow)%lst(jl)
      do il=1,size(arc(j)%lst)

        icol=arc(j)%lst(il)
        if(irow/=icol)then
          ncol=ncol+1
          row(ncol)=icol
        end if
      end do
    end do

    ! sorting is not necessary
!     call isorti(ncol,row(:ncol))
    call hpsorti (ncol, row(:ncol))

    ic=ic+1
    ia(ic)=ia(ic-1)+ncol

    jc=jc+ncol
    ja(jc-ncol+1:jc)=row(:ncol)
  end do
end subroutine

end module
