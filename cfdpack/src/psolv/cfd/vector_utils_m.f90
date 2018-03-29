module vector_utils_m
use constants_m
implicit none

contains


function p_norm(v1,p) result(s)
  real(rp),dimension(:),intent(in)::v1
  integer,intent(in)::p
  real(rp)::s
  s=sum(abs(v1)**p)**(one/p)
end function



function vec_mag(v1) result(s)
  real(rp),dimension(:),intent(in)::v1
  real(rp)::s
  s=sqrt(dot_product(v1,v1))
end function




function unit_vec(v1) result(v)
  real(rp),dimension(:),intent(in)::v1
  real(rp),dimension(size(v1))::v
  v=v1/(vec_mag(v1)+small)
end function




function cross_prod(v1,v2) result(v)
  real(rp),dimension(3),intent(in)::v1,v2
  real(rp),dimension(3)::v
  integer,dimension(3),parameter::i1=(/2,3,1/)
  integer,dimension(3),parameter::i2=(/3,1,2/)
  v=v1(i1)*v2(i2)-v2(i1)*v1(i2)
end function



function vvec_mag(v1) result(v)
  real(rp),dimension(:,:),intent(in)::v1
  real(rp),dimension(size(v1,2))::v
  integer::i
  do i=1,size(v)
    v(i)=vec_mag(v1(:,i))
  end do
end function



function rev_order(v1) result(v)
  integer,dimension(:),intent(in)::v1
  integer,dimension(size(v1))::v
  integer::i,j
  j=0
  do i=size(v1),1,-1
    j=j+1
    v(i)=v1(j)
  end do
end function




function first_loc(v1,s1) result(s)
  integer,dimension(:),intent(in)::v1
  integer,intent(in)::s1
  integer::s,i
  s=0
  do i=1,size(v1)
    select case(v1(i)-s1)
    case(0)
      s=i
      return
    end select
  end do
end function



! C
subroutine setcp (n_lst, lst, iset)
  integer,intent(in)::n_lst
  integer,dimension(0:n_lst-1)::lst
  integer,intent(in)::iset

  integer::i, j, jj, k, l

  i = 0;
  do

    do j = i + 1, n_lst - 1
      jj = j;
      if (lst(j) /= lst(i)) then
        exit;
      end if
      if (j == n_lst - 1) jj = jj + 1;
    end do

    j = jj;

    do k = i + 1, j - 1
      lst(k) = iset;
    end do

    i = j;
    if (i >= n_lst - 1) exit;
  end do
end subroutine



! C
subroutine hpsortf (n, ra)
  integer,intent(in)::n
!   integer,dimension(n)::ra
  real(rp),dimension(n)::ra

  integer::l,ir,i,j
!   integer::rra
  real(rp)::rra

  l=n/2+1
  ir=n

10 continue
  if(l > 1)then
    l=l-1
    rra=ra(l)
  else
    rra=ra(ir)
    ra(ir)=ra(1)
    ir=ir-1
    if(ir.eq.1)then
      ra(1)=rra
      return
    end if
  end if
  i=l
  j=l+l

20 if(j.le.ir)then
    if(j < ir)then
      if(ra(j) < ra(j+1))  j=j+1
    end if
    if(rra < ra(j))then
      ra(i)=ra(j)
      i=j; j=j+j
    else
      j=ir+1
    end if
    goto 20
  end if

  ra(i)=rra
  goto 10
end subroutine



! C
subroutine hpsorti (n, ra)
  integer,intent(in)::n
  integer,dimension(n)::ra
!   real(rp),dimension(n)::ra

  integer::l,ir,i,j
  integer::rra
!   real(rp)::rra

  l=n/2+1
  ir=n

10 continue
  if(l > 1)then
    l=l-1
    rra=ra(l)
  else
    rra=ra(ir)
    ra(ir)=ra(1)
    ir=ir-1
    if(ir.eq.1)then
      ra(1)=rra
      return
    end if
  end if
  i=l
  j=l+l

20 if(j.le.ir)then
    if(j < ir)then
      if(ra(j) < ra(j+1))  j=j+1
    end if
    if(rra < ra(j))then
      ra(i)=ra(j)
      i=j; j=j+j
    else
      j=ir+1
    end if
    goto 20
  end if

  ra(i)=rra
  goto 10
end subroutine



! C
subroutine uniq (m0, lst0, m1, lst1)
  integer, intent(in) :: m0
  integer, dimension(:), intent(in) :: lst0
  integer, intent(out) :: m1
  integer, dimension(:), allocatable :: lst1

  integer :: i, int_max

  allocate (lst1(m0))
  lst1 = lst0

  int_max = huge (i)
  call hpsorti (m0, lst1)
  call setcp (m0, lst1, int_max)
  call hpsorti (m0, lst1)

  m1 = 0
  do i = 1, m0
    if (lst1(i) < int_max) m1 = m1 + 1
  end do
end subroutine




subroutine iuniq (m0, lst0)
  integer, intent(inout) :: m0
  integer, dimension(:), intent(inout) :: lst0

  integer :: i, int_max, m1

  int_max = huge (i)

  call hpsorti (m0, lst0)
  call setcp (m0, lst0, int_max)
  call hpsorti (m0, lst0)

  m1 = 0
  do i = 1, m0
    if (lst0(i) < int_max) m1 = m1 + 1
  end do
  m0 = m1
end subroutine
end module
