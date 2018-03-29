module sort_m
implicit none

contains

subroutine hpsorti (n, ra)
  integer,intent(in)::n
  integer,dimension(n)::ra
!   real(rp),dimension(n)::ra

  integer::l,ir,i,j
  integer::rra
!   real(rp)::rra

  if (n <= 1) return

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
    if (lst0(i) < int_max) then
      m1 = m1 + 1
    end if
  end do
  m0 = m1
end subroutine


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

end module
