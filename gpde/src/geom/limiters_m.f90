module limiters_m
implicit none

! umin umax lim
! real(rp),dimension(:,:),allocatable::wk1,wk2,wk3

contains



!==============================================================================!
function barth(d2,d1min,d1max) result(res)
!------------------------------------------------------------------------------!
  use const_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real(rp):: d2, d1min, d1max
!------------------------------------------------------------------------------!
  real(rp):: num, den, res
!------------------------------------------------------------------------------!
  if (d2 > small) then
    num    = d1max
    den    = d2
    res = min(1.d0,num/(den+small))
  else if (d2 < -small) then
    num    = d1min
    den    = d2
    res = min(1.d0,num/(den+small))
  else
    res = 1.d0
  endif
end function



!==============================================================================!
function venkat(d2,d1min,d1max,eps2) result(res)
!------------------------------------------------------------------------------!
  use const_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  real(rp):: d2, d1min, d1max, eps2
!------------------------------------------------------------------------------!
  real(rp):: num, den, res
!------------------------------------------------------------------------------!
  if (d2 > small) then
    num    = (d1max*d1max+eps2)*d2 + 2.d0*d2*d2*d1max
    den    = d2*(d1max*d1max+2.d0*d2*d2+d1max*d2+eps2)
    res = num/(den+small)
  else if (d2 < -small) then
    num    = (d1min*d1min+eps2)*d2 + 2.d0*d2*d2*d1min
    den    = d2*(d1min*d1min+2.d0*d2*d2+d1min*d2+eps2)
    res = num/(den+small)
  else
    res = 1.d0
  endif
end function



!==============================================================================!
subroutine local_extrema(fce_elm,u,umin,umax)
!------------------------------------------------------------------------------!
  use const_m
  use list_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(lst_t),dimension(:),intent(in)::fce_elm
  real(rp),dimension(:,:),intent(in)::u
  real(rp),dimension(:,:),intent(inout)::umin
  real(rp),dimension(:,:),intent(inout)::umax
!------------------------------------------------------------------------------!
  integer::i,j,k,ib,ic,n,nc,nk
!------------------------------------------------------------------------------!
  n = size(umin,2)
  nc = size(u,1)
  nk = size(fce_elm)

  do i=1,n
    umin(:,i) = u(:,i)
    umax(:,i) = u(:,i)
  end do

!#AD NO-II-LOOP
  do k=1,nk
    i=fce_elm(k)%lst(1)
    j=fce_elm(k)%lst(2)

    if(j>0)then
      do ic=1,nc
        umin(ic,i)=min(umin(ic,i),u(ic,j))
        umin(ic,j)=min(umin(ic,j),u(ic,i))

        umax(ic,i)=max(umax(ic,i),u(ic,j))
        umax(ic,j)=max(umax(ic,j),u(ic,i))
      end do
    else
      ib=abs(j)
      do ic=1,nc
        umin(ic,i)=min(umin(ic,i),u(ic,ib))
        umax(ic,i)=max(umax(ic,i),u(ic,ib))
      end do
    end if
  end do
end subroutine



!==============================================================================!
subroutine grad_limiter(fce_elm,geom,u,gradu,umin,umax,lim,fac,ctrl)
!------------------------------------------------------------------------------!
  use const_m
  use list_m
  use geom_data_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(lst_t),dimension(:),intent(in)::fce_elm
  type(geom_t),intent(in)::geom
  real(rp),dimension(:,:)::u
  real(rp),dimension(:,:,:)::gradu
  real(rp),dimension(:,:)::umin
  real(rp),dimension(:,:)::umax
  real(rp),intent(in)::fac
  real(rp),dimension(:,:)::lim
  integer,intent(in)::ctrl
!------------------------------------------------------------------------------!
  real(rp)::fac3
  real(rp),dimension(2)::eps2
  real(rp),dimension(size(gradu,1))::dx
  real(rp),dimension(size(u,1),2)::d1min,d1max,d2,limval
  integer::i,j,k,nk,ic,nc,ndx
!------------------------------------------------------------------------------!
  lim = 1.d0
  fac3 = fac**3
  ndx = size(gradu,1)
  nc = size(u,1)
  nk = size(fce_elm)

!#AD NO-II-LOOP
  do k=1,nk
    i=fce_elm(k)%lst(1)
    j=fce_elm(k)%lst(2)

    if(j>0)then
      eps2(1) = fac3*geom%vol(i)**(3.d0/ndx)
      eps2(2) = fac3*geom%vol(j)**(3.d0/ndx)

      do ic=1,nc
        d1min(ic,1) = umin(ic,i) - u(ic,i)
        d1max(ic,1) = umax(ic,i) - u(ic,i)

        d1min(ic,2) = umin(ic,j) - u(ic,j)
        d1max(ic,2) = umax(ic,j) - u(ic,j)
      end do

      dx = geom%x_fc(:,k)-geom%x_vc(:,i)
      do ic=1,nc
        d2(ic,1) = sum(gradu(:,ic,i)*dx)
        d2(ic,2) = sum(gradu(:,ic,j)*(-dx))
      end do

      select case(ctrl)
      case (1)
        do ic=1,nc
          limval(ic,1) = barth(d2(ic,1),d1min(ic,1),d1max(ic,1))
          limval(ic,2) = barth(d2(ic,2),d1min(ic,2),d1max(ic,2))
        end do
      case (2)
        do ic=1,nc
          limval(ic,1) = venkat(d2(ic,1),d1min(ic,1),d1max(ic,1),eps2(1))
          limval(ic,2) = venkat(d2(ic,2),d1min(ic,2),d1max(ic,2),eps2(2))
        end do
      end select

      do ic=1,nc
        lim(ic,i) = min(limval(ic,1),lim(ic,i))
        lim(ic,j) = min(limval(ic,2),lim(ic,j))
      end do
    end if
  end do
end subroutine

end module
