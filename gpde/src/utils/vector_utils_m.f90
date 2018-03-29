module vector_utils_m
use const_m
implicit none
contains


function p_norm(v1,p) result(s)
  real(rp),dimension(:),intent(in)::v1
  integer,intent(in)::p
  real(rp)::s1,s
  integer::i,n
  n=size(v1)
  s=0.d0
  do i=1,n
    s1=abs(v1(i))
    s=s+s1**p
  end do
  s1=1.d0/p
  s=s**s1
!   s=sum(abs(v1)**p)**(1.d0/p)
end function


function dot_prod(v1,v2) result(s)
  real(rp),dimension(:),intent(in)::v1
  real(rp),dimension(:),intent(in)::v2
  real(rp)::s
  s=sum(v1*v2)
end function


function vec_mag(v1) result(s)
  real(rp),dimension(:),intent(in)::v1
  real(rp)::s
  s=sqrt(dot_prod(v1,v1))
end function


function unit_vec(v1) result(v)
  real(rp),dimension(:),intent(in)::v1
  real(rp),dimension(size(v1))::v
  real(rp)::s
  s=vec_mag(v1)
  s=s+1.d-20
  v=v1/s
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

end module
