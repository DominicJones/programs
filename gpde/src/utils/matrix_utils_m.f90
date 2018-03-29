module matrix_utils_m
implicit none
contains



subroutine matrix_inv(mat,inv)
  use const_m
  real(rp),dimension(:,:),intent(in)::mat
  real(rp),dimension(:,:),intent(out)::inv
  real(rp)::det
  integer::m,n

  m=size(mat,1)
  n=size(mat,2)

  if(.not.(m==n.and.(m==2.or.m==3)))then
    return
  end if

  select case(m)
  case(2)
    det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
  case(3)
    det=(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))*mat(1,1)
    det=det+(mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3))*mat(1,2)
    det=det+(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))*mat(1,3)
  end select

  if(abs(det) < small)then
    inv = 0.d0
    return
  end if

  select case(m)
  case(2)
    inv(1,1)= mat(2,2)
    inv(1,2)=-mat(1,2)
    inv(2,1)=-mat(2,1)
    inv(2,2)= mat(1,1)
  case(3)
    inv(1,1)=mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
    inv(2,1)=mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
    inv(3,1)=mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
    inv(1,2)=mat(3,2)*mat(1,3)-mat(3,3)*mat(1,2)
    inv(2,2)=mat(3,3)*mat(1,1)-mat(3,1)*mat(1,3)
    inv(3,2)=mat(3,1)*mat(1,2)-mat(3,2)*mat(1,1)
    inv(1,3)=mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
    inv(2,3)=mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
    inv(3,3)=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
  end select

  inv=inv/det
end subroutine

end module
