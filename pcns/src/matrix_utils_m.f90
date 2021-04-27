module matrix_utils_m
use constants_m
implicit none

contains


subroutine matrix_inv(mat,inv)
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
    inv = zero
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


function csr_index(i,j,ja,ia) result(kij)
  integer,intent(in)::i,j
  integer,dimension(:),intent(in)::ja,ia
  integer::kij

  do kij=ia(i),ia(i+1)-1
    select case(ja(kij)-j)
    case(0)
      return
    end select
  end do

  kij=-1
end function



subroutine csr_matvec_prod(a_ij,phi,ja,ia,res)
  integer,dimension(:)::ja,ia
  real(rp),dimension(:)::a_ij
  real(rp),dimension(:)::phi

  real(rp),dimension(:),intent(inout)::res
  integer::i,k,n
  n=size(ia)-1

  do i=1,n
    res(i)=zero
    do k=ia(i),ia(i+1)-1
      res(i)=res(i)+a_ij(k)*phi(ja(k))
    end do
  end do
end subroutine



function csr_diag(a_ij,ja,ia) result(diag)
  integer,dimension(:)::ja,ia
  real(rp),dimension(:)::a_ij

  real(rp),dimension(size(ia)-1)::diag
  integer::i,kii,n
  n=size(ia)-1

!$AD II-LOOP
  do i=1,n
    kii=csr_index(i,i,ja,ia)
    diag(i)=a_ij(kii)
  end do
end function



subroutine csr_transpose(a_ij,ja,ia,at_ij)
  integer,dimension(:)::ja,ia
  real(rp),dimension(:)::a_ij,at_ij
  integer::i,kij,j,kji,n
  n=size(ia)-1

  do i=1,n
    do kij=ia(i),ia(i+1)-1
      j=ja(kij)
      kji=csr_index(j,i,ja,ia)
      at_ij(kji)=a_ij(kij)
    end do
  end do
end subroutine



subroutine csr_to_csc(ja,ia,jao,iao,a,ao)
  integer,dimension(:)::ia,iao,ja,jao
  real(rp),dimension(:),optional::a,ao

  integer i,j,k,n,next,ipos

  n=size(ia)-1
  iao = 0

  do i=1, n
    do k=ia(i), ia(i+1)-1
      j = ja(k)+1
      iao(j) = iao(j)+1
    end do
  end do

  iao(1) = 1
  do i=1,n
    iao(i+1) = iao(i) + iao(i+1)
  end do

  do i=1,n
    do k=ia(i),ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if(present(a).and.present(ao))then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next+1
    end do
  end do

  do i=n,1,-1
    iao(i+1) = iao(i)
  end do
  iao(1) = 1
end subroutine
end module
