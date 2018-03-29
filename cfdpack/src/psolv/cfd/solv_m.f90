module solv_m
use constants_m!, rp => rsp
implicit none

logical,private,save::pri_lin_eqn=.true.

contains


subroutine lin_eqns_solver(a_ij,phi,rhs,norm,ja,ia,n_iter,reduc,isol,id)
  use matrix_utils_m
  use vector_utils_m

  real(rp),dimension(:), intent(in)::a_ij
  real(rp),dimension(:), intent(inout)::phi
  real(rp),dimension(:), intent(inout)::rhs
  real(rp), intent(inout)::norm
  real(rp), intent(in)::reduc
  integer,dimension(:), intent(in)::ja,ia
  integer, intent(in)::n_iter
  integer, intent(in)::isol,id

  select case(isol)
  case(1)
    call cg_solver(a_ij,phi,rhs,norm,ja,ia,n_iter,reduc)
  case (2)
    call bicg_solver(a_ij,phi,rhs,norm,ja,ia,n_iter,reduc)
  end select
end subroutine


subroutine cg_solver(a_ij,phi,res,norm,ja,ia,n_iter,reduc)
  use constants_m
  use matrix_utils_m
  use vector_utils_m

  real(rp),dimension(:), intent(in)::a_ij
  real(rp),dimension(:), intent(inout)::phi
  real(rp),dimension(:), intent(inout)::res
  real(rp), intent(inout)::norm
  real(rp), intent(in)::reduc
  integer,dimension(:), intent(in)::ja,ia
  integer, intent(in)::n_iter

  real(rp),dimension(size(ia)-1,6)::wk
  real(rp),dimension(6)::cof
  real(rp)::norm0
  integer,parameter::p=1
  integer::i,kii,j,kij,kji,iter

  ! calculate residual vector and norm
  call csr_matvec_prod(a_ij,phi,ja,ia,wk(:,1))
  res=res-wk(:,1)
  norm=p_norm(res,p)
  norm0=norm

  ! calculate preconditioning vector: C = M^{-1}
  call preconditioner(a_ij,ja,ia,wk(:,6))

  ! cof: res_n=1 alpha=2 beta=3 s=4
  cof=0; cof(4)=large

  ! wk: res=1 p=2 z=3
  wk(:,1)=res
  wk(:,2:3)=0


  ! iterate:
  do iter=1,n_iter
    cof(1)=0

    ! solve: M z = res
    call direct_solve(a_ij,wk(:,6),wk(:,1),ja,ia,wk(:,3))

    ! beta = 1 / s
    cof(3)=one/(cof(4)+small)

    ! s = p.dot.z
    cof(4)=sum(wk(:,1)*wk(:,3))

    ! beta = beta * s
    cof(3)=cof(3)*cof(4)

    ! p = z + beta.p
    wk(:,2)=wk(:,3)+cof(3)*wk(:,2)

    ! calculate: z = A.p
    call csr_matvec_prod(a_ij,wk(:,2),ja,ia,wk(:,3))

    ! alpha = s / p.dot.z
    cof(2)=sum(wk(:,2)*wk(:,3))
    cof(2)=cof(4)/(cof(2)+small)

    ! phi = phi + alpha.p
    ! res = res - alpha.z
    ! res_n = sum |res|
    phi=phi+cof(2)*wk(:,2)
    wk(:,1)=wk(:,1)-cof(2)*wk(:,3)
    norm=p_norm(wk(:,1),p)

    ! res_n = res_n / res0_n
    cof(1)=norm/(norm0+small)
    if(cof(1)<reduc)exit
  end do
end subroutine




subroutine bicg_solver(a_ij,phi,res,norm,ja,ia,n_iter,reduc)
  use constants_m
  use matrix_utils_m
  use vector_utils_m

  real(rp),dimension(:), intent(in)::a_ij
  real(rp),dimension(:), intent(inout)::phi
  real(rp),dimension(:), intent(inout)::res
  real(rp), intent(inout)::norm
  real(rp), intent(in)::reduc
  integer,dimension(:), intent(in)::ja,ia
  integer, intent(in)::n_iter

  real(rp),dimension(size(ia)-1,6)::wk
  real(rp),dimension(6)::cof
  real(rp)::norm0
  integer,parameter::p=1
  integer::i,kii,j,kij,kji,iter,n

  n=size(ia)-1

  ! calculate residual vector and norm
  call csr_matvec_prod(a_ij,phi,ja,ia,wk(:,1))
  res=res-wk(:,1)
  norm=p_norm(res,p)
  norm0=norm

  ! calculate preconditioning vector: C = M^{-1}
  call preconditioner(a_ij,ja,ia,wk(:,6))

  ! initialise workspace scalars:
  ! res_n=1 alpha=2 beta0=3 beta=4 gam=5 omeg=6
  cof=0; cof(2:3)=1; cof(5)=1

  ! initialise workspace vectors:
  ! res=1 p=2 u=3 v=4 z=5
  wk(:,1)=res; wk(:,2:5)=0

  ! iterate:
  do iter=1,n_iter
    cof(1)=0

    ! beta = res.dot.res0
    ! omega = beta.gamma / alpha.beta0
    ! beta0 = beta
    cof(4)=sum(wk(:,1)*res)
    cof(6)=cof(4)*cof(5)/(cof(2)*cof(3)+small)
    cof(3)=cof(4)

    ! p = res + omega.(p - alpha.u)
    do i=1,n
      wk(i,2)=wk(i,1)+cof(6)*(wk(i,2)-cof(2)*wk(i,3))
    end do

    ! solve: M z = p
    call direct_solve(a_ij,wk(:,6),wk(:,2),ja,ia,wk(:,5))

    ! calculate: u = A.z
    call csr_matvec_prod(a_ij,wk(:,5),ja,ia,wk(:,3))

    ! gamma = beta / u.dot.res0
    cof(5)=sum(wk(:,3)*res)
    cof(5)=cof(4)/(cof(5)+small)

    ! phi = phi + gamma.z
    ! res = res - gamma.u
    phi=phi+cof(5)*wk(:,5)
    wk(:,1)=wk(:,1)-cof(5)*wk(:,3)

    ! solve: M z = res
    call direct_solve(a_ij,wk(:,6),wk(:,1),ja,ia,wk(:,5))

    ! calculate:  v = A.z
    call csr_matvec_prod(a_ij,wk(:,5),ja,ia,wk(:,4))

    ! alpha = v.dot.res / v.dot.v
    cof(1:2)=zero
    do i=1,n
      cof(1)=cof(1)+wk(i,4)*wk(i,1)
      cof(2)=cof(2)+(wk(i,4)+small)**2
    end do
    cof(2)=cof(1)/(cof(2)+small)

    ! phi = phi + alpha.z
    ! res = res - alpha.v
    ! res_n = sum |res|
    phi=phi+cof(2)*wk(:,5)
    wk(:,1)=wk(:,1)-cof(2)*wk(:,4)
    norm=p_norm(wk(:,1),p)

    ! res_n = res_n / res0_n
    cof(1)=norm/(norm0+small)
    if(cof(1)<reduc)exit
  end do
end subroutine




subroutine preconditioner(a_ij,ja,ia,c)
  use constants_m
  use matrix_utils_m
  use vector_utils_m

  integer,dimension(:), intent(in)::ja,ia
  real(rp),dimension(:), intent(in)::a_ij

  real(rp),dimension(:), intent(inout)::c
  real(rp)::summ
  integer::i,kii,j,kij,kji,n

  n=size(ia)-1
  do i=1,n
    kii=csr_index(i,i,ja,ia)
    summ=a_ij(kii)

    do kij=ia(i),ia(i+1)-1
      j=ja(kij)
      if(j<i)then
        kji=csr_index(j,i,ja,ia)
        summ=summ-a_ij(kij)*c(j)*a_ij(kji)
      else
        exit
      end if
    end do

    c(i)=one/(summ+small)
  end do
end subroutine




subroutine direct_solve(a_ij,c,b,ja,ia,x)
  use constants_m
  use matrix_utils_m
  use vector_utils_m

  integer,dimension(:), intent(in)::ja,ia
  real(rp),dimension(:), intent(in)::a_ij
  real(rp),dimension(:), intent(in)::c,b

  real(rp),dimension(:), intent(inout)::x
  real(rp)::cof,summ
  integer::i,kii,j,kij,kji,n

  n=size(ia)-1
  do i=1,n
    summ=b(i)

    do kij=ia(i),ia(i+1)-1
      j=ja(kij)
      if(j<i)then
        summ=summ-a_ij(kij)*x(j)
      else
        exit
      end if
    end do
    x(i)=summ*c(i)
  end do

  do i=1,n
    cof=one/(c(i)+small)
    x(i)=cof*x(i)
  end do

  do i=n,1,-1
    summ=x(i)

    do kij=ia(i+1)-1,ia(i),-1
      j=ja(kij)
      if(j>i)then
        summ=summ-a_ij(kij)*x(j)
      else
        exit
      end if
    end do
    x(i)=summ*c(i)
  end do
end subroutine

end module
