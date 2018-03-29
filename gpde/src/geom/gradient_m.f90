module gradient_m
use const_m
implicit none

real(rp),dimension(:,:,:),allocatable::grad_wk

contains



!==============================================================================!
subroutine ls_grad_matrix(msh,geom)
!------------------------------------------------------------------------------!
  use vector_utils_m
  use matrix_utils_m
  use mesh_data_m
  use geom_data_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(inout)::geom
!------------------------------------------------------------------------------!
  real(rp),dimension(msh%n_dx)::dx
  real(rp),dimension(msh%n_dx,msh%n_dx)::dx_mat
  integer::i,j,k,ib,ic,jc
!------------------------------------------------------------------------------!

  geom%dx_inv=0.d0

  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)

    select case(j)
    case(1:)
      dx=geom%x_vc(:,j)-geom%x_vc(:,i)

      do ic=1,msh%n_dx
        do jc=1,msh%n_dx
          geom%dx_inv(ic,jc,i)=geom%dx_inv(ic,jc,i) + dx(ic)*dx(jc)
          geom%dx_inv(ic,jc,j)=geom%dx_inv(ic,jc,j) + dx(ic)*dx(jc)
        end do
      end do

    case default
      ib=abs(j)
      dx=geom%x_fc(:,k)-geom%x_vc(:,i)

      do ic=1,msh%n_dx
        do jc=1,msh%n_dx
          geom%dx_inv(ic,jc,i)=geom%dx_inv(ic,jc,i) + dx(ic)*dx(jc)
        end do
      end do
    end select
  end do

  do i=1,msh%n_elm
    dx_mat=geom%dx_inv(:,:,i)
    call matrix_inv(dx_mat,geom%dx_inv(:,:,i))
  end do
end subroutine



!==============================================================================!
subroutine ls_gradient(msh,geom,phi,grad,neu_bnd)
!------------------------------------------------------------------------------!
  use vector_utils_m
  use matrix_utils_m
  use mesh_data_m
  use geom_data_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:,:),intent(inout)::phi
  real(rp),dimension(:,:,:),intent(inout)::grad
  logical,intent(in)::neu_bnd
!------------------------------------------------------------------------------!
  real(rp),dimension(size(phi,1))::dphi
  real(rp),dimension(msh%n_dx)::dx
  real(rp),dimension(msh%n_dx,msh%n_dx)::dx_mat
  integer::i,j,k,ib,ic,n_cmp
!------------------------------------------------------------------------------!
  n_cmp=size(phi,1)
  grad=0.d0

  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)

    select case(j)
    case(1:)
      dphi=phi(:,j)-phi(:,i)
      dx=geom%x_vc(:,j)-geom%x_vc(:,i)

      do ic=1,n_cmp
        grad(:,ic,i)=grad(:,ic,i) + dx*dphi(ic)
        grad(:,ic,j)=grad(:,ic,j) + dx*dphi(ic)
      end do

    case default
      ib=abs(j)
      dphi=phi(:,ib)-phi(:,i)
      if(neu_bnd)dphi=0.d0
      dx=geom%x_fc(:,k)-geom%x_vc(:,i)

      do ic=1,n_cmp
        grad(:,ic,i)=grad(:,ic,i) + dx*dphi(ic)
      end do
    end select
  end do

  do i=1,msh%n_elm
    do ic=1,n_cmp
      grad(:,ic,i)=matmul(geom%dx_inv(:,:,i),grad(:,ic,i))
    end do
  end do

  if(neu_bnd)then
    do k=1,msh%n_bfce
      i=msh%fce_elm(k)%lst(1)
      ib=abs(msh%fce_elm(k)%lst(2))
      dx=geom%x_fc(:,k)-geom%x_vc(:,i)

      do ic=1,n_cmp
        phi(ic,ib)=phi(ic,i)+dot_prod(grad(:,ic,i),dx)
      end do
    end do
  end if
end subroutine



!==============================================================================!
subroutine gradient(msh,geom,phi,grad,n_iter,neu_bnd)
!------------------------------------------------------------------------------!
  use vector_utils_m
  use mesh_data_m
  use geom_data_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:,:),intent(inout)::phi
  real(rp),dimension(:,:,:),intent(inout)::grad
  integer,intent(in)::n_iter
  logical,intent(in)::neu_bnd
!------------------------------------------------------------------------------!
!   real(rp),dimension(msh%n_dx,size(phi,1),msh%n_elm)::grad_wk
  real(rp)::w_fp1,w_fn1
  real(rp),dimension(size(phi,1))::fphi,vphi,nphi
  real(rp),dimension(msh%n_dx)::dx,vdx,ndx,norm1,vec
  real(rp),dimension(msh%n_dx,size(phi,1))::fgrad
  integer::i,j,k,ib,ic,iter,n_cmp
!------------------------------------------------------------------------------!
  if(n_iter<1)then
    call ls_gradient(msh,geom,phi,grad,neu_bnd)
    return
  end if

  n_cmp=size(phi,1)
  if(n_cmp>size(grad,2))stop "error calling gradient()"

  k=msh%n_elm
  j=n_cmp
  i=msh%n_dx
  grad_wk(1:i,1:j,1:k)=0.d0

  do iter=1,n_iter
    grad=0.d0


    if(neu_bnd)then
      do k=1,msh%n_bfce
        i=msh%fce_elm(k)%lst(1)
        ib=abs(msh%fce_elm(k)%lst(2))

        vdx=geom%x_fc(:,k)-geom%x_vc(:,i)
        do ic=1,n_cmp
          phi(ic,ib)=phi(ic,i)+dot_prod(grad_wk(:,ic,i),vdx)
        end do
      end do
    end if


    do k=1,msh%n_fce
      i=msh%fce_elm(k)%lst(1)
      j=msh%fce_elm(k)%lst(2)
      norm1=geom%norm(:,k)

      select case(j)
      case(1:)
        w_fp1=geom%w_fp(k)
        w_fn1=1.d0-w_fp1

        do ic=1,n_cmp
          fgrad(:,ic)=w_fp1*grad_wk(:,ic,j)+w_fn1*grad_wk(:,ic,i)
        end do
        dx=w_fp1*geom%x_vc(:,j)+w_fn1*geom%x_vc(:,i)
        dx=geom%x_fc(:,k)-dx
        do ic=1,n_cmp
          fphi(ic)=dot_prod(fgrad(:,ic),dx)
        end do
        fphi=w_fp1*phi(:,j)+w_fn1*phi(:,i)+fphi

        vec=geom%area(k)*norm1*geom%volr(i)
        do ic=1,n_cmp
          grad(:,ic,i)=grad(:,ic,i)+vec*fphi(ic)
        end do

        vec=geom%area(k)*norm1*geom%volr(j)
        do ic=1,n_cmp
          grad(:,ic,j)=grad(:,ic,j)-vec*fphi(ic)
        end do

      case default
        ib=abs(j)
        vec=geom%area(k)*norm1*geom%volr(i)
        do ic=1,n_cmp
          grad(:,ic,i)=grad(:,ic,i)+vec*phi(ic,ib)
        end do
      end select
    end do

    if(iter/=n_iter)then
      k=msh%n_elm
      j=n_cmp
      i=msh%n_dx
      grad_wk(1:i,1:j,1:k)=grad(1:i,1:j,1:k)
    end if
  end do
end subroutine



!==============================================================================!
subroutine test_gradient(msh,geom,pwr)
!------------------------------------------------------------------------------!
  use const_m
  use vector_utils_m
  use mesh_data_m
  use geom_data_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(inout)::msh
  type(geom_t),intent(inout)::geom
  integer,intent(in)::pwr
!------------------------------------------------------------------------------!
  real(rp),dimension(:,:),allocatable::phi
  real(rp),dimension(:,:,:),allocatable::grad
  real(rp)::grad_exact,err,summ
  real(rp),dimension(msh%n_dx)::tmp
  integer::i,k,ix,ib
!------------------------------------------------------------------------------!
  allocate(phi(msh%n_dx,msh%n_ebf)); phi=0
  allocate(grad(msh%n_dx,msh%n_dx,msh%n_elm))


  do i=1,msh%n_elm
    phi(:,i)=geom%x_vc(:,i)**pwr
  end do

  do i=1,msh%n_bfce
    ib=abs(msh%fce_elm(i)%lst(2))
    phi(:,ib)=geom%x_fc(:,i)**pwr
  end do


  print "(1x,a)","gradient error:"
  do k=0,3
    call gradient(msh,geom,phi,grad,k,.false.)
    summ = one - sum (grad)/(msh%n_elm * msh%n_dx)
    print "(a,1x,g16.6)","GG("//trim(itoa(k))//") gradient error:",summ
  end do
end subroutine

end module
