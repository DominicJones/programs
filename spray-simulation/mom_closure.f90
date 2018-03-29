module mom_closure_m
use mom_math_m
implicit none
contains



!******************************************************************************

  subroutine maxent_pdf(mu,mu_pow,r,str,r_ub, &
               x_lb_opt,x_ub_opt,r_norm_opt,urf_opt,ip_opt)
!   use f95_lapack !! requires: f95_lapack
  implicit none

  real(kind=wp),dimension(:),intent(in)::mu
  integer,dimension(:),intent(in)::mu_pow
  real(kind=wp),dimension(:,:),intent(inout)::r
  integer,intent(in)::str
  real(kind=wp),dimension(:,:),intent(inout),optional::r_ub
  real(kind=wp),intent(in),optional::x_lb_opt
  real(kind=wp),intent(in),optional::x_ub_opt
  real(kind=wp),intent(in),optional::r_norm_opt
  real(kind=wp),intent(in),optional::urf_opt
  integer,intent(in),optional::ip_opt

  real(kind=wp),dimension(size(mu,1))::lambda
  real(kind=wp),dimension(12),save::lambda_s
  real(kind=wp),dimension(size(mu,1))::mu_approx
  real(kind=wp),dimension(size(mu,1),size(mu,1))::a
  real(kind=wp),dimension(size(mu,1),1)::b

  real(kind=wp)::r_norm,urf
  real(kind=wp)::x_lb,x_ub
  real(kind=wp)::r_lambda,tol_lambda=0.0001_wp
  integer::i,j,mu_ij,ip
  integer::it,it_max=50


  ! lower x-ordinate
  x_lb=0._wp
  if(present(x_lb_opt))x_lb=x_lb_opt


  ! upper x-ordinate
  x_ub=1._wp
  if(present(x_ub_opt))x_ub=x_ub_opt


  ! mean
  r_norm=1._wp
  if(present(r_norm_opt))r_norm=r_norm_opt


  ! urf
  urf=0.5_wp
  if(present(urf_opt))urf=urf_opt


  ! p(x)
  ip=5
  if(present(ip_opt))ip=ip_opt




  ! initialize data
  lambda=0._wp


  ! calculate x-ordinates
  call loc_x_ord(x_lb,x_ub,r(:,1))


  ! iterate
  do it=1,it_max


    ! calculate y-ordinates
    do i=1,size(r,1)
      r(i,2)=pdf(r(i,1),lambda,x_lb,x_ub,r_norm,ip)
    end do


    ! solve b
    do j=1,size(mu,1)
      mu_ij=mu_pow(j)
      b(j,1)=int_pdf(r,mu_ij)
      b(j,1)=mu(j)-b(j,1)
    end do


    ! solve a
    do i=1,size(mu,1)
      do j=i,size(mu,1)
        mu_ij=(i-1)+mu_pow(j)
        a(i,j)=-int_pdf(r,mu_ij)
      end do
    end do


    ! solve x
!     call la_sysv(a,b) !! requires: f95_lapack


    ! update lambda
    do i=1,size(mu,1)
      lambda(i)=lambda(i)+urf*b(i,1)
    end do


    ! check convergence
    b(:,1)=abs(b(:,1)/(lambda(:)+small))
    r_lambda=maxval(b(:,1))

!     print*,'max ent:',it,r_lambda
    if(r_lambda<tol_lambda)exit
  end do


  ! store lambda
!   if(str/=0)then
!     lambda_s(1:size(mu,1))=lambda
!   end if


  ! line search y-ordinates
  if(present(r_ub))then
    do i=1,size(r_ub,1)
      r_ub(i,2)=pdf(r_ub(i,1),lambda,x_lb,x_ub,r_norm,ip)
    end do
  end if

  contains



!------------------------------------------------------------------------------

    function pdf(x,lambda,x_lb,x_ub,r_norm,ip)

    real(kind=wp),intent(in)::x
    real(kind=wp),dimension(:),intent(in)::lambda
    real(kind=wp),intent(in)::x_lb
    real(kind=wp),intent(in)::x_ub
    real(kind=wp),intent(in)::r_norm
    integer,intent(in)::ip

    real(kind=wp)::pdf

    real(kind=wp)::s,p_x,x_pow_k,exp_sum
    integer::k


    ! p(x)
    if(ip==1)p_x=1._wp
    if(ip==2)p_x=1._wp/(x_ub-x_lb)
    if(ip==3)p_x=exp(-x/r_norm)/r_norm
    if(ip==4)then
      s=r_norm*sqrt(2._wp/pi)
      p_x=exp(-0.5_wp*(x/s)**2)*x/s**2
    end if
    if(ip==5)then
      p_x=16._wp*x/r_norm**2*exp(-4._wp*x/r_norm)
    end if


    ! exp[Sum]
    x_pow_k=1._wp
    exp_sum=lambda(1)

    do k=2,size(lambda,1)
      x_pow_k=x_pow_k*x
      exp_sum=exp_sum+lambda(k)*x_pow_k
    end do
    exp_sum=exp(-exp_sum)


    ! pdf
    pdf=p_x*exp_sum

    end function pdf
  end subroutine maxent_pdf



!******************************************************************************

  subroutine splines_pdf(mu,mu_pow,r,r_ub, &
               x_lb_opt,x_ub_opt,deg_spl_opt, &
               omega_opt,linear_b_opt,rcond_opt,tol_s_opt)
!   use f95_lapack !! requires: f95_lapack
  implicit none

  real(kind=wp),dimension(:),intent(in)::mu
  integer,dimension(:),intent(in)::mu_pow
  real(kind=wp),dimension(:,:),intent(inout)::r
  real(kind=wp),dimension(:,:),intent(inout),optional::r_ub
  real(kind=wp),intent(in),optional::x_lb_opt
  real(kind=wp),intent(in),optional::x_ub_opt
  integer,intent(in),optional::deg_spl_opt
  real(kind=wp),intent(in),optional::omega_opt
  integer,intent(in),optional::linear_b_opt
  real(kind=wp),intent(in),optional::rcond_opt
  real(kind=wp),intent(in),optional::tol_s_opt

  real(kind=wp),dimension(:,:),allocatable::r_spl
  real(kind=wp),dimension(:,:),allocatable::a
  real(kind=wp),dimension(:,:),allocatable::b
  real(kind=wp),dimension(:),allocatable::s

  real(kind=wp),dimension(:,:),allocatable::svd_u
  real(kind=wp),dimension(:),allocatable::svd_s
  real(kind=wp),dimension(:,:),allocatable::svd_vt
  real(kind=wp),dimension(:,:),allocatable::a_copy
  real(kind=wp)::grad_sigma,grad_tp

  real(kind=wp)::x_lb,x_ub
  real(kind=wp)::y_lb,y_ub
  integer::deg_spl
  real(kind=wp)::omega
  integer::linear_b
  real(kind=wp)::rcond
  real(kind=wp)::tol_s

  real(kind=wp)::dx,y_min,y_max
  real(kind=wp),dimension(20),save::r_tp,d_tp
  real(kind=wp),dimension(size(r,1))::grad
  real(kind=wp),dimension(:),allocatable::cof
  integer,dimension(:),allocatable::pow

  integer,save::num_tp
  integer::i,j,iplus,jplus,k,mu_i,pow_mu,loc(1)
  integer::num_spl,num_pts,num_cof


  ! lower x-ordinate
  x_lb=zero
  if(present(x_lb_opt))x_lb=x_lb_opt


  ! upper x-ordinate
  x_ub=one
  if(present(x_ub_opt))x_ub=x_ub_opt


  ! polynomial degree of splines
  deg_spl=3
  if(present(deg_spl_opt))deg_spl=deg_spl_opt


  ! weakening factor
  omega=ten**(-4)
  if(present(omega_opt))omega=omega_opt


  ! linear boundary conditions
  ! [0=neither, 1=lb, 2=ub, 3=both]
  linear_b=3
  if(present(linear_b_opt))linear_b=linear_b_opt


  ! singular value cutoff
  ! [s_min=rcond*s_max]
  rcond=ten**(-12)
  if(present(rcond_opt))rcond=rcond_opt


  ! maximium difference in singular values
  tol_s=0.99_wp
  if(present(tol_s_opt))tol_s=tol_s_opt




  ! number of splines
  num_spl=size(mu,1)+deg_spl


  ! number of points
  num_pts=num_spl+1
  allocate(r_spl(num_pts,2))
  r_spl=zero


  ! number of coefficients (s_ij)
  num_cof=(deg_spl+1)*num_spl


  ! linear algebra arrays
  allocate(a(num_cof,num_cof)); a=zero
  allocate(b(num_cof,1)); b=zero
  allocate(s(num_cof)); s=zero


  ! working arrays
  allocate(cof(0:deg_spl)); cof=zero
  allocate(pow(0:deg_spl)); pow=zero


  ! spline x-ordinates
  call loc_x_ord(x_lb,x_ub,r_spl(:,1))




  ! elements of matrix A
  j=1
  dx=zero
  do i=1,num_cof-size(mu,1),deg_spl

    do k=0,deg_spl
      cof(k)=one
      pow(k)=k
    end do

    k=(i-1)/deg_spl+1
    if(i>1)then
      dx=r_spl(k,1)-r_spl(k-1,1)
    end if

!     if(k==1)b(i,1)=y_lb
!     if(k==num_pts)b(i,1)=y_ub

    do iplus=0,deg_spl-1
      do jplus=0,deg_spl+iplus+1
        if(jplus<=deg_spl)then
          a(i+iplus,j+jplus)=cof(jplus)*dx**pow(jplus)*omega**(iplus)
        else
          if(k>1.and.k<num_pts)then
            a(i+iplus,j+jplus)=-a(i+iplus,j+jplus-deg_spl-1)
          end if
        end if
      end do

      do jplus=0,deg_spl
        cof(jplus)=pow(jplus)*cof(jplus)
        pow(jplus)=max(0,pow(jplus)-1)
      end do
    end do


    ! linear boundaries
    if((k==1.and.(linear_b==1.or.linear_b==3)) &
    .or.(k==num_pts.and.(linear_b==2.or.linear_b==3)))then
      do iplus=0,deg_spl-1
        do jplus=0,deg_spl
          if(jplus>=2)a(i+iplus,j+jplus)=zero
          if(iplus==1.and.jplus==1)a(i+iplus,j+jplus)=zero
          if(iplus==1.and.jplus==2)a(i+iplus,j+jplus)=omega**(iplus)
          if(iplus==2.and.jplus==3)a(i+iplus,j+jplus)=omega**(iplus)
        end do
      end do
    end if

    if(i/=1)j=j+deg_spl+1
  end do




  ! moment equations
  mu_i=1
  do i=num_cof-size(mu,1)+1,num_cof

    k=0
    do j=1,num_spl
      do jplus=0,deg_spl

        k=k+1
        a(i,k)=int_mu(mu_pow(mu_i),jplus,-r_spl(j,1),r_spl(j,1),r_spl(j+1,1))
      end do
    end do

    b(i,1)=mu(mu_i)
    mu_i=mu_i+1
  end do




  ! least-squares solution using svd
  if(present(rcond_opt))then
!     call la_gelsd(a,b,rcond=rcond) !! requires: f95_lapack


  ! singular value decomposition and solution
  else
    allocate(svd_s(size(a,1)))
    allocate(svd_vt(size(a,1),size(a,1)))
!     call la_gesdd(a=a,s=svd_s,vt=svd_vt,job='U') !! requires: f95_lapack

!     open(101,file='svd_0.dat',status='unknown')
!     do i=1,size(a,1)
!       write(101,'(i4,e15.5)')i,svd_s(i)
!     end do
!     close(101)

    do i=1,size(a,1)-1
      grad_sigma=log10(svd_s(i))-log10(svd_s(i+1))
      j=i; if(grad_sigma>tol_s.and.i>0.6_wp*size(a,1))exit
    end do

!     open(101,file='svd_1.dat',status='unknown')
!     do i=1,j
!       write(101,'(i4,e15.5)')i,svd_s(i)
!     end do
!     close(101)

    svd_s(:j)=1._wp/svd_s(:j)
    svd_s(j+1:)=0._wp

    do i=1,size(a,1)
      a(:,i)=a(:,i)*svd_s(i)
    end do

    a=matmul(transpose(svd_vt),transpose(a))
    b(:,1)=matmul(a,b(:,1))
  end if


  ! pdf x-ordinates
  call loc_x_ord(x_lb,x_ub,r(:,1))


  ! pdf y-ordinates
  do i=1,size(r,1)
    r(i,2)=pdf(deg_spl,b(:,1),r_spl(:,1),r(i,1))
  end do


  ! line search y-ordinates
  if(present(r_ub))then
    do i=1,size(r_ub,1)
      r_ub(i,2)=pdf(deg_spl,b(:,1),r_spl(:,1),r_ub(i,1))
    end do
  end if

  contains



!------------------------------------------------------------------------------

    function pdf(k,s_ij,x_i,x)
    implicit none

    integer,intent(in)::k
    real(kind=wp),dimension(:),intent(in)::s_ij
    real(kind=wp),dimension(:),intent(in)::x_i
    real(kind=wp),intent(in)::x
    real(kind=wp)::pdf

    integer::i,ii,j,ij
    real(kind=wp)::dx


    dx=zero
    do ii=ubound(x_i,1)-1,1,-1
      if(x>=x_i(ii))then
        dx=(x-x_i(ii))
        i=ii
        exit
      end if
    enddo

    pdf=zero
    do j=0,k
      ij=(k+1)*i+j-k
      if(j==0)then
        pdf=pdf+s_ij(ij)
      else
        pdf=pdf+s_ij(ij)*dx**j
      end if
    end do

    end function pdf



!------------------------------------------------------------------------------

    function int_mu(l,n,const,x_lb,x_ub)
    implicit none

    integer,intent(in)::l
    integer,intent(in)::n
    real(kind=wp),intent(in)::const
    real(kind=wp),intent(in)::x_lb
    real(kind=wp),intent(in)::x_ub
    real(kind=wp)::int_mu

    integer::k


    int_mu=zero
    do k=0,n
      int_mu=int_mu &
        +ncr(n,k)*const**(n-k)/(k+l+1) &
        *(x_ub**(k+l+1)-x_lb**(k+l+1))
    end do

    end function int_mu
  end subroutine splines_pdf



!******************************************************************************

  subroutine legendre_pdf(mu,r,r_ub)
  implicit none

  real(kind=wp),dimension(:),intent(in)::mu
  real(kind=wp),dimension(:,:),intent(inout)::r
  real(kind=wp),dimension(:,:),intent(inout),optional::r_ub

  integer::i
  real(kind=wp)::x_lb,x_ub,y_max


  ! lower bound
  x_lb=0._wp


  ! upper bound
  x_ub=1._wp


  ! pdf x-ordinates
  call loc_x_ord(x_lb,x_ub,r(:,1))


  ! pdf y-ordinates
  do i=1,ubound(r,1)
    r(i,2)=pdf(mu,r(i,1))
  end do


  ! line search y-ordinates
  if(present(r_ub))then
    do i=1,size(r_ub,1)
      r_ub(i,2)=pdf(mu,r_ub(i,1))
    end do
  end if

  contains



!------------------------------------------------------------------------------

    function pdf(mu,x)
    implicit none


    real(kind=wp),dimension(:),intent(in)::mu
    real(kind=wp),intent(in)::x
    real(kind=wp)::pdf

    integer::i


    pdf=0._wp

    do i=0,size(mu,1)-1
      pdf=pdf+(2*i+1)*lambda(i,mu)*lambda(i,x=x)
    end do

    end function pdf



!------------------------------------------------------------------------------

    function lambda(i,mu,x)
    implicit none

    integer,intent(in)::i
    real(kind=wp),dimension(:),intent(in),optional::mu
    real(kind=wp),intent(in),optional::x
    real(kind=wp)::lambda

    integer::j
    real(kind=wp)::term


    lambda=0._wp

    do j=0,i
      if(present(mu))then
        term=mu(j+1)
      else
        term=(x+small)**(j)
      end if

      lambda=lambda+ &
        (-1)**(j)*ncr(i,j)*ncr(i+j,j)*term
    end do

    end function lambda
  end subroutine legendre_pdf



!******************************************************************************

  subroutine laguerre_pdf(mu,r,r_ub,x_lb_opt,x_ub_opt,k_opt)
  implicit none

  real(kind=wp),dimension(:),intent(in)::mu
  real(kind=wp),dimension(:,:),intent(inout)::r
  real(kind=wp),dimension(:,:),intent(inout),optional::r_ub
  real(kind=wp),intent(in),optional::x_lb_opt
  real(kind=wp),intent(in),optional::x_ub_opt
  integer,intent(in),optional::k_opt

  integer::i,k
  real(kind=wp)::x_lb,x_ub,y_max


  ! lower bound
  x_lb=0._wp
  if(present(x_lb_opt))x_lb=x_lb_opt


  ! upper bound
  x_ub=1._wp
  if(present(x_ub_opt))x_ub=x_ub_opt


  ! parameter k
  k=0
  if(present(k_opt))k=k_opt


  ! pdf x-ordinates
  call loc_x_ord(x_lb,x_ub,r(:,1))


  ! pdf y-ordinates
  do i=1,ubound(r,1)
    r(i,2)=pdf(mu,r(i,1),k)
  end do


  ! line search y-ordinates
  if(present(r_ub))then
    do i=1,size(r_ub,1)
      r_ub(i,2)=pdf(mu,r_ub(i,1),k)
    end do
  end if

  contains



!------------------------------------------------------------------------------

    function pdf(mu,x,k)
    implicit none


    real(kind=wp),dimension(:),intent(in)::mu
    real(kind=wp),intent(in)::x
    integer,intent(in)::k
    real(kind=wp)::pdf

    integer::i,j,l


    pdf=0._wp

    do i=0,size(mu,1)-1
      l=1
      do j=i,i+k
        l=l*j
      end do
      l=max(1,l)

      pdf=pdf+lambda(i,k,mu)*lambda(i,k,x=x)*l
    end do

    pdf=x**(k)*exp(-x)*pdf

    end function pdf



!------------------------------------------------------------------------------

    function lambda(i,k,mu,x)
    implicit none

    integer,intent(in)::i
    integer,intent(in)::k
    real(kind=wp),dimension(:),intent(in),optional::mu
    real(kind=wp),intent(in),optional::x
    real(kind=wp)::lambda

    integer::j
    real(kind=wp)::term


    lambda=0._wp

    do j=0,i
      if(present(mu))then
        term=mu(j+1)
      else
        term=(x+small)**(j)
      end if

      lambda=lambda+ &
        (-1)**(j)*ncr(i+k,i-j)*term/fact(j)
    end do

    end function lambda
  end subroutine laguerre_pdf



!******************************************************************************

  subroutine gamma_pdf(mu,mu_pow,r,r_ub,x_lb_opt,x_ub_opt,num_opt)
  implicit none

  real(kind=wp),dimension(:),intent(in)::mu
  integer,dimension(:),intent(in)::mu_pow
  real(kind=wp),dimension(:,:),intent(inout)::r
  real(kind=wp),dimension(:,:),intent(inout),optional::r_ub
  real(kind=wp),intent(in),optional::x_lb_opt
  real(kind=wp),intent(in),optional::x_ub_opt
  real(kind=wp),intent(out),optional::num_opt

  real(kind=wp)::par_r1,par_r2,par_rr
  real(kind=wp)::par_k,par_t,num,cof
  real(kind=wp)::x_lb,x_ub
  integer::imu1,i

  real(kind=wp),parameter::k_lb=2._wp
  real(kind=wp),parameter::k_ub=25._wp


  ! order of first moment
  imu1=mu_pow(1)


  ! parameter k
  par_r1=mu(2)/mu(1)
  par_r2=mu(3)/mu(2)
  par_rr=par_r2/par_r1
  par_k=(imu1*(1-par_rr)+1)/(par_rr-1)

  if(par_k<k_lb)par_k=k_lb
  if(par_k>k_ub)par_k=k_ub


  ! parameter t
  par_t=par_r1/(par_k+imu1)


  ! coefficient
  cof=par_t**(par_k)*exp(gammln(par_k))
  cof=1._wp/cof




  ! number
  if(present(num_opt))then
    num_opt=mu(1)
    do
      if(imu1==0)exit
      imu1=imu1-1
      num_opt=num_opt/(par_t*(par_k+imu1))
    end do
  end if


  ! lower bound
  x_lb=0._wp
  if(present(x_lb_opt))x_lb=x_lb_opt


  ! upper bound
  x_ub=1._wp
  if(present(x_ub_opt))x_ub=x_ub_opt


  ! pdf x-ordinates
  call loc_x_ord(x_lb,x_ub,r(:,1))


  ! pdf y-ordinates
  do i=1,ubound(r,1)
    r(i,2)=pdf(par_k,par_t,cof,r(i,1))
  end do


  ! line search y-ordinates
  if(present(r_ub))then
    do i=1,size(r_ub,1)
      r_ub(i,2)=pdf(par_k,par_t,cof,r_ub(i,1))
    end do
  end if

  contains



!------------------------------------------------------------------------------

    function pdf(par_k,par_t,cof,x)
    implicit none

    real(kind=wp),intent(in)::par_k
    real(kind=wp),intent(in)::par_t
    real(kind=wp),intent(in)::cof
    real(kind=wp),intent(in)::x
    real(kind=wp)::pdf


    pdf=cof*(x+small)**(par_k-1._wp)*exp(-x/par_t)

    end function pdf
  end subroutine gamma_pdf



!******************************************************************************

end module
