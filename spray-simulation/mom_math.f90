module mom_math_m
use mom_data_m
implicit none
contains



!******************************************************************************

  function fact(n)
  implicit none

  integer,intent(in)::n
  integer::fact
  integer::i


  fact=1

  if(n>1)then
    do i=2,n
      fact=i*fact
    enddo
  endif

  end function fact



!******************************************************************************

  recursive function ncr(n,k)
  implicit none

  integer,intent(in)::n
  integer,intent(in)::k
  integer::ncr
  integer::max_d,i,j,tot_n,tot_d,loc(1)
  integer,dimension(:),allocatable::numer,denom


  ! initialize
  max_d=max(k,n-k)


  ! create numerator and denominator lists
  ncr=1; j=0
  if(k>0)then
    allocate(numer(max_d)); numer=1
    allocate(denom(max_d)); denom=1

    do i=max_d+1,n
      j=j+1
      numer(j)=i
      denom(j)=j
    end do


    ! cyclical reduction of lists
    do i=j,1,-1
      loc=maxloc(numer(:))
      do
        if(mod(numer(loc(1)),denom(i))==0)then
          ncr=ncr*numer(loc(1))/denom(i)
          numer(loc(1))=1
          denom(i)=1
          exit
        else
          loc(1)=loc(1)-1
          if(loc(1)==0)exit
        end if
      end do
    end do


    ! inclusion of what did not reduce
    tot_n=1; tot_d=1
    do i=1,j
      tot_n=tot_n*numer(i)
      tot_d=tot_d*denom(i)
    end do
    deallocate(numer,denom)


    ! result
    ncr=ncr*tot_n
    ncr=ncr/tot_d
  end if

  end function ncr



!******************************************************************************

  recursive function gamma_moment(i,mu_pow,param_v,param_r,param_k,mu) result(mom)
  integer,intent(in)::i
  integer,dimension(:),intent(in),optional::mu_pow
  real(kind=wp),intent(in),optional::param_v,param_r,param_k
  real(kind=wp),dimension(:),intent(in),optional::mu
  real(kind=wp)::mom
  real(kind=wp)::param_r1,param_r2,param_rr
  real(kind=wp)::param_kk,param_tt,number
  integer::j,lb
  real(kind=wp),parameter::k_lb=2._wp
  real(kind=wp),parameter::k_ub=25._wp


  if(present(param_v))then

  if(i==3)then
    mom=param_v/(4._wp/3._wp*pi)
  else
    j=3
    do
      if(j==3)then
        mom=gamma_moment(j,mu_pow,param_v,param_r,param_k)
      end if
      if(i<3)then
        mom=mom/param_r*(param_k+2._wp)/(param_k+j-1)
        j=j-1
      else if(i>3)then
        mom=mom*param_r*(param_k+j)/(param_k+2._wp)
        j=j+1
      end if
      if(j==i)exit
    end do
  end if


  else if(present(mu))then

    ! first moment power
    lb=mu_pow(1)

    ! parameter k
    param_r1=mu(2)/mu(1)
    param_r2=mu(3)/mu(2)
    param_rr=param_r2/param_r1
    param_kk=(lb*(1-param_rr)+1)/(param_rr-1)
    if(param_kk<k_lb)param_kk=k_lb
    if(param_kk>k_ub)param_kk=k_ub

    ! parameter t
    param_tt=param_r1/(param_kk+lb)

    ! number
    number=mu(1)
    do
      if(lb==0)exit
      lb=lb-1
      number=number/(param_tt*(param_kk+lb))
    end do

    ! integration
    mom=number*param_tt**(i) &
      *exp(gammln(param_kk+i)-gammln(param_kk))

!       mom=mom &
!         *(gammp(param_kk+i,r_ub/param_tt) &
!         -gammp(param_kk+i,r_lb/param_tt))
  end if
  end function



!******************************************************************************

    function int_pdf(r,k)

    real(kind=wp),dimension(:,:),intent(in)::r
    integer,intent(in)::k
    real(kind=wp)::int_pdf

    real(kind=wp)::dx,x,y,x_pow_k
    integer::i


    ! integrate
    int_pdf=0._wp

    do i=2,size(r,1)
      dx=r(i,1)-r(i-1,1)
      x=0.5_wp*(r(i-1,1)+r(i,1))
      y=0.5_wp*(r(i-1,2)+r(i,2))

      x_pow_k=x**(k)
      int_pdf=int_pdf+dx*y*x_pow_k
    end do

    end function int_pdf



!******************************************************************************

  subroutine integrate_pdf(num,i1,r,x_b,mu_i)
  implicit none

  real(kind=wp),intent(in)::num
  integer,intent(in)::i1
  real(kind=wp),dimension(:,:),intent(in)::r
  real(kind=wp),dimension(:,:),intent(inout)::x_b
  real(kind=wp),dimension(:,:),intent(inout)::mu_i

  type dat
    type(dat),pointer::before,after,prev,next
    real(kind=wp),dimension(2)::r
    real(kind=wp),dimension(:),allocatable::mu_i
    integer::b,bl,id
  end type dat
  type(dat),pointer::root,temp,prev,next

  integer::i,j,id
  real(kind=wp)::x_lb,x_ub,del_x,dx,x,x_pow_i,dy,y,grad
  logical,dimension(size(x_b,1))::l_summ


  ! initialise
  id=0
  nullify(root,temp,prev,next)


  ! create coordinate list
  do i=1,size(r,1)
    id=id+1

    allocate(temp)
    nullify(temp%before,temp%after)

    temp%r(:)=r(i,:)
    temp%b=0
    temp%bl=0
    temp%id=id
    call add_node(root,temp)
  end do

  x_lb=minval(r(:,1))
  x_ub=maxval(r(:,1))
  del_x=x_ub*ten**(-4)



  ! add boundary limits to the list
  nullify(temp)

  do i=1,size(x_b,1)
    do j=1,size(x_b,2)

      if(x_b(i,j)<x_lb)x_b(i,j)=x_lb+del_x
      if(x_b(i,j)>x_ub)x_b(i,j)=x_ub-del_x

      id=id+1
      allocate(temp)
      nullify(temp%before,temp%after)

      temp%r(1)=x_b(i,j)
      temp%b=i
      temp%bl=j
      temp%id=id
      call add_node(root,temp)

      nullify(prev,next)
      call find_previous(root,temp,prev)
      call find_next(root,temp,next)

      if(associated(prev).and.associated(next))then
        dy=next%r(2)-prev%r(2)
        dx=next%r(1)-prev%r(1)+del_x
        grad=dy/dx

        dx=temp%r(1)-prev%r(1)
        y=prev%r(2)
        temp%r(2)=grad*dx+y
      end if
    end do
  end do


  ! find first node
  temp=>root
  do
    call find_previous(root,temp,prev)
    if(.not.associated(prev))exit
    temp=>prev
  end do


  ! calculate moments
  mu_i=0._wp
  l_summ(:)=.false.

  nullify(prev,next)
  do
    if(.not.associated(prev))then
      call find_previous(root,temp,prev)
    end if

    if(associated(prev))then
      allocate(temp%mu_i(size(mu_i,2)))

      dx=temp%r(1)-prev%r(1)
      x=0.5_wp*(prev%r(1)+temp%r(1))
      y=0.5_wp*(prev%r(2)+temp%r(2))

      do j=1,size(mu_i,2)
        x_pow_i=max(x,del_x)**(i1+j-1)
        temp%mu_i(j)=y*dx*x_pow_i

        do i=1,size(x_b,1)
          if(l_summ(i))then
            mu_i(i,j)=mu_i(i,j)+temp%mu_i(j)
          end if
        end do
      end do

      if(temp%b>0)then
        if(temp%bl==1)l_summ(temp%b)=.true.
        if(temp%bl==2)l_summ(temp%b)=.false.
      end if
    end if

    call find_next(root,temp,next)
    if(.not.associated(next))exit
    prev=>temp
    temp=>next
  end do


  ! multiply by number
  mu_i=num*mu_i

  contains



!------------------------------------------------------------------------------

    recursive subroutine add_node(root,temp)
    type(dat),pointer::root,temp

    if(.not.associated(root))then
      root=>temp

    else if(temp%r(1)<root%r(1))then
      if(associated(root%before))then
        call add_node(root%before,temp)
      else
        root%before=>temp
      end if

    else
      if(associated(root%after))then
        call add_node(root%after,temp)
      else
        root%after=>temp
      end if
    end if

    end subroutine add_node



!------------------------------------------------------------------------------

    recursive subroutine find_before(root,ref,search)
    type(dat),pointer::root,ref,search

    if((ref%r(1)>root%r(1)) &
      .or.((ref%r(1)==root%r(1)).and.(ref%id>root%id)))then
      search=>root

    else
      if(associated(root%before))then
        call find_before(root%before,ref,search)
      else
        nullify(search)
      end if
    end if

    end subroutine find_before



!------------------------------------------------------------------------------

    recursive subroutine find_after(root,ref,search)
    type(dat),pointer::root,ref,search

    if((ref%r(1)<root%r(1)) &
      .or.((ref%r(1)==root%r(1)).and.(ref%id<root%id)))then
      search=>root

    else
      if(associated(root%after))then
        call find_after(root%after,ref,search)
      else
        nullify(search)
      end if
    end if

    end subroutine find_after



!------------------------------------------------------------------------------

    recursive subroutine find_between(root,before,after,between)
    type(dat),pointer::root,before,after,between


    if((before%r(1)<root%r(1).and.after%r(1)>root%r(1)) &
    .or.(before%r(1)==root%r(1).and.before%id<root%id &
    .and.after%r(1)>root%r(1)) &
    .or.(after%r(1)==root%r(1).and.after%id>root%id &
    .and.before%r(1)<root%r(1)) &
    .or.((before%r(1)==root%r(1).and.before%id<root%id) &
    .and.(after%r(1)==root%r(1).and.after%id>root%id)))then
      between=>root
      if(associated(between,before).or.associated(between,after))then
        nullify(between)
        return
      end if

    else
      if(before%r(1)>=root%r(1))then
        if(associated(root%after))then
          call find_between(root%after,before,after,between)
        end if

      else if(after%r(1)<=root%r(1))then
        if(associated(root%before))then
          call find_between(root%before,before,after,between)
        end if
      end if
    end if

    end subroutine find_between



!------------------------------------------------------------------------------

    subroutine find_previous(root,ref,search)
    type(dat),pointer::root,ref,search, between


    nullify(search)
    call find_before(root,ref,search)

    if(.not.associated(search))return
    do
      nullify(between)
      call find_between(root,search,ref,between)
      if(.not.associated(between))return
      search=>between
    end do

    end subroutine find_previous



!------------------------------------------------------------------------------

    subroutine find_next(root,ref,search)
    type(dat),pointer::root,ref,search, between


    nullify(search)
    call find_after(root,ref,search)

    if(.not.associated(search))return
    do
      nullify(between)
      call find_between(root,ref,search,between)
      if(.not.associated(between))return
      search=>between
    end do

    end subroutine find_next
  end subroutine integrate_pdf



!******************************************************************************

  subroutine product_difference(mu,w_mu,r_mu)
!   use f95_lapack
  implicit none

  real(kind=wp),dimension(:),intent(inout)::mu
  real(kind=wp),dimension(:),intent(inout)::w_mu
  real(kind=wp),dimension(:),intent(inout)::r_mu

  real(kind=wp),dimension(:,:),allocatable::p
  real(kind=wp),dimension(:),allocatable::alpha
  real(kind=wp),dimension(:,:),allocatable::a
  real(kind=wp),dimension(:),allocatable::w
  real(kind=wp)::delta
  integer::i,j,j_max,n,lb,ub
  logical::update,mod_mu


  ! initialize
  allocate(p(size(mu,1),size(mu,1)+1))
  allocate(alpha(size(mu,1)))

  p(:,1)=0._wp
  p(1,1)=1._wp


  ! matrix p(i,j)
  j=2
  do i=1,size(p,1)
    p(i,j)=(-1)**(i-1)*mu(i)
  end do

  do j=3,size(p,2)
    do i=1,size(p,1)-(j-2)
      p(i,j)=p(1,j-1)*p(i+1,j-2)-p(1,j-2)*p(i+1,j-1)
    end do
  end do


  ! vector alpha(n)
  alpha(1)=zero
  do n=2,size(p,1)
    alpha(n)=p(1,n+1)/(p(1,n)*p(1,n-1))
  end do


  ! vector a(i)
  allocate(a(size(mu,1)/2,size(mu,1)/2))
  a=0._wp

  do n=1,size(a,1)
    a(n,n)=alpha(2*n)+alpha(2*n-1)
  end do


  ! vector b(i)
  allocate(w(size(mu,1)/2))

  do n=1,size(a,1)-1
    a(n,n+1)=sqrt(abs(alpha(2*n+1)*alpha(2*n)))
  end do

  deallocate(alpha)


  ! solve abscissas and weights
!   call la_syevd(a, w, 'V', 'U') !!requires:f90_lapack

  w_mu(:)=a(1,:)**2*mu(1)
  r_mu(:)=w(:)

  end subroutine product_difference



!******************************************************************************

  function mu_quadr(w_mu,r_mu,pow)
  implicit none

  real(kind=wp),dimension(:),intent(in)::w_mu
  real(kind=wp),dimension(:),intent(in)::r_mu
  real(kind=wp),intent(in)::pow
  real(kind=wp)::mu_quadr


  ! moment using quadrature
  mu_quadr=sum(w_mu(:)*r_mu(:)**(pow))

  end function mu_quadr



!******************************************************************************

  function mu_interp(mu_i,ista,pow)
  implicit none

  real(kind=wp),dimension(ista:),intent(in)::mu_i
  integer,intent(in)::ista
  real(kind=wp),intent(in)::pow
  real(kind=wp)::mu_interp

  integer::l,u,n
  real(kind=wp)::delta


  ! moment locations
  n=nint(pow)
  l=floor(pow)
  u=ceiling(pow)


  ! difference
  delta=pow-l


  ! interpolated moment
  mu_interp=mu_i(l)**(1._wp-delta)*mu_i(u)**delta

  end function mu_interp



!******************************************************************************

  function gammp(a,x)
  implicit none

  ! incomplete gamma function p(a,x)
  ! [num. res. in f77, ch 6.2]
  ! given a and x, returns the incomplete gamma function p(a,x).
  ! p(a,0.)=0., p(a,infinity)=1.

  real(kind=wp)::a,x
  real(kind=wp)::gammp,gammcf,gamser,gln


  if(x<a+1)then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=one-gammcf
  endif

  if(gammp<small)gammp=zero

  end function gammp



!******************************************************************************

  function gammq(a,x)
  implicit none

  ! incomplete gamma function q(a,x)
  ! [num. res. in f77, ch 6.2]
  ! given a and x, returns the incomplete gamma function q(a,x).
  ! q(a,x)=1-p(a,x)
  ! q(a,0.)=1., q(a,infinity)=0.

  real(kind=wp)::a,x
  real(kind=wp)::gammq,gammcf,gamser,gln


  if(x<a+1)then
    call gser(gamser,a,x,gln)
    gammq=one-gamser
  else
    call gcf(gammcf,a,x,gln)
    gammq=gammcf
  endif

  if(gammq<small)gammq=zero

  end function gammq



!******************************************************************************

  subroutine gser(gamser,a,x,gln)
  implicit none

!   series representation of the incomplete gamma function
!   [num. res. in f77, ch 6.2]
! returns the incomplete gamma function p(a,x)

  real(kind=wp)::gamser,a,x,gln
  real(kind=wp)::ap,del,summ!,gammln
  integer::n
  integer,parameter::itmax=10
  real(kind=wp),parameter::eps=3.e-7


  gln=gammln(a)

  if(x<small)then
    gamser=zero
    return
  endif

  ap=a
  summ=one/a
  del=summ

  do n=1,itmax
    ap=ap+one
    del=del*x/ap
    summ=summ+del
    if(abs(del)<abs(summ)*eps)exit
  enddo

  gamser=summ*exp(-x+a*log(x)-gln)
  if(gamser<small)gamser=zero

  end subroutine gser



!******************************************************************************

  subroutine gcf(gammcf,a,x,gln)
  implicit none

  ! continued fraction representation of the incomplete gamma function
  ! [num. res. in f77, ch 6.2]
  ! returns the incomplete gamma function q(a,x)

  real(kind=wp)::gammcf,a,x,gln
  real(kind=wp)::an,b,c,d,del,h!,gammln
  integer::i

  integer,parameter::itmax=10
  real(kind=wp),parameter::eps=3.e-7
  real(kind=wp),parameter::fpmin=1.e-30


  gln=gammln(a)
  b=x+one-a
  c=one/fpmin
  d=one/b
  h=d

  do i=1,itmax
    an=-i*(i-a)
    b=b+2.
    d=an*d+b
    if(abs(d)<fpmin)d=fpmin
    c=b+an/c
    if(abs(c)<fpmin)c=fpmin
    d=one/d
    del=d*c
    h=h*del
    if(abs(del-one)<eps)exit
  enddo

  gammcf=exp(-x+a*log(x)-gln)*h

  end subroutine gcf



!******************************************************************************

  function gammln(xx)
  implicit none

  ! natural logarithm of the gamma function
  ! [num. res. in f77, ch 6.1]
  ! given xx, where xx > 0, ln(g(xx)) is returned

  real(kind=wp)::gammln
  real(kind=wp)::xx
  integer::j
  real(kind=wp)::ser,tmp,x,y
  real(kind=wp),save::stp=zero
  real(kind=wp),dimension(1:6),save::cof


  if(stp<small)then
    stp=2.5066282746310005_wp
    cof(1)=76.18009172947146_wp
    cof(2)=-86.50532032941677_wp
    cof(3)=24.01409824083091_wp
    cof(4)=-1.231739572450155_wp
    cof(5)=0.1208650973866179_wp*ten**(-2)
    cof(6)=-0.5395239384953_wp*ten**(-5)
  end if

  if(xx<ten**(-2))then
    xx=ten**(-2)
  end if

  x=xx
  y=x
  tmp=x+5.5_wp
  tmp=(x+0.5_wp)*log(tmp)-tmp
  ser=1.0000000001900_wp

  do j=1,6
    y=y+one
    ser=ser+cof(j)/y
  enddo

  gammln=tmp+log(stp*ser/x)

  end function gammln



!******************************************************************************

  subroutine loc_x_ord(x_lb,x_ub,x,y,c)
  implicit none

  real(kind=wp),intent(in)::x_lb
  real(kind=wp),intent(in)::x_ub
  real(kind=wp),dimension(:),intent(inout)::x
  real(kind=wp),dimension(:),intent(in),optional::y
  real(kind=wp),intent(in),optional::c

  integer::i,num_pts
  real(kind=wp)::cof,dx,grad,const,r(2,2)


  num_pts=size(x,1)
  dx=(x_ub-x_lb)/(num_pts-1)

  x(1)=x_lb

  do i=2,num_pts
    x(i)=x(i-1)+dx
  end do

  end subroutine loc_x_ord

!******************************************************************************

end module
