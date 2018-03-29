module mom_solve_m
use mom_closure_m
implicit none
contains



!******************************************************************************

  recursive subroutine solve_pdf(mu,mu_pow,num,r_norm,r, &
                         r_max_opt,mth_opt,it_max_opt)
  implicit none

  real(kind=wp),dimension(:),intent(in)::mu
  integer,dimension(:),intent(in)::mu_pow
  real(kind=wp),intent(inout)::num
  real(kind=wp),intent(inout)::r_norm
  real(kind=wp),dimension(:,:),intent(inout)::r
  real(kind=wp),intent(in),optional::r_max_opt
  integer,intent(in),optional::mth_opt
  integer,intent(in),optional::it_max_opt

  real(kind=wp)::r_max
  integer::mth
  integer::it_max

  real(kind=wp)::r_max_c
  real(kind=wp),dimension(:),allocatable::r_max_s,y_rat_s
  real(kind=wp),dimension(size(mu,1))::mu_norm

  integer::it,i,j,k,loc(3)
  integer::i_lb,i_mb,i_ub
  real(kind=wp)::y_min,y_max,y_rat,tol_y=0.0001_wp


  ! closure method
  mth=5
  if(present(mth_opt))mth=mth_opt


  ! number of iterations
  it_max=25
  if(present(it_max_opt))it_max=it_max_opt


  ! upper bound
  r_max=3.5_wp*r_norm
  if(present(r_max_opt))r_max=r_max_opt


  ! number
  if(mth==5)then
    num=gamma_moment(0,mu=mu,mu_pow=mu_pow)
  end if


  ! for iterative solution
  if(it_max>0)then

    ! store r_max
    allocate(r_max_s(it_max))
    r_max_s=0._wp


    ! store y_rat
    allocate(y_rat_s(it_max))
    y_rat_s=0._wp
  end if




  ! iterate
  loc=0_wp
  do it=1,max(1,it_max)


    ! normalize moments
    do i=1,size(mu,1)
      mu_norm(i)=mu(i)/(r_max**(mu_pow(i))*num)
    end do


    ! store r_max
    if(it_max>0)r_max_s(it)=r_max




    ! method
    if(mth==1)then
      call maxent_pdf(mu_norm,mu_pow,r,str=it_max,r_norm_opt=r_norm/r_max)

    else if(mth==2)then
      call splines_pdf(mu_norm,mu_pow,r,deg_spl_opt=3)

    else if(mth==3)then
      call legendre_pdf(mu_norm,r)

    else if(mth==4)then
      call laguerre_pdf(mu_norm,r)

    else if(mth==5)then
      call gamma_pdf(mu_norm,mu_pow,r)
    end if




    ! details
    if(it_max>0)then
      loc(1:1)=maxloc(r(:,1),r(:,1)<0.10_wp); i_lb=loc(1)
      loc(2:2)=minloc(r(:,1),r(:,1)>0.50_wp); i_mb=loc(2)
      loc(3:3)=minloc(r(:,1),r(:,1)>0.95_wp); i_ub=loc(3)

      loc(1:1)=minloc(r(i_lb:i_ub,2))+i_lb-1
      y_min=r(loc(1),2)

      loc(2:2)=maxloc(r(i_lb:i_ub,2))+i_lb-1
      y_max=r(loc(2),2)

      y_rat=y_min/y_max
      y_rat_s(it)=y_rat

      loc(3:3)=minloc(abs(y_rat_s(:it)))
      if(loc(3)==it)loc(3)=0

!       print*,'it:',it,r_max,y_min,y_max,y_rat,loc(1:3)
    else
      exit
    end if



    ! update upper bound
    if(mth==1.or.mth==5)then
      loc(3:3)=minloc(r(i_mb:i_ub,2),r(i_mb:i_ub,2)>tol_y*y_max)+i_mb-1
      r_max=r(loc(3),1)*r_max
      exit

    else
      if(it<it_max.and.y_rat<-tol_y.and.abs(y_rat)>tol_y.and.abs(y_rat)<1.1_wp)then
        r_max=0.95_wp*r_max
      else
        if(loc(3)==0.or.y_rat>0._wp)loc(3)=it
!         print*,'exit:',it,loc(3),r_max_s(loc(3))

        r_max=r_max_s(loc(3))
        exit
      end if
    end if

  end do




  ! solve with optimum r_max
  if(it_max/=0)then
!     print*,'optimum r_max:',r_max
    call solve_pdf(mu,mu_pow,num,r_norm,r, &
           r_max_opt=r_max,it_max_opt=0,mth_opt=mth)
  end if


  ! true pdf
  if(it_max==0)then
    do i=1,size(r,1)
      r(i,1)=r(i,1)*r_max
      r(i,2)=r(i,2)/r_max
      r(i,2)=max(0._wp,r(i,2))
    end do
  end if

  end subroutine solve_pdf

!******************************************************************************

end module
