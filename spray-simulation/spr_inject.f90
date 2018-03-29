module spr_inject_m
use mom_solve_m
use grd_data_m
use spr_data_m
use cfd_constr_m
implicit none
contains
!
! injection inlet conditions: moments / velocities
!
  subroutine spr_inject(igrd,ifce,ivol,iph,ieq)
    integer,intent(in)::igrd,ifce,ivol,iph,ieq
    integer::i,j,k,lb,ub
    real(wp)::ref_ts,ref_tf,ref_dt,ref_pow
    real(wp)::pr(4),vel(n_dx)
    real(wp)::r_centre,grad,r_outer,r_inner
    real(wp),save::mom(5)=zero
    type(phs_t),dimension(:),pointer::phs_p


    ! initialise
    phs_p=>grd(igrd)%fce(ifce)%phs
    call cp_id(igrd,iph,ieq,ist=inj_st,lb=lb,ub=ub)


    ! orifice radii
    r_outer=r_orif
    r_inner=r_orif-sheet_th
    r_centre=grd(igrd)%fce(ifce)%x(2)
    if(r_centre<r_inner.or.r_centre>r_outer)return


    ! profile 1 (volume; inj on/off)
    ref_ts=0.
    ref_tf=inj_dur
    ref_dt=0.025e-3*t_sc
    ref_pow=1.
    pr(1)=pr_inj(time, &
      ref_ts,ref_ts+ref_dt,ref_pow, &
      ref_tf-ref_dt,ref_tf,one/ref_pow)


    ! profile 2 (velocity; inj on/off)
    ref_tf=ref_tf+0.*ref_dt
    pr(2)=pr_inj(time, &
      ref_ts,ref_ts+ref_dt,ref_pow, &
      ref_tf-ref_dt,ref_tf,one/ref_pow)


    ! profile 3 (radial variation; volume)
    grad=one/r_outer
    pr(3)=(r_centre*grad)**(0.3)


    ! profile 4 (radial variation; velocity)
    grad=one/r_outer
    pr(4)=(r_centre*grad)**(0.7)


    ! moments
    if(mom(1)<small)then
      do i=1,size(mom_pow)
        mom(i)=gamma_moment(mom_pow(i), &
          mom_pow,vof_inj,r32_inj,skw_inj)
      end do
      print*,'mu:',mom(lb:ub)
    end if

    ! reference moments
    grd(igrd)%iph(2)%ieq(2)%phi_max(lb:ub)=mom(lb:ub)


    ! velocity
    vel(1)=spd_inj*cos(pr(4)*ang_inj)
    vel(2)=spd_inj*sin(pr(4)*ang_inj)


    ! moments and velocities
    do i=lb,ub
      phs_p(2)%eqn(2)%phi(i,1,1)=pr(3)*mom(i)*pr(1)
      phs_p(2)%eqn(1)%phi(i,:,1)=vel*pr(2) !! *pr(2)
    end do

    ! continuum velocity
    phs_p(1)%eqn(1)%phi(1,:,1)= &
      phs_p(2)%eqn(1)%phi(i,:,1)



    contains


    ! injection profile
    function pr_inj(t,ts,tsp,bs,tfp,tf,bf)
      real(wp),intent(in)::t,ts,tsp,bs,tfp,tf,bf
      real(wp)::pr_inj,grad,const

      pr_inj=zero
      if(t>ts.and.t<tf)then
        if(t<tsp)then
          grad=tsp-ts
          pr_inj=(t-ts)/(tsp-ts)
          pr_inj=pr_inj**(bs)
        else if(t>tfp)then
          grad=tfp-tf
          pr_inj=one-(t-tfp)/(tf-tfp)
          pr_inj=pr_inj**(bf)
        else
          pr_inj=one
        end if
      end if
    end function
  end subroutine
end module
