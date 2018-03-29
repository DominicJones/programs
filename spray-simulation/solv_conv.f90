module spr_conv_m
use grd_math_m
use mom_solve_m
use grd_data_m
use spr_data_m
use cfd_constr_m
implicit none
contains
  subroutine spr_conv(igrd,ifce,ip,in,iph,ist)
    integer,intent(in)::igrd,ifce,ip,in,iph,ist

    type(grd_t),pointer::grd_p
    type(fce_t),pointer::fce_p
    type(vol_t),pointer::vol_p,vol_n
    type(eqn_t),pointer::vel_p,vel_n,vel_f
    type(eqn_t),pointer::mom_p,mom_n,mom_f

    integer::lb,ub,icp,k,ik,i3
    real(wp)::mu
    real(wp),dimension(n_dx)::vp,vn,vf,srhs


    if(.not.hyd)return
    return
    grd_p=>grd(igrd)
    fce_p=>grd_p%fce(ifce)

    vol_p=>grd_p%vol(ip)
    vel_p=>vol_p%phs(iph)%eqn(1)
    mom_p=>vol_p%phs(iph)%eqn(2)

    if(in>0)then
      vol_n=>grd_p%vol(in)
      vel_n=>vol_n%phs(iph)%eqn(1)
      mom_n=>vol_n%phs(iph)%eqn(2)
    else
      vel_f=>fce_p%phs(iph)%eqn(1)
      mom_f=>fce_p%phs(iph)%eqn(2)
    end if


    call cp_id(igrd,iph,2,ist,lb=lb,ub=ub)
    call cp_id(igrd,iph,2,ist,icp=mom_id(3),id=i3)


    k=0
    do icp=lb,ub
      k=k+1; ik=mom_pow(k)

      if(ik/=3)then
        srhs=zero

        if(in>0)then
          vp=vel_p%phi(i3,:,1)-vel_p%phi(icp,:,1)
          vn=vel_n%phi(i3,:,1)-vel_n%phi(icp,:,1)
          vf=fce_p%w_fp*vn+(one-fce_p%w_fp)*vp
          mu=fce_p%w_fp*mom_n%phi(icp,1,1) &
            +(one-fce_p%w_fp)*mom_p%phi(icp,1,1)
        else
          vf=vel_f%phi(i3,:,1)-vel_f%phi(icp,:,1)
          mu=mom_f%phi(icp,1,1)
        end if

        srhs=-mu*fce_p%s &
          *dot_product(vf,fce_p%n)*vf

        if(vol_p%phs(iph)%eqn(2)%f_st(ist))then
          vel_p%rhs(icp,:)=vel_p%rhs(icp,:)+srhs
        end if
        if(in>0)then
          if(vol_n%phs(iph)%eqn(2)%f_st(ist))then
            vel_n%rhs(icp,:)=vel_n%rhs(icp,:)-srhs
          end if
        end if
      end if
    end do
  end subroutine
end module
