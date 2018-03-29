module spr_impact_m
use grd_math_m
use mom_solve_m
use grd_data_m
use spr_data_m
use cfd_constr_m
implicit none
contains
  subroutine spr_impact(igrd,ifce,ivol,iph,num,pdf,vdf)
    integer,intent(in)::igrd,ifce,ivol,iph
    real(wp),intent(in)::num
    real(wp),dimension(:,:),intent(in)::pdf
    real(wp),dimension(:,:),intent(in)::vdf
    real(wp),dimension(n_dx)::vel_in,vel_out,vel_n,vel_t,vel_spl,vel_reb

    real(wp)::theta,cof_n,cof_t,mu_in,mu_out
    real(wp)::c_inj,c_spl,c_reb,c,dr,we_0,grad_we
    real(wp)::la_cof,la_num,we_num,num_int

    integer::i,k,ik,icp
    integer::lb_inj,ub_inj,iinj
    integer::lb_reb,ub_reb,ireb
    integer::lb_spl,ub_spl,ispl
    logical::wet_wall

    type(grd_t),pointer::grd_p
    type(fce_t),pointer::fce_p
    type(vol_t),pointer::vol_p
    type(iph_t),dimension(:),pointer::iph_p
    type(phs_t),dimension(:),pointer::phs_p
    type(ipr_t),dimension(:),pointer::ipr_p
    type(prp_t),dimension(:),pointer::prp_p

    type(eqn_t),pointer::vel,vel_d,vel_d_v,mom_d,mom_d_v,liq_flm
    real(wp),pointer::dns,vis,dns_d,vis_d,stn_d



    ! set-up
    grd_p=>grd(igrd)
    fce_p=>grd_p%fce(ifce)
    vol_p=>grd_p%vol(ivol)
    iph_p=>grd_p%iph
    phs_p=>fce_p%phs

    vel_d=>phs_p(2)%eqn(1)
    mom_d=>phs_p(2)%eqn(2)
    vel_d_v=>vol_p%phs(2)%eqn(1)
    mom_d_v=>vol_p%phs(2)%eqn(2)
    liq_flm=>vol_p%phs(2)%eqn(3)

    dns=>phs_p(1)%prp(2)%rho(1)
    vis=>phs_p(1)%prp(3)%rho(1)
    dns_d=>phs_p(2)%prp(4)%rho(1)
    vis_d=>phs_p(2)%prp(5)%rho(1)
    stn_d=>phs_p(2)%prp(6)%rho(1)

    call cp_id(igrd,iph,2,inj_st,lb=lb_inj,ub=ub_inj)
    call cp_id(igrd,iph,2,reb_st,lb=lb_reb,ub=ub_reb)
    call cp_id(igrd,iph,2,spl_st,lb=lb_spl,ub=ub_spl)



    ! zero boundary values
    vel_d%phi(lb_spl:ub_spl,:,1)=0.
    mom_d%phi(lb_spl:ub_spl,1,1)=0.
    vel_d%phi(lb_reb:ub_reb,:,1)=0.
    mom_d%phi(lb_reb:ub_reb,1,1)=0.

    vel_spl=0.
    vel_reb=0.



    ! wet wall condition
    wet_wall=.false.
    grad_we=(wall_we(1)-wall_we(2))/vof_wall
    we_0=wall_we(1)-grad_we*vof_wall
    la_cof=grad_we*liq_flm%phi(1,1,1)+we_0

    if(la_cof<wall_we(1))wet_wall=.true.
    la_cof=max(wall_we(1),la_cof)



    ! boundary conditions
    dr=pdf(2,1)-pdf(1,1)
    do i=1,size(pdf,1)
      if(pdf(i,1)<small)cycle

      vel_in=vdf(i,:)
      vel_n=fce_p%n*dot_product(vel_in,fce_p%n)
      vel_t=vel_in-vel_n

      theta=acos(dot_product(vel_in,vel_t) &
        /(vec_mag(vel_in)*vec_mag(vel_t)+small))

      we_num=2.*pdf(i,1)*dns_d &
        *vec_mag(vel_n)**2/stn_d
      la_num=2.*pdf(i,1)*dns_d*stn_d/vis_d**2


      k=0
      do icp=lb_inj,ub_inj
        k=k+1; ik=mom_pow(k)
        iinj=lb_inj-1+k
        ispl=lb_spl-1+k
        ireb=lb_reb-1+k
        num_int=num*pdf(i,2)*pdf(i,1)**(ik)*dr


        ! splashing spray
        if(we_num>la_cof*la_num**(-0.18))then
          cof_n=0.25
          cof_t=0.25
          vel_out=-cof_n*vel_n+cof_t*vel_t

          mom_d%phi(ispl,1,1)=mom_d%phi(ispl,1,1) &
            +n_spl**((3-ik)/3.)*num_int

          if(ik==3)then
            vel_spl=vel_spl+vel_out*num_int &
              /(mom_d_v%phi(iinj,1,1)+small)
          end if
        end if


        ! rebounding spray
        if(wet_wall.and.we_num<5.0)then
          cof_n=0.993-1.76*theta+1.56*theta**2-0.49*theta**3
          cof_t=0.7143
          vel_out=-cof_n*vel_n+cof_t*vel_t

          mom_d%phi(ireb,1,1)=mom_d%phi(ireb,1,1) &
            +num_int

          if(ik==3)then
            vel_reb=vel_reb+vel_out*num_int &
              /(mom_d_v%phi(iinj,1,1)+small)
          end if
        end if
      end do
    end do



    k=0
    do icp=lb_inj,ub_inj
      k=k+1; ik=mom_pow(k)
      iinj=lb_inj-1+k
      ireb=lb_reb-1+k
      ispl=lb_spl-1+k


      ! injected spray passing through the wall
      vel_d%phi(iinj,:,1)=vel_d_v%phi(iinj,:,1)
      mom_d%phi(iinj,:,1)=mom_d_v%phi(iinj,:,1)


      ! splashing and rebounding spray velocities
      vel_d%phi(ispl,:,1)=vel_spl
      vel_d%phi(ireb,:,1)=vel_reb


      ! liquid film volume fraction source term
      if(ik==3)then
        c_inj=mom_d%phi(iinj,1,1) &
          *dot_product(vel_d%phi(iinj,:,1),fce_p%n)

        c_spl=mom_d%phi(ispl,1,1) &
          *dot_product(vel_d%phi(ispl,:,1),fce_p%n)

        c_reb=mom_d%phi(ireb,1,1) &
          *dot_product(vel_d%phi(ireb,:,1),fce_p%n)

        c=fce_p%s*c_sph*(c_inj+c_spl+c_reb)
        liq_flm%rhs=liq_flm%rhs+c
      end if
    end do
  end subroutine
end module
