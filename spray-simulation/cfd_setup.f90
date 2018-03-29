module cfd_setup_m
  use solv_bicg_m
  use solv_cg_m
  use cfd_define_m
  implicit none
  contains



  ! zero speed
  subroutine zero_spd(igrd)
    integer,intent(in)::igrd
    integer::ifce
    do ifce=1,size(ng)
      grd(igrd)%fce(ifce)%spd=zero
    end do
  end subroutine



  ! construct equation
  subroutine cfd_setup_eqn(igrd,iph,ieq,isc,it,vof)
    integer,intent(in)::igrd,iph,ieq,isc,it
    logical,intent(in)::vof

    integer::ivol,ifce,ip,in,ib,ic,id,ivc,itl,n_vc,bnd_ty
    integer::n_bicg
    real(wp)::urf,urf_r,reduc,dvdt,mx_res0_n,frc_vis
    real(wp)::c,d,frc,dns,vis,c_bl,t_bl,dt_r
    real(wp),dimension(n_tl)::frc_dns
    real(wp), &
      dimension(grd(igrd)%iph(iph)%ieq(ieq)%n_vc):: &
      c_exp,c_imp,d_exp,d_imp,cor,t,res0_n,srhs

    real(wp),dimension(size(res0_n))::phi_fce
    real(wp),dimension(:),pointer::phi_p_aux,phi_n_aux,phi_bnd

    real(wp),dimension(:,:,:),allocatable::ddx
    real(wp),dimension(size(res0_n))::phi_u
    real(wp),dimension(:),pointer::phi_c,phi_d
    type(vol_t),pointer::vol_c
    real(wp),dimension(n_dx)::dx_pn
    real(wp)::phi_cn,phi_fn,co_num,theta,gamma

    real(wp),dimension(size(pl),size(res0_n))::rhs
    logical,dimension(size(pl))::set
    logical,dimension(size(ng))::neu

    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p,eqn_n,eqn_fce
    type(eqn_cof_t),pointer::cof_p,cof_n,cof_fce
    type(phi_t),dimension(grd(igrd)%n_vbf)::phi
    real(wp)::tmp,phi_mx,rhs_mx


if(dbg)print*,'initialise'

    n_bicg=5
    grd_p=>grd(igrd)
    ieq_p=>grd_p%iph(iph)%ieq(ieq)

    ieq_p%c_in(isc)=zero
    ieq_p%c_out(isc)=zero
    if(.not.ieq_p%solv)return

    urf=ieq_p%urf
    urf_r=one/urf
    c_bl=ieq_p%c_bl
    t_bl=ieq_p%t_bl
    if(n_tl==2)t_bl=zero
    dt_r=one/dt
    n_vc=size(res0_n)

    allocate(phi_p_aux(n_vc),phi_n_aux(n_vc))
    allocate(ddx(size(pl),n_vc,n_dx))

    pl=zero; rhs=zero
    set=.false.; neu=.false.

    do ivol=1,size(pl)
      vol_p=>grd_p%vol(ivol)
      vol_p%cs=zero
      if(n_tl>1.and.n_vc>1)then
        vol_p%co_num=zero
      end if
    end do


if(dbg)print*,'point phi()%vc / coefficients'


    tmp=zero; phi_mx=zero
    do ivol=1,size(pl)
      vol_p=>grd_p%vol(ivol)
      eqn_p=>vol_p%phs(iph)%eqn(ieq)
      phi(ivol)%vc=>eqn_p%phi(isc,:,1)
      if(isc==1)then
        call set_cof(igrd,0,ivol,iph,ieq)
      end if
      tmp=max(tmp,maxval(abs(phi(ivol)%vc)))
      phi_mx=max(phi_mx,tmp)
    end do

    do ifce=1,size(ng)
      if(grd_p%fce(ifce)%in==0)then
        fce_p=>grd_p%fce(ifce)
        ip=fce_p%ip; ib=fce_p%ib
        eqn_fce=>fce_p%phs(iph)%eqn(ieq)
        phi(ib)%vc=>eqn_fce%phi(isc,:,1)
        if(isc==1)then
          call set_cof(igrd,ifce,ip,iph,ieq)
        end if
      end if
    end do


if(dbg)print*,'gradient / Courant num.'


    if(phi_mx>small)then
      call solv_grad(igrd,phi,ddx,n_grad,neu=0)
    else
      ddx=zero
    end if


if(dbg)print*,'boundary conditions'


    do ifce=1,size(ng)
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip
      in=fce_p%in

      if(in>0)then
        call calc_srhs(igrd,ifce,ip,in,iph,ieq,isc)
      else
        call calc_srhs(igrd,ifce,ip,0,iph,ieq,isc)

        ib=fce_p%ib
        call set_bnd(igrd,ifce,ip,iph,ieq,isc, &
          ddx(ip,:,:),pl(ip),rhs(ip,:),set(ip),neu(ifce))
        tmp=max(tmp,maxval(abs(phi(ib)%vc)))
        phi_mx=max(phi_mx,tmp)
      end if
    end do


if(dbg) &
print*,'source terms'


    do ivol=1,size(pl)
      vol_p=>grd_p%vol(ivol)
      eqn_p=>vol_p%phs(iph)%eqn(ieq)
      cof_p=>eqn_p%cof(isc)
      vol_p%frc_vol=cof_p%frc(1)*vol_p%v

      call calc_rhs(igrd,ivol,iph,ieq,isc, &
        ddx(ivol,:,:),pl(ivol),rhs(ivol,:),set(ivol),neu)

      if(ieq_p%rhs)then
        rhs(ivol,:)=rhs(ivol,:)+eqn_p%rhs(isc,:)
        eqn_p%rhs(isc,:)=zero
      end if
    end do
    rhs_mx=maxval(abs(rhs))
    if(phi_mx<small.and.rhs_mx<small)return


if(dbg) &
print*,'face integration'


    do ifce=1,size(ng)
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip
      vol_p=>grd_p%vol(ip)
      eqn_p=>vol_p%phs(iph)%eqn(ieq)
      cof_p=>eqn_p%cof(isc)

      do ivc=1,n_vc
        phi_p_aux(ivc)=phi(ip)%vc(ivc) &
          +dot_product(ddx(ip,ivc,:),fce_p%dx_p)
      end do

      if(fce_p%in>0)then
        in=fce_p%in
        vol_n=>grd_p%vol(in)
        eqn_n=>vol_n%phs(iph)%eqn(ieq)
        cof_n=>eqn_n%cof(isc)

        do ivc=1,n_vc
          phi_n_aux(ivc)=phi(in)%vc(ivc) &
            +dot_product(ddx(in,ivc,:),fce_p%dx_n)
        end do

        phi_fce=fce_p%w_fp*phi_n_aux &
          +(one-fce_p%w_fp)*phi_p_aux

        if(n_vc>1)then
          fce_p%spd=dot_product(phi_fce,fce_p%n)
        end if

        frc=fce_p%w_fp*cof_n%frc(1) &
          +(one-fce_p%w_fp)*cof_p%frc(1)
        dns=fce_p%w_fp*cof_n%dns(1) &
          +(one-fce_p%w_fp)*cof_p%dns(1)
        vis=fce_p%w_fp*cof_n%vis(1) &
          +(one-fce_p%w_fp)*cof_p%vis(1)

        fce_p%frc_dns=frc*dns
        frc_vis=frc*vis

        c=zero; d=zero; cor=zero

        if(fce_p%frc_dns>small.and.abs(fce_p%spd)>small)then
          c=fce_p%frc_dns*fce_p%s*fce_p%spd
          vol_p%cs=vol_p%cs+c
          vol_n%cs=vol_n%cs-c

          ! UDS
          if(c>zero)c_imp=phi(ip)%vc
          if(c<zero)c_imp=phi(in)%vc

          ! CDS
          if(ieq_p%c_sch==0)then
            c_exp=phi_fce

          else if(ieq_p%c_sch>0)then
            if(c>zero)then
              ic=ip; vol_c=>vol_p
              phi_c=>phi_p_aux
              phi_d=>phi_n_aux
              dx_pn=(vol_n%x+fce_p%dx_n)-(vol_p%x+fce_p%dx_p)
            else
              ic=in; vol_c=>vol_n
              phi_c=>phi_n_aux
              phi_d=>phi_p_aux
              dx_pn=(vol_p%x+fce_p%dx_p)-(vol_n%x+fce_p%dx_n)
            end if

            do ivc=1,n_vc
              phi_u(ivc)=phi_d(ivc)-2._wp &
                *dot_product(ddx(ic,ivc,:),dx_pn)
            end do


            ! TVD
            if(ieq_p%c_sch<=3)then
              c_exp=(phi_c-phi_u) &
                /(phi_d-phi_c+small)

              ! min-mod
              if(ieq_p%c_sch==1)then
                c_exp=max(zero,min(one,c_exp))

              ! van leer
              else if(ieq_p%c_sch==2)then
                c_exp=(c_exp+abs(c_exp)) &
                  /(one+abs(c_exp)+small)

              ! superbee
              else if(ieq_p%c_sch==3)then
                c_exp=max(zero,min(one,2._wp*c_exp), &
                  min(2._wp,c_exp))
              end if

              c_exp=phi_c+0.5_wp*c_exp &
                *(phi_d-phi_c)


            ! HRIC
            else if(ieq_p%c_sch>3)then
              do ivc=1,n_vc
                phi_cn=(phi_c(ivc)-phi_u(ivc)) &
                  /(phi_d(ivc)-phi_u(ivc)+small)

                phi_fn=phi_cn
                if(phi_cn>0._wp.and.phi_cn<0.5_wp)then
                  phi_fn=2._wp*phi_cn
                else if(phi_cn>0.5_wp.and.phi_cn<1._wp)then
                  phi_fn=1._wp
                end if

                co_num=vol_c%co_num
                if(co_num>0.3_wp.and.co_num<0.7_wp)then
                  phi_fn=phi_cn+(phi_fn-phi_cn) &
                    *(0.7_wp-co_num)/(0.7_wp-0.3_wp)
                else if(co_num>0.7_wp)then
                  phi_fn=phi_cn
                end if

                theta=dot_product(ddx(ic,ivc,:),fce_p%n) &
                  /(vec_mag(ddx(ic,ivc,:))+small)
                theta=sqrt(abs(theta))

                phi_fn=phi_fn*theta+phi_cn*(one-theta)

                gamma=(phi_d(ivc)-phi_u(ivc)) &
                  /(phi_d(ivc)-phi_c(ivc)+small)
                gamma=gamma*(one-phi_fn)
                gamma=max(min(gamma,one),zero)

                c_exp(ivc)=gamma*phi_c(ivc) &
                  +(one-gamma)*phi_d(ivc)
              end do
            end if
          end if

          cor=c_bl*c*(c_imp-c_exp)
        end if

        if(frc_vis>small)then
          d=frc_vis*fce_p%s/fce_p%l_pn
          d_imp=phi(in)%vc-phi(ip)%vc
          d_exp=phi_n_aux-phi_p_aux
          cor=cor+d*(d_exp-d_imp)
        end if

        ng(ifce)%p=-d+min(c,zero)
        ng(ifce)%n=-d-max(c,zero)
        if(set(ip))ng(ifce)%p=zero
        if(set(in))ng(ifce)%n=zero

        pl(ip)=pl(ip)-ng(ifce)%p
        pl(in)=pl(in)-ng(ifce)%n
        rhs(ip,:)=rhs(ip,:)+cor
        rhs(in,:)=rhs(in,:)-cor

!         srhs=zero
!         call calc_srhs(igrd,ifce,ip,in,iph,ieq,isc, &
!           frc,dns,srhs)
!         rhs(ip,:)=rhs(ip,:)+srhs
!         rhs(in,:)=rhs(in,:)-srhs

      else
        ib=fce_p%ib
        bnd_ty=fce_p%bnd_ty-mod(fce_p%bnd_ty,10)
        eqn_fce=>fce_p%phs(iph)%eqn(ieq)
        cof_fce=>eqn_fce%cof(isc)

        if(n_vc>1)then
          fce_p%spd=dot_product(phi(ib)%vc,fce_p%n)

          if(bnd_ty==20.and. &
            dot_product(phi(ip)%vc,fce_p%n)<-small)then
            t=phi(ip)%vc
            t=t-fce_p%n*dot_product(t,fce_p%n)
            t=t/(vec_mag(t)+small)
            phi(ip)%vc=t*dot_product(phi(ip)%vc,t)
            fce_p%spd=dot_product(phi(ip)%vc,fce_p%n)
          end if
        end if

        frc=cof_fce%frc(1)
        dns=cof_fce%dns(1)
        vis=cof_fce%vis(1)

        fce_p%frc_dns=frc*dns
        frc_vis=frc*vis

        c=zero; d=zero; cor=zero

        if(fce_p%frc_dns>small.and.abs(fce_p%spd)>small)then
          c=fce_p%frc_dns*fce_p%s*fce_p%spd
          vol_p%cs=vol_p%cs+c

          if(bnd_ty==10.and.vof)then
            if(iph==1)c_in=c_in-frc*fce_p%s*fce_p%spd
            if(iph==2)c_in=c_in-c_sph*dns*fce_p%s*fce_p%spd
          end if
          if(bnd_ty==20.and.vof)then
            if(iph==1)c_out=c_out+frc*fce_p%s*fce_p%spd
            if(iph==2)c_out=c_out+c_sph*dns*fce_p%s*fce_p%spd
          end if
        end if

        if(frc_vis>small)then
          d=frc_vis*fce_p%s/fce_p%l_pn
          d_imp=phi(ib)%vc-phi(ip)%vc
          d_exp=phi(ib)%vc-phi_p_aux
        end if

        ng(ifce)%p=-d+min(c,zero)
        if(set(ip))ng(ifce)%p=zero

        if(bnd_ty==40)then
          if(n_vc==1)then
            neu(ifce)=.true.
          else if(n_vc>1)then
            phi(ib)%vc=phi(ip)%vc-fce_p%n &
              *dot_product(phi(ip)%vc,fce_p%n)
            if(d>small)then
              d_exp=-2._wp*fce_p%n &
                *dot_product(phi(ip)%vc,fce_p%n)
            end if
          end if
        end if

        if(bnd_ty==30.and.n_vc>1.and.d>small)then
          t=phi(ip)%vc+small
          if(vec_mag(phi(ib)%vc)>small)t=phi(ib)%vc
          t=t-fce_p%n*(dot_product(t,fce_p%n))
          t=t/(vec_mag(t)+small)
          phi(ib)%vc=t*(dot_product(phi(ib)%vc,t))
          d_exp=phi(ib)%vc-t*dot_product(phi(ip)%vc,t)
        end if

        if(bnd_ty==20.or.neu(ifce))then
          phi(ib)%vc=phi(ip)%vc
          ng(ifce)%p=zero
        end if

        if(d>small)then
          cor=d*(d_exp-d_imp)
        end if

        pl(ip)=pl(ip)-ng(ifce)%p
        rhs(ip,:)=rhs(ip,:)+cor
        rhs(ip,:)=rhs(ip,:)-ng(ifce)%p*phi(ib)%vc

!         srhs=zero
!         call calc_srhs(igrd,ifce,ip,0,iph,ieq,isc, &
!           frc,dns,srhs)
!         rhs(ip,:)=rhs(ip,:)+srhs
      end if
    end do


if(dbg)print*,'volume integration / source terms'


    do ivol=1,size(pl)
      vol_p=>grd_p%vol(ivol)
      eqn_p=>vol_p%phs(iph)%eqn(ieq)
      cof_p=>eqn_p%cof(isc)
      vol_p%frc_vol=cof_p%frc(1)*vol_p%v


      ! Spray Equations
      if(iph==2.and.n_vc==1.and.n_tl>1)then
        pl(ivol)=pl(ivol)+vol_p%cs
      end if
      if(iph==2.and.n_vc>1.and.pl(ivol)<small)then
        rhs(ivol,:)=zero
      end if
      if(iph==2.and.cof_p%frc(1)*cof_p%dns(1)<small)then
        rhs(ivol,:)=zero
      end if


      ! Continuum Equations
      if(iph==1.and.n_vc>1)then
        rhs(ivol,:)=rhs(ivol,:) &
          -vol_p%dpdx*vol_p%frc_vol
        if(n_dx==2.and.cyl)then
          rhs(ivol,2)=rhs(ivol,2) &
            -2._wp*cof_p%vis(1)/vol_p%x(2)**2 &
            *phi(ivol)%vc(2)*vol_p%frc_vol
        end if
      end if

!       if(ieq_p%rhs)then
!         rhs(ivol,:)=rhs(ivol,:)+eqn_p%rhs(isc,:)
!         eqn_p%rhs(isc,:)=zero
!       end if


      if(n_tl>1)then
        if(it==1)then
          do itl=n_tl-1,1,-1
            eqn_p%phi(isc,:,itl+1)= &
              eqn_p%phi(isc,:,itl)
          end do
        end if

        if(size(cof_p%frc)==1)then
          frc_dns=cof_p%frc(1)
        else
          frc_dns=cof_p%frc
        end if
        if(size(cof_p%dns)==1)then
          frc_dns=frc_dns*cof_p%dns(1)
        else
          frc_dns=frc_dns*cof_p%dns
        end if

        frc_dns(1)=frc_dns(1)*(one+0.5_wp*t_bl)
        frc_dns(2)=-frc_dns(2)*(one+t_bl)
        if(n_tl>2)frc_dns(3)=frc_dns(3)*0.5_wp*t_bl
        frc_dns=frc_dns*vol_p%v*dt_r

        if(iph==1.or.(iph==2.and.n_vc==1))then
          pl(ivol)=pl(ivol)+frc_dns(1)
          if(n_vc>1)vol_p%cs=vol_p%cs+sum(frc_dns)
        else
          pl(ivol)=pl(ivol)-sum(frc_dns(2:n_tl))
        end if
        do ivc=1,n_vc
          rhs(ivol,ivc)=rhs(ivol,ivc) &
            -sum(frc_dns(2:n_tl) &
            *eqn_p%phi(isc,ivc,2:n_tl))
        end do
      end if


      if(set(ivol))then
        rhs(ivol,:)=phi(ivol)%vc
        pl(ivol)=one
      end if

      pl(ivol)=urf_r*pl(ivol)
      rhs(ivol,:)=rhs(ivol,:) &
        +(one-urf)*pl(ivol)*phi(ivol)%vc
    end do

    deallocate(phi_p_aux,phi_n_aux)
    deallocate(ddx)


if(dbg)print*,'residual / solve'


    call solv_bicg(2,0.1_wp,n_bicg,phi,rhs,res0_n,.false.)

    if(it==1)ieq_p%res_n(isc)=zero
    mx_res0_n=maxval(res0_n)
    if(mx_res0_n>small)then
      ieq_p%res_n(isc)=max(ieq_p%res_n(isc),mx_res0_n)
      mx_res0_n=mx_res0_n/ieq_p%res_n(isc)
      mx_res_n=max(mx_res_n,mx_res0_n)
    end if


! if(dbg) &
print*,'reduc:',iph,ieq,isc,mx_res0_n


    if(abs(mx_res0_n)>large)stop
  end subroutine



  ! construct pressure correction equation
  subroutine cfd_setup_pp(igrd,iph,ieq,it,n_pp,urf)
    integer,intent(in)::igrd,iph,ieq,it,n_pp
    real(wp),intent(in)::urf

    real(wp),dimension(size(pl))::pl_r
    real(wp),dimension(size(pl),1,n_dx)::ddx
    real(wp),dimension(size(pl),1)::rhs,rhs_p
    real(wp),dimension(1)::res0_n

    integer::ivol,ifce,ip,in,ib,ipp,ivf,bnd_ty,n_cg
    real(wp)::mx_res0_n,m_sum,delta,cs,dc
    real(wp)::c,d,frc,dns,p_r,l_pn,v,dp,c_rat,t_bl
    real(wp)::p_aux,n_aux,phi_p_aux,phi_n_aux,phi_fce
    real(wp),save::pp_res_n=0.
    real(wp),dimension(n_dx)::ddx_fce,x_p_aux

    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p,eqn_n,eqn_fce
    type(eqn_cof_t),pointer::cof_p,cof_n,cof_fce
    type(phi_t),dimension(grd(igrd)%n_vbf)::phi
    logical::simpler


if(dbg)print*,'initialise'


    grd_p=>grd(igrd)

    n_cg=20
    simpler=.true.
!     simpler=.false.
    pl_r=one/(pl+small); pl=zero

    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      rhs(ivol,1)=-vol_p%cs
    end do

    delta=ten**(-10)
    c_rat=c_out/(c_in+delta)
    if(c_rat<delta)c_rat=one
    c_rat=one/c_rat
    if(abs(c_in)<delta)c_rat=zero


if(dbg) &
print*,'c in/out/rat:',c_in,c_out,c_rat,sum(rhs)


if(dbg)print*,'face integration'


    do ifce=1,size(ng)
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip
      vol_p=>grd_p%vol(ip)

      if(fce_p%in>0)then
        in=fce_p%in
        vol_n=>grd_p%vol(in)

        phi_p_aux=vol_p%p(1) &
          +dot_product(vol_p%dpdx,fce_p%dx_p)
        phi_n_aux=vol_n%p(1) &
          +dot_product(vol_n%dpdx,fce_p%dx_n)

        dp=phi_n_aux-phi_p_aux
        ddx_fce=fce_p%w_fp*vol_n%dpdx &
          +(one-fce_p%w_fp)*vol_p%dpdx
        dp=dp-dot_product(ddx_fce,vol_n%x-vol_p%x)

        p_r=fce_p%w_fp*pl_r(in) &
          +(one-fce_p%w_fp)*pl_r(ip)
        v=fce_p%w_fp*vol_n%v &
          +(one-fce_p%w_fp)*vol_p%v
        d=fce_p%frc_dns*p_r*fce_p%s*v/fce_p%l_pn
        c=-d*dp

        ng(ifce)%p=-d
        ng(ifce)%n=-d

        pl(ip)=pl(ip)-ng(ifce)%p
        pl(in)=pl(in)-ng(ifce)%n

        rhs(ip,1)=rhs(ip,1)-c
        rhs(in,1)=rhs(in,1)+c
      else
        bnd_ty=fce_p%bnd_ty-mod(fce_p%bnd_ty,10)
        if(bnd_ty==20)then
          c=fce_p%frc_dns*fce_p%s*fce_p%spd
          c=c*(c_rat-one)
          rhs(ip,1)=rhs(ip,1)-c
        end if
      end if
    end do


if(dbg)print*,'solve pp'


    do ivol=1,size(pl)
      allocate(phi(ivol)%vc(1))
    end do

    do ipp=1,n_pp
      do ivol=1,grd_p%n_vol
        vol_p=>grd_p%vol(ivol)
        phi(ivol)%vc=zero
      end do

      call solv_cg(1,0.1_wp,n_cg,phi,rhs,res0_n,.true.)

      mx_res0_n=res0_n(1)+small
      if(it==1)pp_res_n=zero
      pp_res_n=max(pp_res_n,mx_res0_n)
      mx_res0_n=mx_res0_n/pp_res_n
      mx_res_n=max(mx_res_n,mx_res0_n)


! if(dbg) &
print*,'redpp:',iph,ieq,ipp,mx_res0_n,sum(rhs)


      if(.not.simpler)then
if(dbg)print*,'update pressure'

        do ivol=1,grd_p%n_vol
          vol_p=>grd_p%vol(ivol)
          vol_p%p(1)=vol_p%p(1) &
            +urf*(phi(ivol)%vc(1)-phi(1)%vc(1))
        end do
      end if


if(dbg)print*,'solve grad(pp)'


      call solv_grad(igrd,phi,ddx,n_grad,neu=2)


if(dbg)print*,'update velocity'


      do ivol=1,grd_p%n_vol
        vol_p=>grd_p%vol(ivol)
        eqn_p=>vol_p%phs(iph)%eqn(ieq)
        eqn_p%phi(1,:,1)=eqn_p%phi(1,:,1) &
          -ddx(ivol,1,:)*pl_r(ivol)*vol_p%frc_vol
      end do


if(dbg)print*,'correct mass flux'


      if(simpler)rhs_p=rhs
      rhs=zero
      do ifce=1,grd_p%n_fce
        fce_p=>grd_p%fce(ifce)
        if(fce_p%in>0)then
          ip=fce_p%ip
          in=fce_p%in

          p_aux=dot_product(ddx(ip,1,:),fce_p%dx_p)
          phi_p_aux=phi(ip)%vc(1)+p_aux

          n_aux=dot_product(ddx(in,1,:),fce_p%dx_n)
          phi_n_aux=phi(in)%vc(1)+n_aux

          dp=ng(ifce)%p*(n_aux-p_aux)
          rhs(ip,1)=rhs(ip,1)-dp
          rhs(in,1)=rhs(in,1)+dp

          dp=ng(ifce)%p*(phi_n_aux-phi_p_aux) &
            /(fce_p%frc_dns*fce_p%s)
          fce_p%spd=fce_p%spd+dp
        end if
      end do
      if(simpler)rhs_p=rhs_p+rhs
    end do


if(dbg)print*,'solve grad(p)'


    do ivol=1,size(pl)
      deallocate(phi(ivol)%vc)
      phi(ivol)%vc=>grd_p%vol(ivol)%p
    end do

    if(simpler)then
      call solv_cg(1,0.1_wp,n_cg,phi,rhs_p,res0_n,.true.)

      mx_res0_n=res0_n(1)+small
      if(it==1)pp_res_n=zero
      pp_res_n=max(pp_res_n,mx_res0_n)
      mx_res0_n=mx_res0_n/pp_res_n
      mx_res_n=max(mx_res_n,mx_res0_n)

! if(dbg) &
print*,'red_p:',iph,ieq,0,mx_res0_n,sum(rhs_p)
    end if

    call solv_grad(igrd,phi,ddx,n_grad,neu=2)

    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      grd_p%vol(ivol)%dpdx=ddx(ivol,1,:)
    end do
  end subroutine



  ! solve gradient
  subroutine solv_grad(igrd,phi,grad,n_iter,neu)
    type(phi_t),dimension(:),intent(inout)::phi
    real(wp),dimension(:,:,:),intent(inout)::grad
    integer,intent(in)::igrd,n_iter,neu
    real(wp),dimension(:,:,:),allocatable::grad0
    integer::ivol,ip,in,ib,ifce,iter,ivc,n_vc
    real(wp),dimension(size(grad,2))::phi_p,phi_n,phi_f
    type(grd_t),pointer::grd_p
    type(fce_t),pointer::fce_p
    type(vol_t),pointer::vol_p,vol_n
    real(wp)::co_num

    n_vc=size(grad,2)
    if(n_iter>1)then
      allocate(grad0(size(grad,1),size(grad,2),size(grad,3)))
      grad0=zero
    end if


    grd_p=>grd(igrd)
    do iter=1,n_iter
      grad=zero

      do ifce=1,grd_p%n_fce
        fce_p=>grd_p%fce(ifce)
        ip=fce_p%ip
        vol_p=>grd_p%vol(ip)

        if(fce_p%in>0)then
          in=fce_p%in
          vol_n=>grd_p%vol(in)

          if(iter>1)then
            do ivc=1,n_vc
              phi_p(ivc)=phi(ip)%vc(ivc) &
                +dot_product(grad0(ip,ivc,:),fce_p%dx_p)
              phi_n(ivc)=phi(in)%vc(ivc) &
                +dot_product(grad0(in,ivc,:),fce_p%dx_n)

              phi_f(ivc)=fce_p%w_fp*phi_n(ivc) &
                +(one-fce_p%w_fp)*phi_p(ivc)
            end do
          else
            phi_f=fce_p%w_fp*phi(in)%vc &
              +(one-fce_p%w_fp)*phi(ip)%vc
          end if

          do ivc=1,n_vc
            grad(ip,ivc,:)=grad(ip,ivc,:) &
              +phi_f(ivc)*fce_p%s*fce_p%n*vol_p%v_r
            grad(in,ivc,:)=grad(in,ivc,:) &
              -phi_f(ivc)*fce_p%s*fce_p%n*vol_n%v_r
          end do

          if(iter==n_iter.and.n_tl>1.and.n_vc>1)then
            co_num=abs(fce_p%s*dt &
              *dot_product(phi_f,fce_p%n))
            vol_p%co_num=max(vol_p%co_num,co_num*vol_p%v_r)
            vol_n%co_num=max(vol_n%co_num,co_num*vol_n%v_r)
          end if
        else
          ib=fce_p%ib
          if(neu>0)then
            phi_f=phi(ip)%vc
            if(iter>1.and.neu==2)then
              do ivc=1,n_vc
                phi_f(ivc)=phi_f(ivc) &
                  +dot_product(grad0(ip,ivc,:),fce_p%x-vol_p%x)
              end do
            end if
          else
            phi_f=phi(ib)%vc
          end if

          do ivc=1,n_vc
            grad(ip,ivc,:)=grad(ip,ivc,:) &
              +phi_f(ivc)*fce_p%s*fce_p%n*vol_p%v_r
          end do

          if(iter==n_iter.and.n_tl>1.and.n_vc>1)then
            co_num=abs(fce_p%s*dt &
              *dot_product(phi_f,fce_p%n))
            vol_p%co_num=max(vol_p%co_num,co_num*vol_p%v_r)
          end if
        end if
      end do

      if(n_dx==2.and.cyl)then
        do ivol=1,grd_p%n_vol
          vol_p=>grd_p%vol(ivol)
          do ivc=1,n_vc
            grad(ivol,ivc,2)=grad(ivol,ivc,2) &
              -phi(ivol)%vc(ivc)*vol_p%s_f*vol_p%v_r
          end do
        end do
      end if

      if(iter/=n_iter)then
        grad0=grad
      end if
    end do
  end subroutine



  ! test convergence of solv_grad
  subroutine test_grad(igrd)
    integer,intent(in)::igrd
    type(phi_t),dimension(grd(igrd)%n_vbf)::phi
    real(wp),dimension(grd(igrd)%n_vol,1,n_dx)::grad
    integer::ivol,ifce,ip,ib
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p
    real(wp),dimension(n_dx)::err

    grd_p=>grd(igrd)
    do ivol=1,size(pl)
      vol_p=>grd_p%vol(ivol)
      allocate(phi(ivol)%vc(1))
      phi(ivol)%vc=product(vol_p%x)
    end do

    do ifce=1,size(ng)
      fce_p=>grd_p%fce(ifce)
      if(fce_p%in==0)then
        ib=fce_p%ib
        allocate(phi(ib)%vc(1))
        phi(ib)%vc=product(fce_p%x)
      end if
    end do

    do n_grad=1,6
      call solv_grad(igrd,phi,grad,n_grad,neu=0)

      err=zero
      do ivol=1,size(pl)
        vol_p=>grd_p%vol(ivol)
        err(1)=err(1)+(grad(ivol,1,1)-vol_p%x(2))**2
        err(2)=err(2)+(grad(ivol,1,2)-vol_p%x(1))**2
      end do
      err=sqrt(err/grd_p%n_vol)
      print*,'grad err',n_grad,err
    end do

    do ivol=1,size(pl)
      deallocate(phi(ivol)%vc)
    end do

    do ifce=1,size(ng)
      fce_p=>grd_p%fce(ifce)
      if(fce_p%in==0)then
        ib=fce_p%ib
        deallocate(phi(ib)%vc)
      end if
    end do
  end subroutine
end module
