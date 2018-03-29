module cfd_define_m
  use grd_math_m
  use cfd_constr_m
  implicit none
  contains



  ! set coefficients: frc, dns, vis
  subroutine set_cof(igrd,ifce,ivol,iph,ieq)
    integer,intent(in)::igrd,ifce,ivol,iph,ieq
    integer::isc,ipr
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p
    type(eqn_cof_t),dimension(:),pointer::cof_p
    real(wp),pointer::vis,tvis
    real(wp)::pr

    grd_p=>grd(igrd)
    iph_p=>grd_p%iph(iph)
    ieq_p=>iph_p%ieq(ieq)

    vol_p=>grd_p%vol(ivol)
    phs_p=>vol_p%phs(iph)
    if(ifce>0)then
      fce_p=>grd_p%fce(ifce)
      phs_p=>fce_p%phs(iph)
    end if

    eqn_p=>phs_p%eqn(ieq)
    cof_p=>eqn_p%cof



    if(iph==1)then

      ! volume fraction
      if(.not.associated(eqn_p%cof(1)%frc))then
        do isc=1,ieq_p%n_sc
          eqn_p%cof(isc)%frc=>phs_p%prp(1)%rho
        end do
      end if

      ! density
      if(.not.associated(eqn_p%cof(1)%dns))then
        do isc=1,ieq_p%n_sc
          eqn_p%cof(isc)%dns=>phs_p%prp(2)%rho
        end do
      end if

      ! viscosity
      if(.not.associated(eqn_p%cof(1)%vis))then
        do isc=1,ieq_p%n_sc
          allocate(eqn_p%cof(isc)%vis(1))
        end do
      end if
      vis=>phs_p%prp(3)%rho(1)
      do isc=1,ieq_p%n_sc
        eqn_p%cof(isc)%vis=vis
      end do

      if(turb)then
        tvis=>phs_p%prp(4)%rho(1)
        if(ieq==1)then
          do isc=1,ieq_p%n_sc
            eqn_p%cof(isc)%vis= &
              eqn_p%cof(isc)%vis+tvis
          end do
        end if
        if(ieq==2)then
          eqn_p%cof(1)%vis=max(vis,tvis)
          eqn_p%cof(2)%vis=max(vis,tvis*0.77)
        end if
      end if
    end if



    if(iph==2)then

      ! volume fraction
      if(.not.associated(eqn_p%cof(1)%frc))then
        do isc=1,ieq_p%n_sc
          eqn_p%cof(isc)%frc=>phs_p%prp(1)%rho
        end do
      end if

      ! density
      if(.not.associated(eqn_p%cof(1)%dns))then
        do isc=1,ieq_p%n_sc
          if(ieq==1)then
            eqn_p%cof(isc)%dns=>phs_p%eqn(2)%phi(isc,1,:)
          else
            eqn_p%cof(isc)%dns=>phs_p%prp(1)%rho !( = one )
          end if
        end do
      end if

      ! viscosity
      if(.not.associated(eqn_p%cof(1)%vis))then
        do isc=1,ieq_p%n_sc
          eqn_p%cof(isc)%vis=>phs_p%prp(2)%rho !( = zero )
        end do
      end if
    end if
  end subroutine



  ! set boundary conditions
  subroutine set_bnd(igrd,ifce,ivol,iph,ieq,isc,ddx,pole,rhs,set,neu)
    use spr_inject_m
    integer,intent(in)::igrd,ifce,ivol,iph,ieq,isc
    real(wp),dimension(:,:),intent(in)::ddx
    real(wp),intent(inout)::pole
    real(wp),dimension(:),intent(inout)::rhs
    logical,intent(inout)::set
    logical,intent(inout)::neu
    integer::ist,i3
    real(wp)::c
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p,phs_fce
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p,eqn_fce
    type(eqn_cof_t),dimension(:),pointer::cof_fce
    real(wp),pointer::dns,vis,tvis,def,ke,ed
    real(wp),dimension(:),pointer::vel
    real(wp),dimension(:,:),pointer::phi_p,phi_fce
    real(wp)::ch_spd,ch_len,tlev
    real(wp)::n_plus,u_plus
    integer::bnd_ty

    grd_p=>grd(igrd)
    vol_p=>grd_p%vol(ivol)
    phs_p=>vol_p%phs(iph)
    eqn_p=>phs_p%eqn(ieq)
    phi_p=>eqn_p%phi(:,:,1)

    fce_p=>grd_p%fce(ifce)
    phs_fce=>fce_p%phs(iph)
    eqn_fce=>phs_fce%eqn(ieq)
    cof_fce=>eqn_fce%cof
    phi_fce=>eqn_fce%phi(:,:,1)

    iph_p=>grd_p%iph(iph)
    ieq_p=>iph_p%ieq(ieq)

    bnd_ty=fce_p%bnd_ty-mod(fce_p%bnd_ty,10)



    if(iph==1)then
      if(ieq==1)then
        if(fce_p%bnd_ty==11)then
!           phi_fce(:,1)=50.*l_sc/t_sc
        end if
      end if


      if(ieq==2)then
        if(fce_p%bnd_ty==11)then
          if(isc==1)then
            ke=>phi_fce(1,1)
            ed=>phi_fce(2,1)

            vel=>phs_p%eqn(1)%phi(1,:,1)
            ch_spd=vec_mag(vel)
            ch_len=0.005*l_sc
            tlev=0.05

            ke=phi_p(1,1)
            ke=max(ke,1.5*(ch_spd*tlev)**2)
            ke=1.5*(ch_spd*tlev)**2
            ed=2.35*ke**1.5/ch_len
          end if
        end if

        if(bnd_ty==30)then
          neu=.true.

          if(isc==2)then
            dns=>phs_fce%prp(2)%rho(1)
            vis=>phs_fce%prp(3)%rho(1)
            tvis=>phs_fce%prp(4)%rho(1)
            ke=>phi_p(1,1)
            ed=>phi_p(2,1)

            ch_spd=0.548*sqrt(max(ke,small))
            ch_len=fce_p%l_pn
            n_plus=ch_spd*ch_len*dns/vis

            tvis=0.
            if(n_plus>11.93)then
              u_plus=2.44*log(9.80*n_plus)
              tvis=ch_spd*dns*ch_len/u_plus
            end if

            ed=ch_spd**3/(0.41*ch_len)
            set=.true.
          end if
        end if
      end if
    end if



    if(iph==2)then
      if(ieq==1.and.isc==1)then
        if(fce_p%bnd_ty==11)then
          call spr_inject(igrd,ifce,ivol,iph,ieq)
        end if
      end if
    end if
  end subroutine



  ! calculate source terms
  subroutine calc_rhs(igrd,ivol,iph,ieq,isc,ddx,pole,rhs,set,neu)
    use spr_hydro_m
    integer,intent(in)::igrd,ivol,iph,ieq,isc
    real(wp),dimension(:,:),intent(in)::ddx
    real(wp),intent(inout)::pole
    real(wp),dimension(:),intent(inout)::rhs
    logical,intent(inout)::set
    logical,dimension(:),intent(inout)::neu
    integer::i,j,ist
    real(wp)::tvis0,tvis1,tvisr,urf,pk,ek
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p
    type(eqn_cof_t),dimension(:),pointer::cof_p
    real(wp),pointer::dns,vis,tvis,def,ke,ed
    real(wp),dimension(:,:),pointer::phi_p


    grd_p=>grd(igrd)
    vol_p=>grd_p%vol(ivol)
    iph_p=>grd_p%iph(iph)
    phs_p=>vol_p%phs(iph)
    ieq_p=>iph_p%ieq(ieq)
    eqn_p=>phs_p%eqn(ieq)
    phi_p=>eqn_p%phi(:,:,1)
    cof_p=>eqn_p%cof



    if(iph==1)then
      if(ieq==1)then
        if(turb)then
          def=>phs_p%prp(5)%rho(1)
          def=0.
          do i=1,n_dx
            do j=1,n_dx
              def=def+0.5*ddx(i,j) &
                *(ddx(i,j)+ddx(j,i))
            end do
            def=def-0.333*ddx(i,i)
          end do
        end if
      end if

      if(ieq==2)then
        dns=>phs_p%prp(2)%rho(1)
        tvis=>phs_p%prp(4)%rho(1)
        def=>phs_p%prp(5)%rho(1)
        ke=>eqn_p%phi(1,1,1)
        ed=>eqn_p%phi(2,1,1)
        ek=max(ed,0.)/max(ke,small)
        pk=2.*tvis*def

        if(isc==1)then
          rhs=rhs+pk*vol_p%frc_vol
          pole=pole+ek*dns*vol_p%frc_vol
        end if

        if(isc==2)then
          rhs=rhs+1.44*ek*pk*vol_p%frc_vol
          pole=pole+1.92*ek*dns*vol_p%frc_vol
        end if
      end if
    end if



    if(iph==2)then
      if(ieq==1.and.isc==1)then
        do ist=1,ieq_p%n_st
          if(phs_p%eqn(2)%f_st(ist))then
            call spr_hydro(igrd,ivol,iph,ist)
          end if
        end do
      end if
    end if
  end subroutine



  ! calculate surface source terms
  subroutine calc_srhs(igrd,ifce,ip,in,iph,ieq,isc)
    use spr_conv_m
    integer,intent(in)::igrd,ifce,ip,in,iph,ieq,isc

    integer::ist
    type(grd_t),pointer::grd_p
    type(iph_t),pointer::iph_p
    type(ieq_t),pointer::ieq_p

    grd_p=>grd(igrd)
    iph_p=>grd_p%iph(iph)
    ieq_p=>iph_p%ieq(ieq)

    if(iph==2)then
      if(ieq==1.and.isc==1)then
        do ist=1,ieq_p%n_st
          call spr_conv(igrd,ifce,ip,in,iph,ist)
        end do
      end if
    end if
  end subroutine
end module
