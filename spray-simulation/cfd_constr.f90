module cfd_constr_m
  use grd_data_m
  implicit none
  contains



  ! construct standard data
  subroutine c_std(igrd)
    integer,intent(in)::igrd
    integer::ivol
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p

    grd_p=>grd(igrd)
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      allocate(vol_p%p(1)); vol_p%p=0.
      allocate(vol_p%dpdx(n_dx)); vol_p%dpdx=0.
    end do
  end subroutine



  ! construct phases
  subroutine c_phs(igrd,n_ph)
    integer,intent(in)::igrd,n_ph
    integer::ivol,ifce
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p

    grd_p=>grd(igrd)
    grd_p%n_ph=n_ph
    allocate(grd_p%iph(n_ph))

    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      allocate(vol_p%phs(n_ph))
    end do

    do ifce=1,grd(igrd)%n_fce
      fce_p=>grd_p%fce(ifce)
      if(fce_p%in==0)then
        allocate(fce_p%phs(n_ph))
      end if
    end do
  end subroutine



  ! construct properties
  subroutine c_prp(igrd,iph,n_pr)
    integer,intent(in)::igrd,iph,n_pr
    integer::ivol,ifce
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p

    grd_p=>grd(igrd)
    iph_p=>grd_p%iph(iph)
    iph_p%n_pr=n_pr
    allocate(iph_p%ipr(n_pr))

    do ivol=1,grd_p%n_vol
      phs_p=>grd(igrd)%vol(ivol)%phs(iph)
      allocate(phs_p%prp(n_pr))
    end do

    do ifce=1,grd_p%n_fce
      fce_p=>grd(igrd)%fce(ifce)
      if(fce_p%in==0)then
        phs_p=>fce_p%phs(iph)
        allocate(phs_p%prp(n_pr))
      end if
    end do
  end subroutine



  ! construct property storage
  subroutine c_prp_rho(igrd,iph,ipr,rho,s_tl)
    integer,intent(in)::igrd,iph,ipr
    real(kind=wp),intent(in)::rho
    logical,intent(in)::s_tl
    integer::ivol,ifce
    character(len=1)::ipr_ch
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p

    grd_p=>grd(igrd)
    iph_p=>grd_p%iph(iph)
    ipr_p=>iph_p%ipr(ipr)

    ipr_p%s_tl=s_tl
    if(s_tl)ipr_p%n_tl=n_tl
    allocate(ipr_p%rho(1))
    ipr_p%rho=rho

    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      phs_p=>vol_p%phs(iph)
      prp_p=>phs_p%prp(ipr)
      allocate(prp_p%rho(ipr_p%n_tl))
      prp_p%rho=rho
    end do

    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)
      if(fce_p%in==0)then
        phs_p=>fce_p%phs(iph)
        prp_p=>phs_p%prp(ipr)
        allocate(prp_p%rho(ipr_p%n_tl))
        prp_p%rho=rho
      end if
    end do
  end subroutine



  ! construct equations
  subroutine c_eqn(igrd,iph,n_eq)
    integer,intent(in)::igrd,iph,n_eq
    integer::ivol,ifce
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p

    grd_p=>grd(igrd)
    iph_p=>grd_p%iph(iph)
    iph_p%n_eq=n_eq
    allocate(iph_p%ieq(n_eq))

    do ivol=1,grd_p%n_vol
      phs_p=>grd(igrd)%vol(ivol)%phs(iph)
      allocate(phs_p%eqn(n_eq))
    end do

    do ifce=1,grd_p%n_fce
      fce_p=>grd(igrd)%fce(ifce)
      if(fce_p%in==0)then
        phs_p=>fce_p%phs(iph)
        allocate(phs_p%eqn(n_eq))
      end if
    end do
  end subroutine



  ! construct equation storage
  subroutine c_eqn_phi(igrd,iph,ieq,n_sc,n_vc,phi,n_cp, &
    ddx,rhs,urf,c_bl,c_sch,lb_cut)
    integer,intent(in)::igrd,iph,ieq,n_sc,n_vc
    real(wp),intent(in),optional::urf,c_bl,lb_cut
    real(kind=wp),dimension(:,:),intent(in),optional::phi
    integer,intent(in),optional::n_cp,c_sch
    logical,intent(in),optional::ddx,rhs
    integer::ivol,ifce,itl,n_st,isc
    character(len=1)::ieq_ch
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p

    grd_p=>grd(igrd)
    iph_p=>grd_p%iph(iph)
    ieq_p=>iph_p%ieq(ieq)

    ieq_p%n_sc=n_sc
    ieq_p%n_vc=n_vc
    ieq_p%n_tl=n_tl
    allocate(ieq_p%res_n(n_sc)); ieq_p%res_n=0.
    allocate(ieq_p%c_in(n_sc)); ieq_p%c_in=0.
    allocate(ieq_p%c_out(n_sc)); ieq_p%c_out=0.
    allocate(ieq_p%phi_max(n_sc)); ieq_p%phi_max=0.

    if(present(urf))ieq_p%urf=urf
    if(present(c_bl))ieq_p%c_bl=c_bl
    if(present(c_sch))ieq_p%c_sch=c_sch
    if(present(lb_cut))ieq_p%lb_cut=lb_cut
    if(present(ddx))ieq_p%ddx=ddx
    if(present(rhs))ieq_p%rhs=rhs

    if(present(n_cp))then
      allocate(ieq_p%st_id(n_sc))
      ieq_p%n_cp=n_cp
      n_st=0
      do isc=1,n_sc
        if(mod(isc,n_cp)==0)then
          n_st=n_st+1
          ieq_p%st_id(isc-n_cp+1:isc)=n_st
        end if
      end do
      ieq_p%n_st=n_st
!       do isc=1,n_sc
!         print*,'ist',iph,ieq,isc,ieq_p%st_id(isc)
!       end do
    end if

    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      phs_p=>vol_p%phs(iph)
      eqn_p=>phs_p%eqn(ieq)
      allocate(eqn_p%phi(n_sc,n_vc,n_tl))
      eqn_p%phi=0.
      if(present(phi))then
        do itl=1,n_tl
          eqn_p%phi(:,:,itl)=phi
        end do
      end if
      if(ieq_p%ddx)then
        allocate(eqn_p%ddx(n_sc,n_vc,n_dx))
        eqn_p%ddx=0.
      end if
      if(ieq_p%rhs)then
        allocate(eqn_p%rhs(n_sc,n_vc))
        eqn_p%rhs=0.
      end if
      if(present(n_cp))then
        allocate(eqn_p%f_st(n_st))
        eqn_p%f_st=.false.
      end if
      allocate(eqn_p%cof(n_sc))
    end do

    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)
      if(fce_p%in==0)then
        phs_p=>fce_p%phs(iph)
        eqn_p=>phs_p%eqn(ieq)
        allocate(eqn_p%phi(n_sc,n_vc,1))
        eqn_p%phi=0.
        allocate(eqn_p%cof(n_sc))
      end if
    end do
  end subroutine



  ! look up component index
  subroutine cp_id(igrd,iph,ieq,ist,icp,lb,ub,id)
    integer,intent(in)::igrd,iph,ieq,ist
    integer,intent(in),optional::icp
    integer,intent(out),optional::lb,ub,id
    integer::n_cp

    n_cp=grd(igrd)%iph(iph)%ieq(ieq)%n_cp
    if(present(lb))lb=(ist-1)*n_cp+1
    if(present(ub))ub=ist*n_cp
    if(present(id).and.present(icp))then
      id=(ist-1)*n_cp+icp
    end if
  end subroutine
end module
