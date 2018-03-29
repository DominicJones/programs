module cfd_solve_m
  use grd_geom_m
  use cfd_setup_m
  use cfd_outp_m
  use spr_data_m
  use mom_math_m
  implicit none
  contains


  ! solve flow
  subroutine cfd_solve(igrd,outp_dir)
    integer,intent(in)::igrd
    character(len=*),intent(in)::outp_dir

    integer::n_ph,iph,n_sc,isc,n_cp,icp,n_st,ist
    integer::n_ts,its,ict,iprt,n_iter,iter
    real(wp)::frc,dns,vis,tvis,stn,ch_spd,urfp
    real(wp)::ch_len,tlev,urft,tm(2,1),reduc
    logical::vof,prt
    character::yn



    ! SETTINGS
    m_sc=ten**(3)
    t_sc=ten**(3)

    n_ph=2
    dbg=.false.
    turb=.true.
    spr=.true.

    iprt=50
    n_ts=250; dt=2.e-6*t_sc
    n_iter=15; reduc=1.e-3
    call test_grad(igrd)
    n_grad=2

    call c_phs(igrd,n_ph)



    ! CONTINUUM
    iph=1
    frc=1.; dns=5.5; vis=1.846e-5
    dns=dns*m_sc/l_sc**3
    vis=vis*m_sc/(l_sc*t_sc)

    urfp=0.1 ! not required

    if(turb)then
      ch_spd=40.*l_sc/t_sc
      ch_len=0.01*l_sc
      tlev=0.02
      tm(1,1)=1.5*(ch_spd*tlev)**2
      tm(2,1)=2.35*tm(1,1)**(1.5)/ch_len
      tvis=0.09*dns*tm(1,1)**2/tm(2,1)
      urft=0.65
      print*,'ke,ed,tvis/vis:',tm(:,1),tvis/vis
    end if

    call c_prp(igrd,iph,5)
    call c_prp_rho(igrd,iph,ipr=1,rho=frc,s_tl=.true.)  ! vol. frac.
    call c_prp_rho(igrd,iph,ipr=2,rho=dns,s_tl=.false.) ! density
    call c_prp_rho(igrd,iph,ipr=3,rho=vis,s_tl=.false.) ! molec. visc
    call c_prp_rho(igrd,iph,ipr=4,rho=tvis,s_tl=.false.) ! turb. visc
    call c_prp_rho(igrd,iph,ipr=5,rho=zero,s_tl=.false.) ! deformation

    call c_eqn(igrd,iph,2)
    call c_eqn_phi(igrd,iph,ieq=1,n_sc=1,n_vc=n_dx, &
      urf=0.7_wp,c_bl=0.9_wp,c_sch=1,rhs=.true.)
    call c_eqn_phi(igrd,iph,ieq=2,n_sc=2,n_vc=1, &
      urf=urft,c_bl=0.5_wp,c_sch=1,phi=tm)



    ! SPRAY
    if(spr)then
      iph=2
      frc=0.; dns=702.; vis=5.65e-4; stn=0.0226
      dns=dns*m_sc/l_sc**3
      vis=vis*m_sc/(l_sc*t_sc)
      stn=stn*m_sc/t_sc**2

      bre_mod=2
      bre=.true.
      col=.false.
      drg=.true.
      imp=.true.

      r_orif=0.5e-3*l_sc !! 0.25mm
      sheet_th=r_orif
      inj_dur=0.5e-3*t_sc
      vof_inj=0.99
      skw_inj=7.
      r32_inj=25.e-6*l_sc
      spd_inj=120.*l_sc/t_sc
      ang_inj=20.*pi/180.

      vel_exp=0.4
      vdf_mth=1

      n_cp=3; mp1=1; mpr=3; pdf_mth=5
      allocate(mom_pow(n_cp))
      mom_pow=(/2,3,4/)
      n_st=1; n_sc=n_st*n_cp

      call c_prp(igrd,iph,n_pr=7)
      call c_prp_rho(igrd,iph,ipr=1,rho=one,s_tl=.false.)  ! 1
      call c_prp_rho(igrd,iph,ipr=2,rho=zero,s_tl=.false.) ! 0
      call c_prp_rho(igrd,iph,ipr=3,rho=frc,s_tl=.false.) ! vol. frac.
      call c_prp_rho(igrd,iph,ipr=4,rho=dns,s_tl=.false.) ! density
      call c_prp_rho(igrd,iph,ipr=5,rho=vis,s_tl=.false.) ! viscosty
      call c_prp_rho(igrd,iph,ipr=6,rho=stn,s_tl=.false.) ! surf. tens.
      call c_prp_rho(igrd,iph,ipr=7,rho=zero,s_tl=.false.) ! SMR

      call c_eqn(igrd,iph,n_eq=3)
      call c_eqn_phi(igrd,iph,ieq=1,n_sc=n_sc,n_vc=n_dx,n_cp=n_cp, &
        urf=0.7_wp,c_bl=0.0_wp,c_sch=1,rhs=.true.)
      call c_eqn_phi(igrd,iph,ieq=2,n_sc=n_sc,n_vc=1,n_cp=n_cp, &
        urf=0.7_wp,c_bl=0.9_wp,c_sch=4,rhs=.true.,lb_cut=ten**(-7))
      call c_eqn_phi(igrd,iph,ieq=3,n_sc=1,n_vc=1, &
        urf=0.9_wp,rhs=.true.)
    end if



    ! ITERATE
    ict=0
    do its=1,n_ts
      time=its*dt

      do iter=1,n_iter
        mx_res_n=zero
        c_in=zero; c_out=zero


        if(n_tl>1.and.iter==1)then
          call stor_prop(igrd,iph=1,ipr=1)
        end if


        if(spr)then
          hyd=.true.
          if(iter<=3)hyd=.false.
          iph=2
          isc=0
          do ist=1,n_st
            do icp=1,n_cp
              isc=isc+1
              vof=.false.
              if(mom_pow(icp)==3)vof=.true.
              call cfd_setup_eqn(igrd,iph,1,isc,iter,vof)
              call cfd_setup_eqn(igrd,iph,2,isc,iter,.false.)
            end do
          end do
          if(imp)then
            call zero_spd(igrd)
            call cfd_setup_eqn(igrd,iph,3,1,iter,.true.)
          end if
          call crop_spr(igrd,iter)
        end if


        iph=1
        call cfd_setup_eqn(igrd,iph,1,1,iter,.true.)
        call cfd_setup_pp(igrd,iph,1,iter,n_pp=1,urf=urfp)
        if(turb)then
          call cfd_setup_eqn(igrd,iph,2,1,iter,.false.)
          call cfd_setup_eqn(igrd,iph,2,2,iter,.false.)
          call mod_tvis(igrd,urft)
        end if


        write(*,'(1x,a6,2i6,es14.4)') &
          'REDUC:',its,iter,mx_res_n
        if(mx_res_n<reduc)exit
      end do



      ! OUTPUT
      prt=.false.
      if(n_tl==1)prt=.true.
      if(n_tl>1.and.mod(its,iprt)==0)then
        ict=ict+1
        prt=.true.
      end if

      if(prt)then
        iph=1
        call eqn_outp(igrd,iph,1,ict,trim(outp_dir))
        if(turb)then
!           call eqn_outp(igrd,iph,2,ict,trim(outp_dir))
!           call prp_outp(igrd,iph,4,ict,trim(outp_dir))
        end if

        if(spr)then
          iph=2
          call eqn_outp(igrd,iph,1,ict,trim(outp_dir))
          call eqn_outp(igrd,iph,2,ict,trim(outp_dir))
!           call eqn_outp(igrd,iph,3,ict,trim(outp_dir))
          call prp_outp(igrd,iph,3,ict,trim(outp_dir))
          call prp_outp(igrd,iph,7,ict,trim(outp_dir))
        end if
      end if

      if(n_tl==1)exit
    end do
  end subroutine



  ! store properties
  subroutine stor_prop(igrd,iph,ipr)
    integer,intent(in)::igrd,iph,ipr
    integer::ivol,itl

    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p

    grd_p=>grd(igrd)
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      prp_p=>vol_p%phs(iph)%prp(ipr)
      do itl=n_tl-1,1,-1
        prp_p%rho(itl+1)=prp_p%rho(itl)
      end do
    end do
  end subroutine



  ! modify turbulent viscosity
  subroutine mod_tvis(igrd,urf)
    integer,intent(in)::igrd
    real(wp),intent(in)::urf
    integer::ivol
    real(wp)::tvis1
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(phs_t),pointer::phs_p
    real(wp),pointer::dns,tvis,ke,ed

    grd_p=>grd(igrd)
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      phs_p=>vol_p%phs(1)

      dns=>phs_p%prp(2)%rho(1)
      tvis=>phs_p%prp(4)%rho(1)
      ke=>phs_p%eqn(2)%phi(1,1,1)
      ed=>phs_p%eqn(2)%phi(2,1,1)

      tvis1=max(ke,zero)**2 &
        /max(ed,small)*0.09*dns
      tvis=urf*tvis1+(one-urf)*tvis
    end do
  end subroutine



  ! crop spray and calc vof & r32
  subroutine crop_spr(igrd,iter)
    integer,intent(in)::igrd,iter
    integer::ivol,ifce,iph,ieq,isc,ipr,itl,ist
    integer::lb,ub,icp,icpm,i3,i2
    real(wp)::sfrc,frc,r32,phi_min,mom2,mom3

    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(fce_t),pointer::fce_p
    type(iph_t),pointer::iph_p
    type(phs_t),pointer::phs_p
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p


    grd_p=>grd(igrd)
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)

      iph=2
      ieq=2
      iph_p=>grd_p%iph(iph)
      ieq_p=>iph_p%ieq(ieq)
      phs_p=>vol_p%phs(iph)
      eqn_p=>phs_p%eqn(ieq)


      r32=zero
      frc=zero
      do ist=1,ieq_p%n_st
        call cp_id(igrd,iph,ieq,ist=ist,lb=lb,ub=ub)
        eqn_p%f_st(ist)=.true.
        icpm=0
        do icp=lb,ub
          icpm=icpm+1
          phi_min=ieq_p%lb_cut*ieq_p%phi_max(icpm)
          if(eqn_p%phi(icp,1,1)<phi_min)then
            eqn_p%f_st(ist)=.false.
            exit
          end if
        end do

        if(eqn_p%f_st(ist))then
          call cp_id(igrd,iph,ieq,ist,icp=mom_id(3),id=i3)
          sfrc=c_sph*eqn_p%phi(i3,1,1)
          frc=frc+sfrc

          if(mom_id(2)>0)then
            call cp_id(igrd,iph,ieq,ist,icp=mom_id(2),id=i2)
            r32=r32+eqn_p%phi(i3,1,1) &
              /eqn_p%phi(i2,1,1)*sfrc
          else
            mom2=gamma_moment(2,mom_pow,mu=eqn_p%phi(lb:ub,1,1))
            mom3=gamma_moment(3,mom_pow,mu=eqn_p%phi(lb:ub,1,1))
            r32=r32+mom3/(mom2+small)*sfrc
          end if
        else
          phs_p%eqn(1)%phi(lb:ub,:,1)=zero
          phs_p%eqn(2)%phi(lb:ub,1,1)=zero
        end if
      end do

      phs_p%prp(3)%rho(1)=frc
      phs_p%prp(7)%rho(1)=r32/(frc+small)

      frc=frc+phs_p%eqn(3)%phi(1,1,1)
      vol_p%phs(1)%prp(1)%rho(1)=one-frc
    end do


    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)

      if(fce_p%in==0)then
        frc=zero
        do ist=1,ieq_p%n_st
          call cp_id(igrd,2,2,ist=ist,icp=mom_id(3),id=i3)
          frc=frc+fce_p%phs(2)%eqn(2)%phi(i3,1,1)*c_sph
        end do
        fce_p%phs(1)%prp(1)%rho(1)=one-frc
      end if
    end do
  end subroutine
end module
