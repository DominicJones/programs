module spr_hydro_m
use spr_impact_m
use grd_math_m
use mom_solve_m
use grd_data_m
use spr_data_m
use cfd_constr_m
! use f95_lapack !! requires: f95_lapack
implicit none
contains
  subroutine spr_hydro(igrd,ivol,iph,ist)
    integer,intent(in)::igrd,ivol,iph,ist
    integer::icp
    integer,parameter::n_seg=25
    real(wp),dimension(n_seg,2)::pdf
    real(wp),dimension(n_seg,n_dx)::vdf
    real(wp),dimension(n_dx)::tvel,dvel
    real(wp),dimension(n_dx,2)::cvel

    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(iph_t),dimension(:),pointer::iph_p
    type(phs_t),dimension(:),pointer::phs_p
    type(ipr_t),dimension(:),pointer::ipr_p
    type(prp_t),dimension(:),pointer::prp_p

    type(eqn_t),pointer::vel,vel_d,mom_d
    real(wp),pointer::dns,vis,ke
    real(wp),pointer::dns_d,vis_d,stn_d

    integer,parameter::n_poly=4
    real(wp)::r_norm,num,mean,vari,skew,ratio,dr
    real(wp)::vel_fl(n_dx),dx
    real(wp)::vdf_x(size(mom_pow)+1),vdf_y(n_dx,size(mom_pow)+1)
    real(wp)::mat(n_poly,size(mom_pow)),rhs(max(n_poly,size(mom_pow)),n_dx)
    real(wp)::sou(n_seg,3),sor(n_dx)
    real(wp)::spd,spd_mn,r_st(2),n_bre(2),c_bre(2)
    real(wp)::re_num,oh_num,we_num,c_drg
    real(wp)::ts,we_crit,pdf_mn
    real(wp)::f_gam,we_l,we_b,s_mn,s_mx
    real(wp)::p_col,p_coa,e_col,n_col
    real(wp),dimension(2)::rn
    integer::n_coa
    real,dimension(:),allocatable::mu_norm
    integer::i3,i2,lb,ub,lbp,mthd,i,j,k,ik,ib
    integer::ivf,ifce



    ! set-up
    if(.not.hyd)return

    grd_p=>grd(igrd)
    vol_p=>grd_p%vol(ivol)

    iph_p=>grd_p%iph
    phs_p=>vol_p%phs

    vel=>phs_p(1)%eqn(1)
    vel_d=>phs_p(2)%eqn(1)
    mom_d=>phs_p(2)%eqn(2)

    ke=>phs_p(1)%eqn(2)%phi(1,1,1)
    dns=>phs_p(1)%prp(2)%rho(1)
    vis=>phs_p(1)%prp(3)%rho(1)
    dns_d=>phs_p(2)%prp(4)%rho(1)
    vis_d=>phs_p(2)%prp(5)%rho(1)
    stn_d=>phs_p(2)%prp(6)%rho(1)


!     mom_d%phi(1,1,1)=514.707
!     mom_d%phi(2,1,1)=8.74173
!     mom_d%phi(3,1,1)=0.172381
!     mom_d%phi(4,1,1)=0.00381252
!     vel_d%phi(1,:,1)=92.3584
!     vel_d%phi(2,:,1)=93.0596
!     vel_d%phi(3,:,1)=94.4267
!     vel_d%phi(4,:,1)=96.1894
!     vel%phi(1,:,1)=23.8605
! 
!     mom_d%phi(1,1,1)=394.432!514.707
!     mom_d%phi(2,1,1)=5.88195!8.74173
!     mom_d%phi(3,1,1)=0.101923!0.172381
!     mom_d%phi(4,1,1)=0.00196866!0.00381252
!     vel_d%phi(1,:,1)=68.0059!92.3584
!     vel_d%phi(2,:,1)=69.3555!93.0596
!     vel_d%phi(3,:,1)=72.696!94.4267
!     vel_d%phi(4,:,1)=77.0512!96.1894
!     vel%phi(1,:,1)=5.25923!23.8605
! 
!     mom_d%phi(1,1,1)=23.3382!394.432!514.707
!     mom_d%phi(2,1,1)=0.169361!5.88195!8.74173
!     mom_d%phi(3,1,1)=0.003112!0.101923!0.172381
!     mom_d%phi(4,1,1)=0.000155959!0.00196866!0.00381252
!     vel_d%phi(1,:,1)=55.1771!68.0059!92.3584
!     vel_d%phi(2,:,1)=51.5851!69.3555!93.0596
!     vel_d%phi(3,:,1)=52.341!72.696!94.4267
!     vel_d%phi(4,:,1)=58.2868!77.0512!96.1894
!     vel%phi(1,:,1)=0.589269!5.25923!23.8605


!     mom_d%phi(1,1,1)=18467.9
!     mom_d%phi(2,1,1)=169.204
!     mom_d%phi(3,1,1)=2.03169
!     mom_d%phi(4,1,1)=0.0294586
!     vel_d%phi(1,:,1)=101.471
!     vel_d%phi(2,:,1)=102.005
!     vel_d%phi(3,:,1)=102.318
!     vel_d%phi(4,:,1)=102.437
!     vel%phi(1,:,1)=37.4231
! 
! 
!     mom_d%phi(1,1,1)=6201.11
!     mom_d%phi(2,1,1)=48.2283
!     mom_d%phi(3,1,1)=0.482322
!     mom_d%phi(4,1,1)=0.00561636
!     vel_d%phi(1,:,1)=74.447
!     vel_d%phi(2,:,1)=68.5326
!     vel_d%phi(3,:,1)=65.2815
!     vel_d%phi(4,:,1)=67.5938
!     vel%phi(1,:,1)=11.1468
! 
! 
!     mom_d%phi(1,1,1)=684.43
!     mom_d%phi(2,1,1)=0.894604
!     mom_d%phi(3,1,1)=0.00125258
!     mom_d%phi(4,1,1)=2.70287E-05
!     vel_d%phi(1,:,1)=66.7343
!     vel_d%phi(2,:,1)=57.5367
!     vel_d%phi(3,:,1)=50.7085
!     vel_d%phi(4,:,1)=58.2257
!     vel%phi(1,:,1)=1.2158


    call cp_id(igrd,iph,2,ist,lb=lb,ub=ub)
    call cp_id(igrd,iph,2,ist,icp=mom_id(mpr),id=i3)
    call cp_id(igrd,iph,2,ist,icp=mom_id(mpr-1),id=i2)
    r_norm=mom_d%phi(i3,1,1)/mom_d%phi(i2,1,1)
!     print*,'rn',r_norm



    ! instaneous velocity
    call random_number(vel_fl)
    vel_fl=(vel_fl-0.5)*sqrt(2./3.*ke)
    vel_fl=vel%phi(1,:,1)+vel_fl
    spd_mn=1.*l_sc/t_sc



    ! PDF
    num=zero; lbp=0
    mthd=pdf_mth

    if(mthd/=5)then
      allocate(mu_norm(size(mom_pow)))
      num=mom_d%phi(lb,1,1)
      mu_norm=mom_d%phi(lb:ub,1,1)
      mu_norm=mu_norm/mu_norm(1)

      mean=mu_norm(2)/mu_norm(1)
      vari=mu_norm(3)-mean**2
      skew=mu_norm(4)-3*mean*mu_norm(3)+2*mean**3
      ratio=mean/vari/(vari/skew)

      if(.not.(vari>small.and.mean/sqrt(vari)>1.1 &
        .and.abs(ratio)<2.0))then
        mthd=5; lbp=mp1
      end if
    end if

    call solve_pdf(mom_d%phi(lb+lbp:ub,1,1), &
      mom_pow(1+lbp:),num,r_norm,pdf,mth_opt=mthd)

    dr=pdf(2,1)-pdf(1,1)



    ! VDF
    if(vdf_mth==1)then
      cvel(:,2)=vel_d%phi(i2,:,1)-vel_fl
      cvel(:,2)=cvel(:,2)*mom_d%phi(i2,1,1) &
        /(mom_d%phi(i2,1,1)**(1.-vel_exp) &
        *mom_d%phi(i3,1,1)**vel_exp)

      do i=1,size(pdf,1)
        vdf(i,:)=vel_fl &
          +cvel(:,2)*(pdf(i,1)+small)**vel_exp
      end do


    else if(vdf_mth==2)then
      k=0
      do j=lb,ub
        k=k+1
        rhs(k,:)=mom_d%phi(j,1,1)*vel_d%phi(j,:,1)
      end do
      rhs=rhs/mom_d%phi(lb,1,1)

      do j=1,size(mom_pow)
        do i=1,n_poly
          k=mom_pow(j)+(i-1)
          mat(i,j)=int_pdf(pdf,k)
        end do
      end do

!       call la_gelsd(mat,rhs) !! requires: f95_lapack

      do i=1,size(pdf,1)
        vdf(i,:)=zero
        do j=1,n_poly
          vdf(i,:)=vdf(i,:)+rhs(j,:)*pdf(i,1)**mom_pow(j)
        end do
      end do


    else if(vdf_mth==3)then
      dx=pdf(n_seg,1)/size(mom_pow)
      vdf_x=zero
      do i=2,size(vdf_x)
        vdf_x(i)=vdf_x(i-1)+dx
      end do

      vdf_y=zero
      vdf_y(:,1)=cvel(:,1)
      do i=2,size(vdf_x)
        vdf_y(:,i)=vel_d%phi(lb+i-2,:,1)
      end do

      j=1
      do i=1,size(pdf,1)
        cvel(:,1)=(vdf_y(:,j+1)-vdf_y(:,j))/dx
        cvel(:,2)=vdf_y(:,j)-cvel(:,1)*vdf_x(j)

        vdf(i,:)=cvel(:,1)*pdf(i,1)+cvel(:,2)

        if(pdf(i,1)>vdf_x(j+1))j=j+1
        j=min(j,size(vdf_x)-1)
      end do
    end if

!     open(101,file='distr.dat')
!     do i=1,size(pdf,1)
!       write(101,*),pdf(i,:),vdf(i,1)
!       write(*,*),pdf(i,:),vdf(i,1)
!     end do
!     close(101)
!     sou=zero; sor=zero



    ! Source Terms
    do i=2,size(pdf,1)
      if(pdf(i,1)<0.01*r_norm)cycle


      if(col)then
        do j=2,i-1
          if(pdf(j,1)<0.01*r_norm)cycle

          call random_number(rn(1:1))
          if(rn(1)<0.5)cycle

          e_col=spd*(2.*pdf(j,1))**2 &
            *dns_d/(9.*vis*2.*pdf(i,1))
          if(e_col>1.125)then
            e_col=(1.+0.75*log(2.*e_col) &
              /(e_col-1.214))**(-2)
          else
            e_col=0.
          end if

          spd=vec_mag(vdf(i,:)-vdf(j,:))
          n_col=pi*(pdf(i,1)+pdf(j,1))**2 &
            *spd*num*pdf(j,2)*dr*dt*e_col
          p_col=exp(-n_col)
          call random_number(rn)

          n_coa=0
          if(rn(1)>p_col)then
            f_gam=pdf(i,1)/pdf(j,1)
            f_gam=f_gam**3-2.4*f_gam**2+2.7*f_gam
            we_l=dns_d*spd**2*pdf(j,1)/stn_d
            we_b=2.4*f_gam*rn(2)**3
            p_coa=min(1.,2.4*f_gam/we_l)

            if(p_coa>rn(2).and.we_l>we_b)then
              s_mn=p_col; s_mx=p_col
              n_coa=1
              do
                s_mx=s_mx &
                  +p_col*n_col**(n_coa)/fact(n_coa)
                if(rn(1)>s_mn.and.rn(1)<s_mx)exit
                s_mn=s_mx
                if(n_coa>=25)exit
                n_coa=n_coa+1
              end do
            end if
          end if

          if(n_coa>0)then
            k=0
            do icp=lb,ub
              k=k+1; ik=mom_pow(k)
              sor(1)=((pdf(i,1)**3+n_coa*pdf(j,1)**3)**(ik/3.) &
                -pdf(i,1)**ik-n_coa*pdf(j,1)**ik) &
                *num*pdf(i,2)*dr/dt*vol_p%v

              mom_d%rhs(icp,1)=mom_d%rhs(icp,1)+sor(1)
!               if(ik==0)sou(i,2)=sou(i,2)+sor(1)
            end do
          end if
        end do
      end if


      dvel=vdf(i,:)-vel_fl
      spd=vec_mag(dvel)


      if(bre.and.spd>spd_mn)then
        r_st=zero; n_bre=zero; c_bre=zero


        if(bre_mod==1)then ! R&D
          we_num=pdf(i,1)*dns*spd**2/stn_d
          we_crit=8.4 ! 3.6 - 8.4
          if(we_num>we_crit)then
            r_st(1)=we_crit/(dns*spd**2/stn_d)
            c_bre(1)=pi/4.*(2.*pdf(i,1))**(1.5) &
              *sqrt(dns_d)/sqrt(stn_d)
          end if

          re_num=2.*pdf(i,1)*dns*spd/vis
          if(we_num/sqrt(re_num)>0.5)then
            r_st(2)=(0.5/(sqrt(vis)*sqrt(dns)*spd**(1.5)))**2
            r_st(2)=0.5*r_st(2)
            c_bre(2)=20. ! 2 - 20
            c_bre(2)=c_bre(2)*pdf(i,1)/spd*sqrt(dns_d/dns)
          end if


        else if(bre_mod==2)then ! P&E
          oh_num=vis_d/sqrt(2.*pdf(i,1)*dns_d*stn_d)
          we_num=2.*pdf(i,1)*dns*spd**2/stn_d
          we_crit=12.*(1.+1.077*oh_num**(1.6))

          if(we_num>we_crit)then
            if(we_num>12..and.we_num<18.)then
              ts=6.*(we_num-12.)**(-0.25)
            else if(we_num>18..and.we_num<45.)then
              ts=2.45*(we_num-12.)**(0.25)
            else if(we_num>45..and.we_num<350.)then
              ts=14.1*(we_num-12.)**(-0.25)
            else if(we_num>350..and.we_num<2670.)then
              ts=0.766*(we_num-12.)**(0.25)
            else if(we_num>2760.)then
              ts=5.5
            end if

            r_st(1)=sqrt(dns/dns_d)*(0.375*ts+0.2274*ts**2)
            r_st(1)=we_crit*(one-r_st(1))**(-2)*stn_d/(dns*spd**2)
            r_st(1)=0.5*r_st(1)
            c_bre(1)=ts*2.*pdf(i,1)/spd*sqrt(dns_d/dns)
          end if


        else if(bre_mod==3)then ! H&F
          we_num=pdf(i,1)*dns*spd**2/stn_d
          if(we_num>6..and.we_num<1000.)then
            oh_num=vis_d/sqrt(2.*pdf(i,1)*dns_d*stn_d)
            r_st(1)=6.2*pdf(i,1)*sqrt(sqrt(dns_d/dns)) &
              *sqrt(vis_d/(2.*pdf(i,1)*dns_d*spd))
            c_bre(1)=5.*2.*pdf(i,1)*sqrt(dns_d/dns) &
              /(spd*(one-oh_num/7.))
          end if
        end if


        do j=1,2
          if(r_st(j)>0.001*r_norm)then
            n_bre(j)=max(one,(pdf(i,1)/r_st(j))**(3))
            c_bre(j)=num*pdf(i,2)*dr/(c_bre(j)+small)
          end if
        end do
      end if


      if(drg)then
        c_drg=zero
        re_num=2.*pdf(i,1)*dns_d*spd/vis
        if(re_num<50.)re_num=50.
        if(re_num>1000.)re_num=1000.
        c_drg=24./re_num+3.48/re_num**(0.313)
        c_drg=3./8.*dns/dns_d &
          *c_drg*spd/(pdf(i,1)+small)
        c_drg=c_drg*num*pdf(i,2)*dr
      end if



      k=0
      do icp=lb,ub
        k=k+1; ik=mom_pow(k)


        if(bre)then
          do ib=1,size(n_bre)
            if(n_bre(ib)>one)then
              sor(1)=(n_bre(ib)**((3-ik)/3.)-one) &
                *c_bre(ib)*pdf(i,1)**(ik)*vol_p%v
              mom_d%rhs(icp,1)=mom_d%rhs(icp,1)+sor(1)
!               if(ik==0)sou(i,1)=sor(1)
            end if
          end do
        end if


        if(drg)then
          sor=-c_drg*pdf(i,1)**(ik)*dvel*vol_p%v
          vel_d%rhs(icp,:)=vel_d%rhs(icp,:)+sor
          if(ik==3)vel%rhs(1,:)=vel%rhs(1,:)-sor*dns_d
!           if(ik==3)sou(i,3)=sor(1)
        end if
      end do
    end do

!     open(101,file='source.dat')
!     do i=1,size(pdf,1)
!       write(101,*)pdf(i,:),sum(sou(:i,1)),sum(sou(:i,2)),sum(sou(:i,3))
!       write(*,*)pdf(i,:),sum(sou(:i,1)),sum(sou(:i,2)),sum(sou(:i,3))
!     end do
!     close(101)
!     stop



    ! wall impaction
    if(imp.and.ist==inj_st.and.vol_p%bnd_v)then
      do ivf=1,vol_p%n_fce
        ifce=vol_p%fce_lst(ivf)
        if(grd_p%fce(ifce)%bnd_ty==32)then
          call spr_impact(igrd,ifce,ivol,iph,num,pdf,vdf)
        end if
      end do
    end if
  end subroutine
end module
