module transport_m
implicit none
contains



!==============================================================================!
subroutine spd_conv_diff(msh,geom,pde,spd,cc,cd,a_ij,rhs)
!------------------------------------------------------------------------------!
  use vector_utils_m
  use matrix_utils_m
  use mesh_format_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::pde
  real(rp),dimension(:),intent(in)::spd
  real(rp),dimension(:),intent(in)::cc
  real(rp),dimension(:),intent(in)::cd
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),dimension(:,:),intent(inout)::rhs
!------------------------------------------------------------------------------!
  real(rp),dimension(msh%n_dx)::dx,tang,norm1,fvel,delta
  real(rp),dimension(msh%n_dx,pde%n_cmp)::fgrad,tr
  real(rp),dimension(pde%n_cmp)::c_imp,c_exp,c_cor
  real(rp),dimension(pde%n_cmp)::d_imp,d_exp,d_cor
  real(rp),dimension(pde%n_cmp)::vphi,nphi,fphi

  real(rp),dimension(msh%n_dx)::dx_cu
  real(rp),dimension(msh%n_dx,pde%n_cmp)::gradc
  real(rp),dimension(pde%n_cmp)::phiu,phic,phid,tmp,tmp2,tmp3

  real(rp)::w_fp1,w_fn1,c,c_ij,c_ji,fcc,d,d_ij,fcd,urfr
  integer::i,j,k,kii,kij,kji,kjj,ic,ib,ph_ty,ix
!------------------------------------------------------------------------------!

!$AD II-LOOP
  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)

!     kii=coo2csr(i,i,msh%elm_ja,msh%elm_ia)
    kii=msh%fce_tag(k)%lst(5)

    norm1=geom%norm(:,k)

!     dx=geom%x_pc(:,k)-geom%x_vc(:,i)
! !$AD NO-II-LOOP
!     do ic=1,pde%n_cmp
!       vphi(ic)=pde%phi(ic,i) &
!         +dot_prod(pde%grad(:,ic,i),dx)
!     end do

    if(j>0)then
!       kjj=coo2csr(j,j,msh%elm_ja,msh%elm_ia)
!       kij=coo2csr(i,j,msh%elm_ja,msh%elm_ia)
!       kji=coo2csr(j,i,msh%elm_ja,msh%elm_ia)
      kjj=msh%fce_tag(k)%lst(6)
      kij=msh%fce_tag(k)%lst(7)
      kji=msh%fce_tag(k)%lst(8)

!       dx=geom%x_nc(:,k)-geom%x_vc(:,j)
! !$AD NO-II-LOOP
!       do ic=1,pde%n_cmp
!         nphi(ic)=pde%phi(ic,j) &
!           +dot_prod(pde%grad(:,ic,j),dx)
!       end do

      w_fp1=geom%w_fp(k)
      w_fn1=1.d0-w_fp1

      delta=geom%x_vc(:,j)-geom%x_vc(:,i)
      delta=unit_vec(delta)
      delta=norm1-delta

      fgrad=w_fp1*pde%grad(:,:,j)+w_fn1*pde%grad(:,:,i)
      dx=w_fp1*geom%x_vc(:,j)+w_fn1*geom%x_vc(:,i)
      dx=geom%x_fc(:,k)-dx
      do ic=1,pde%n_cmp
        fphi(ic)=dot_prod(fgrad(:,ic),dx)
      end do
      fphi=w_fp1*pde%phi(:,j)+w_fn1*pde%phi(:,i)+fphi

      fcc=pde%lcc*(w_fp1*cc(j)+w_fn1*cc(i))
      fcd=pde%lcd*(w_fp1*cd(j)+w_fn1*cd(i))

      c=fcc*geom%area(k)*spd(k)
      d=fcd*geom%area(k)/geom%l_pn(k)

      c_ij=min(c,0.d0)
      c_ji=-max(c,0.d0)
      d_ij=-d


      if(c>0.d0)then
        c_imp=pde%phi(:,i)
      else
        c_imp=pde%phi(:,j)
      end if


      select case(pde%lc_sch)
      case(1:3)
        if(c>0.d0)then
          dx_cu=geom%x_vc(:,i)-geom%x_vc(:,j)
          gradc=pde%grad(:,:,i)
          phid=pde%phi(:,j)
          phic=pde%phi(:,i)
        else
          dx_cu=geom%x_vc(:,j)-geom%x_vc(:,i)
          gradc=pde%grad(:,:,j)
          phid=pde%phi(:,i)
          phic=pde%phi(:,j)
        end if

        do ic=1,pde%n_cmp
          phiu(ic)=phic(ic)+sum(gradc(:,ic)*dx_cu)
        end do

!         c_exp=(phic-phiu)/(phid-phic+small)
        tmp=(phic-phiu)/(phid-phic+small)

        select case(pde%lc_sch)
        case(1) ! min-mod
!           c_exp=min(one,c_exp)
!           c_exp=max(zero,c_exp)
          tmp2=min(one,tmp)
          c_exp=max(zero,tmp2)
        case(2) ! van leer
!           tmp=abs(c_exp)
!           c_exp=(c_exp+tmp)/(one+tmp+small)
          tmp2=abs(tmp)
          c_exp=(tmp+tmp2)/(one+tmp2+small)
        case(3) ! superbee
!           tmp=min(one,2*c_exp)
!           tmp2=min(2.d0,c_exp)
!           c_exp=max(tmp,tmp2)
!           c_exp=max(zero,c_exp)
          tmp2=min(one,2*tmp)
          tmp3=min(2.d0,tmp)
          c_exp=max(tmp2,tmp3)
          c_exp=max(zero,c_exp)
        end select

        c_exp=phic+half*c_exp*(phid-phic)

      case default ! central differencing
        c_exp=fphi
      end select


      c_cor=c*(c_imp-c_exp)
      c_cor=c_cor*pde%cbl


      do ic=1,pde%n_cmp
        d_cor(ic)=dot_prod(fgrad(:,ic),delta)
      end do
      d_cor=d_cor*fcd*geom%area(k)

      if(pde%lspd==1.and.pde%lcd==1)then
        do ic=1,msh%n_dx
          d_exp(ic)=dot_prod(fgrad(ic,:),norm1)
        end do
        d_exp=d_exp*fcd*geom%area(k)
        d_cor=d_cor+d_exp
      end if


!       d_imp=pde%phi(:,j)-pde%phi(:,i)
!       if(pde%lspd==1.and.pde%lcd==1)then
!         tr=transpose(fgrad)
!         d_exp=matmul(fgrad+tr,norm1)
!       else
!         d_exp=(nphi-vphi)/geom%l_pn(k)
!         do ic=1,pde%n_cmp
!           d_exp(ic)=dot_prod(fgrad(:,ic),norm1)
!         end do
!       end if
!       d_cor=fcd*geom%area(k)*d_exp-d*d_imp



      a_ij(kij)=a_ij(kij)+c_ij+d_ij
      a_ij(kji)=a_ij(kji)+c_ji+d_ij

      a_ij(kii)=a_ij(kii)-c_ij-d_ij
      a_ij(kjj)=a_ij(kjj)-c_ji-d_ij

      rhs(:,i)=rhs(:,i)+c_cor+d_cor
      rhs(:,j)=rhs(:,j)-c_cor-d_cor

    else

      ib=abs(j)
      ph_ty=msh%fce_tag(k)%lst(4)

      fgrad=pde%grad(:,:,i)
      fphi=pde%phi(:,ib)
      vphi=pde%phi(:,i) !!

!       dx=geom%x_pc(:,k)-geom%x_vc(:,i)
! !$AD NO-II-LOOP
!       do ic=1,pde%n_cmp
!         vphi(ic)=pde%phi(ic,i) &
!           +dot_prod(pde%grad(:,ic,i),dx)
!       end do

      fcc=pde%lcc*cc(ib)
      fcd=pde%lcd*cd(ib)


      ! "slip wall" momentum boundary condition on BC=39
      select case(ph_ty)
      case(wall(2))
        if(pde%lspd==1)fcd=0.d0
      end select


      c=fcc*geom%area(k)*spd(k)
      d=fcd*geom%area(k)/geom%l_pn(k)

      c_ij=min(c,0.d0)
      d_ij=-d

      d_imp=fphi-pde%phi(:,i)
!       d_exp=fphi-vphi


      if(pde%lspd==1.and.pde%lcd==1)then
        select case(ph_ty)

        ! wall boundary condition for momentum
        case(wall(1):wall(2))
          d_exp=fphi-(vphi-norm1*dot_prod(vphi,norm1))
          d_exp=d_exp/geom%l_pn(k)

        ! symmetry boundary condition for momentum
        case(symm(1):symm(2))
          d_exp=2*(0.d0-norm1*dot_prod(vphi,norm1))
          d_exp=d_exp/geom%l_pn(k)
        case default
          tr=transpose(fgrad)
          d_exp=matmul(fgrad+tr,norm1)
        end select
      else
        do ic=1,pde%n_cmp
          d_exp(ic)=sum(fgrad(:,ic)*norm1)
        end do
      end if

!       d_cor=d*(d_exp-d_imp)
      d_cor=fcd*geom%area(k)*d_exp-d*d_imp


      ! standard boundary conditions
      if(pde%lsbc==1)then
        select case(ph_ty)
        case(outlt(1):outlt(2))
          fphi=vphi
          c_ij=0; d_ij=0
        case(symm(1):symm(2))
          if(pde%lspd==0)then
            fphi=vphi
            c_ij=0; d_ij=0
          end if
          if(pde%lspd==1.and.pde%lcd==1)then
            fphi=vphi-norm1*dot_prod(vphi,norm1)
          end if
        end select
      end if


      ! distance to nearest wall eqn.
      if(pde%ldwall==1)then
        if(ph_ty<wall(1).or.ph_ty>wall(2))then
          fphi=vphi
          c_ij=0; d_ij=0
        end if
      end if


      a_ij(kii)=a_ij(kii)-c_ij-d_ij

      rhs(:,i)=rhs(:,i)-(c_ij+d_ij)*fphi
      rhs(:,i)=rhs(:,i)+d_cor

      pde%phi(:,ib)=fphi
    end if
  end do


  urfr=1.d0/pde%urf

!$AD II-LOOP
  do i=1,msh%n_elm

!     kii=coo2csr(i,i,msh%elm_ja,msh%elm_ia)
    kii=msh%elm_tag(i)%lst(5)

    vphi=pde%phi(:,i)
    a_ij(kii)=urfr*a_ij(kii)
    rhs(:,i)=rhs(:,i)+(1.d0-pde%urf)*a_ij(kii)*vphi
  end do
end subroutine



!==============================================================================!
subroutine flux_out_corr(msh,area,cc,spd,c_lim,c_rat)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use mesh_format_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  real(rp),dimension(:)::area
  real(rp),dimension(:),intent(in)::cc
  real(rp),dimension(:),intent(in)::spd
  real(rp),intent(in)::c_lim
  real(rp),intent(out)::c_rat
!------------------------------------------------------------------------------!
  real(rp)::c,c_in,c_out
  integer::i,j,k,ib,ph_ty
!------------------------------------------------------------------------------!
  c_in=0; c_out=0

!$AD II-LOOP
  do k=1,msh%n_bfce
    j=msh%fce_elm(k)%lst(2)

    if(j<0)then
      ib=abs(j)
      ph_ty=msh%fce_tag(k)%lst(4)

      select case(ph_ty)
      case(inlet(1):inlet(2))
        c=cc(ib)*area(k)*spd(k)
        c_in=c_in+c
      case(outlt(1):outlt(2))
        c=cc(ib)*area(k)*spd(k)
        c_out=c_out+c
      end select
    end if
  end do

  c_rat=c_in/(c_out+small)
  c_rat=abs(c_rat)
  c_rat=min(c_lim,c_rat)
end subroutine



!==============================================================================!
subroutine pp_setup(msh,geom,vel,pde,spd,cc,apr,a_ij,rhs)
!------------------------------------------------------------------------------!
  use vector_utils_m
  use matrix_utils_m
  use mesh_format_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::vel,pde
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),dimension(:),intent(in)::cc
  real(rp),dimension(:),intent(inout)::apr,spd
  real(rp),dimension(:,:),intent(inout)::rhs
!------------------------------------------------------------------------------!
  real(rp)::fvel(vel%n_cmp),vel_aux(vel%n_cmp,2)
  real(rp)::norm1(msh%n_dx),c_rat
  real(rp)::fgrad(msh%n_dx,msh%n_dx),dx(msh%n_dx),phi_aux(2)
  real(rp)::w_fp1,w_fn1,fcc,fapr,fvol,c,d_ij,d_cor
  integer::i,j,k,kii,kij,kji,kjj,ib,ph_ty,ic
!------------------------------------------------------------------------------!

  ! create matrix and rhs
!$AD II-LOOP
  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)

!     kii=coo2csr(i,i,msh%elm_ja,msh%elm_ia)
    kii=msh%fce_tag(k)%lst(5)

    norm1=geom%norm(:,k)

    if(j>0)then
!       kjj=coo2csr(j,j,msh%elm_ja,msh%elm_ia)
!       kij=coo2csr(i,j,msh%elm_ja,msh%elm_ia)
!       kji=coo2csr(j,i,msh%elm_ja,msh%elm_ia)
      kjj=msh%fce_tag(k)%lst(6)
      kij=msh%fce_tag(k)%lst(7)
      kji=msh%fce_tag(k)%lst(8)

      w_fp1=geom%w_fp(k)
      w_fn1=1.d0-w_fp1

      fgrad=w_fp1*vel%grad(:,:,j)+w_fn1*vel%grad(:,:,i)
      dx=w_fp1*geom%x_vc(:,j)+w_fn1*geom%x_vc(:,i)
      dx=geom%x_fc(:,k)-dx
      do ic=1,msh%n_dx
        fvel(ic)=dot_prod(fgrad(:,ic),dx)
      end do
      fvel=w_fp1*vel%phi(:,j)+w_fn1*vel%phi(:,i)+fvel
      spd(k)=dot_prod(fvel,norm1)

      fcc=w_fp1*cc(j)+w_fn1*cc(i)
      fapr=w_fp1*apr(j)+w_fn1*apr(i)
      d_ij=-fcc*fapr*geom%area(k)**2

      c=-fcc*geom%area(k)*spd(k)

      fgrad(:,1)=w_fp1*pde%grad(:,1,j)+w_fn1*pde%grad(:,1,i)
      dx=geom%x_vc(:,j)-geom%x_vc(:,i)
      d_cor=pde%phi(1,j)-pde%phi(1,i)
      d_cor=d_cor-dot_prod(fgrad(:,1),dx)
      d_cor=-d_cor*d_ij

      a_ij(kij)=a_ij(kij)+d_ij
      a_ij(kji)=a_ij(kji)+d_ij

      a_ij(kii)=a_ij(kii)-d_ij
      a_ij(kjj)=a_ij(kjj)-d_ij

      rhs(1,i)=rhs(1,i)+c+d_cor
      rhs(1,j)=rhs(1,j)-c-d_cor

    else

      ib=abs(j)
      ph_ty=msh%fce_tag(k)%lst(4)

      fvel=vel%phi(:,ib)
      spd(k)=dot_prod(fvel,norm1)
      c=-cc(ib)*geom%area(k)*spd(k)

      select case(ph_ty)
      case(outlt(1):outlt(2))
        if(spd(k)<0.d0)spd(k)=0.d0
        c=0.d0
      end select

      rhs(1,i)=rhs(1,i)+c
    end if
  end do


  call flux_out_corr(msh,geom%area,cc,spd,1.d+3,c_rat)


!$AD II-LOOP
  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)

!     kii=coo2csr(i,i,msh%elm_ja,msh%elm_ia)
    kii=msh%fce_tag(k)%lst(5)

    if(j<0)then
      ib=abs(j)
      ph_ty=msh%fce_tag(k)%lst(4)

      select case(ph_ty)
      case(outlt(1):outlt(2))
        spd(k)=spd(k)*c_rat
        c=-cc(ib)*geom%area(k)*spd(k)
        rhs(1,i)=rhs(1,i)+c
      end select
    end if
  end do
end subroutine

end module
