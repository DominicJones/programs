module pde_utils_m
implicit none
contains



!==============================================================================!
subroutine print_bound(msh,x_fc,phi,ic,ibnd)
!------------------------------------------------------------------------------!
  use mesh_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  real(rp),dimension(:,:),intent(in)::x_fc
  real(rp),dimension(:,:),intent(in)::phi
  integer::ic,ibnd
!------------------------------------------------------------------------------!
  integer::i,ib,k,l1,ph_ty
!------------------------------------------------------------------------------!

  print "(1x,a,i4,a,i4)","Bnd:",ibnd,", Cmp:",ic
  print "(1x,a)","( i,ib,l1,phi(i),phi(ib),x[] )"

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)
    l1=msh%fce_vrt(k)%lst(1)

    if(ph_ty==ibnd)then
      print "(1x,3i6,5g12.4)",i,ib,l1,phi(ic,i),phi(ic,ib),x_fc(:,k)
    end if
  end do
end subroutine




!==============================================================================!
subroutine objective(msh,geom,cc,cd,vel,pres,bnd,obj)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use vector_utils_m
  use pde_data_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:),intent(inout)::cc,cd
  type(pde_t),intent(inout)::vel,pres
  type(bnd_t),intent(inout)::bnd
  real(rp),intent(out)::obj
!------------------------------------------------------------------------------!
  real(rp),dimension(msh%n_dx)::xtmp
!------------------------------------------------------------------------------!

  if(bnd%qnty=="wall force")then
    call wall_force(msh,geom,cd,vel%phi,pres%phi(1,:),bnd%id,bnd%val)
    print "(1x,a,i2,a,3es16.8)","wall ",bnd%id," force:",bnd%val

    xtmp=bnd%val-bnd%ref
    obj=vec_mag(xtmp)
  end if

  if(bnd%qnty=="power loss")then
    call power_loss(msh,geom,cc,vel%phi,pres%phi(1,:),bnd%val(1))
    print "(1x,a,es16.8)","power loss:",bnd%val(1)

    obj=bnd%val(1)
  end if

  if(bnd%qnty=="pressure loss")then
    call pressure_loss(msh,geom,cc,vel%phi,pres%phi(1,:),bnd%val(1))
    print "(1x,a,es16.8)","pressure loss:",bnd%val(1)

    obj=bnd%val(1)
  end if

  if(bnd%qnty=="flow uniformity")then
    call flow_uniformity(msh,geom,vel%phi,bnd%id,bnd%val(1))
    print "(1x,a,g16.8)",trim(bnd%qnty)//":",bnd%val(1)

    obj=1.d0-bnd%val(1)
  end if
end subroutine



!==============================================================================!
subroutine wall_force(msh,geom,cd,vel,pres,iw,f_wall)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:),intent(in)::cd,pres
  real(rp),dimension(:,:),intent(in)::vel
  integer,intent(in)::iw
  real(rp),dimension(:),intent(out)::f_wall
!------------------------------------------------------------------------------!
  real(rp)::d,fpres
  real(rp),dimension(msh%n_dx)::n,t,f_shear,f_pres,vvel,fvel
  integer::i,ib,k,ph_ty
!------------------------------------------------------------------------------!
  f_wall=0; f_pres=0; f_shear=0

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    if(ph_ty==iw)then
      d=cd(ib)*geom%area(k)/geom%l_pn(k)

      vvel=vel(:,i)
      fvel=vel(:,ib)
      fpres=pres(ib)

      n=geom%norm(:,k)
      t=fvel+vvel
      t=t-n*dot_prod(t,n)
      t=t/(vec_mag(t)+small)
!       t=unit_vec(t)

      f_shear=-d*(fvel-t*dot_prod(vvel,t))
      f_pres=fpres*n*geom%area(k)

      f_wall=f_wall+(f_shear+f_pres)
    end if
  end do
end subroutine



!==============================================================================!
subroutine power_loss(msh,geom,cc,vel,pres,p_loss)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use vector_utils_m
  use mesh_format_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:),intent(in)::cc,pres
  real(rp),dimension(:,:),intent(in)::vel
  real(rp),intent(out)::p_loss
!------------------------------------------------------------------------------!
  real(rp)::p_tot,area,spd,fpres,fcc
  real(rp),dimension(msh%n_dx)::norm,fvel
  real(rp)::pres_in,pres_out,pwr_in,pwr_out,spd_in,spd_out
  integer::i,ib,k,ph_ty
!------------------------------------------------------------------------------!
  pres_in=0.d0
  pres_out=0.d0
  pwr_in=0.d0
  pwr_out=0.d0
  spd_in=0.d0
  spd_out=0.d0
  p_loss=0.d0

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    select case(ph_ty)
    case(inlet(1):inlet(2),outlt(1):outlt(2))
      norm=geom%norm(:,k)
      area=geom%area(k)

!       fvel=vel(:,ib) !boundary
      fvel=vel(:,i) !near-cell

!       fpres=pres(ib) !boundary
      fpres=pres(i) !near-cell

!       fcc=cc(ib) !boundary
      fcc=cc(i) !near-cell

      spd=sum(fvel*norm)

      p_tot=0.d0
      p_tot=p_tot+fpres
      p_tot=p_tot+half*fcc*sum(fvel**2)
      p_tot=p_tot*area
    end select

    select case(ph_ty)
    case(inlet(1):inlet(2))
      pres_in=pres_in+p_tot
      if(spd>zero)spd=zero
      pwr_in=pwr_in - p_tot*spd
!       spd_in=spd_in+spd*area
    end select

    select case(ph_ty)
    case(outlt(1):outlt(2))
      if(spd<zero)spd=zero
      pres_out=pres_out+p_tot
      pwr_out=pwr_out + p_tot*spd
!       spd_out=spd_out+spd*area
    end select
  end do

!   p_loss=pwr_out/spd_out+pwr_in/spd_in
  p_loss=pwr_in-pwr_out
!   p_loss=pres_out-pres_in
end subroutine




!==============================================================================!
subroutine pressure_loss(msh,geom,cc,vel,pres,p_loss)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use vector_utils_m
  use mesh_format_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:),intent(in)::cc,pres
  real(rp),dimension(:,:),intent(in)::vel
  real(rp),intent(out)::p_loss
!------------------------------------------------------------------------------!
  real(rp)::p_tot,area,spd,fpres,fcc
  real(rp),dimension(msh%n_dx)::norm,fvel
  real(rp)::pres_in,pres_out,pwr_in,pwr_out,spd_in,spd_out
  integer::i,ib,k,ph_ty
!------------------------------------------------------------------------------!
  pres_in=0.d0
  pres_out=0.d0
  p_loss=0.d0

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    select case(ph_ty)
    case(inlet(1):inlet(2),outlt(1):outlt(2))
      norm=geom%norm(:,k)
      area=geom%area(k)

!       fvel=vel(:,ib) !boundary
      fvel=vel(:,i) !near-cell

!       fpres=pres(ib) !boundary
      fpres=pres(i) !near-cell

!       fcc=cc(ib) !boundary
      fcc=cc(i) !near-cell

      spd=sum(fvel*norm)

      p_tot=0.d0
      p_tot=p_tot+fpres
      p_tot=p_tot+half*fcc*sum(fvel**2)
      p_tot=p_tot*area
    end select

    select case(ph_ty)
    case(inlet(1):inlet(2))
      pres_in=pres_in+p_tot
!       pwr_in=pwr_in+p_tot*spd
!       spd_in=spd_in+spd*area
    end select

    select case(ph_ty)
    case(outlt(1):outlt(2))
      pres_out=pres_out+p_tot
!       pwr_out=pwr_out+p_tot*spd
!       spd_out=spd_out+spd*area
    end select
  end do

!   p_loss=pwr_out/spd_out+pwr_in/spd_in
  p_loss=pres_in-pres_out
end subroutine




!==============================================================================!
subroutine flow_uniformity(msh,geom,vel,ibnd,unif)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:,:),intent(in)::vel
  integer,intent(in)::ibnd
  real(rp),intent(out)::unif
!------------------------------------------------------------------------------!
  real(rp)::spdn,dspd,area_sum,spd_ave
  real(rp),dimension(msh%n_dx)::n,fvel
  integer::i,ib,k,ph_ty,nw
!------------------------------------------------------------------------------!

! average velocity.dot.normal:
  nw=0
  spd_ave=0.d0
  area_sum=0.d0

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    if(ph_ty==ibnd)then
      nw=nw+1
      fvel=vel(:,ib)
      n=geom%norm(:,k)
      spd_ave=spd_ave+dot_prod(fvel,n)
      area_sum=area_sum+geom%area(k)
    end if
  end do

  spd_ave=spd_ave/nw


! boundary velocity uniformity index (1 = uniform)
  unif=0.d0

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    if(ph_ty==ibnd)then
      n=geom%norm(:,k)
      fvel=vel(:,ib)
      spdn=dot_prod(fvel,n)
      dspd=abs(spdn-spd_ave)
      unif=unif+dspd*geom%area(k)
    end if
  end do

  unif=1.d0-unif/(2.d0*area_sum*spd_ave)
end subroutine



!==============================================================================!
subroutine approximate_dwall(msh,geom,dwall,d_wall)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
  use mesh_format_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(inout)::geom
  type(pde_t),intent(in)::dwall
  real(rp),dimension(:),intent(inout)::d_wall
!------------------------------------------------------------------------------!
  real(rp)::phi,gd_prod,gd_mag
  real(rp),dimension(msh%n_dx)::grad
  integer::i,ib,k,ph_ty
!------------------------------------------------------------------------------!

  do i=1,msh%n_elm
    phi=dwall%phi(1,i)
    grad=dwall%grad(:,1,i)
    gd_prod=sum(grad*grad)
    gd_mag=sqrt(gd_prod)
    d_wall(i)=-gd_mag + sqrt(gd_prod+2.d0*phi) ! +/- sqrt()
  end do

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    select case(ph_ty)
    case(wall(1):wall(2))
      d_wall(ib)=0.d0
    case default
      phi=dwall%phi(1,ib)
      grad=dwall%grad(:,1,i)
      gd_prod=sum(grad*grad)
      gd_mag=sqrt(gd_prod)
      d_wall(ib)=-gd_mag + sqrt(gd_prod+2.d0*phi) ! +/- sqrt()
!       d_wall(ib)=d_wall(i)
    end select
  end do
end subroutine
end module
