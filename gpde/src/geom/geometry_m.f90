module geometry_m
implicit none
contains



!==============================================================================!
subroutine geometry(msh,geom)
!------------------------------------------------------------------------------!
  use const_m
  use list_m
  use mesh_data_m
  use geom_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(inout)::geom
!------------------------------------------------------------------------------!
  real(rp)::fvol,area1,l_pn1
  real(rp),dimension(msh%n_dx)::dx,x1,x_fc1,x_pc1,x_nc1,norm1
  real(rp),dimension(3)::edge1,edge,cprod,unit_z
  integer::i,j,k,l,kk,ll,l1,n_dx,n_fvrt,n_evrt
!------------------------------------------------------------------------------!
  n_dx=msh%n_dx
  unit_z=0
  unit_z(3)=1


  ! volume centre
!$AD II-LOOP
  do i=1,msh%n_elm
    x1=0.d0
    n_evrt=size(msh%elm_vrt(i)%lst)
!$AD NO-II-LOOP
    do ll=1,n_evrt
      l=msh%elm_vrt(i)%lst(ll)
      x1=x1+geom%x(:,l)
    end do
    x1=x1/n_evrt
    geom%x_vc(:,i)=x1
!     print *, "elm_vrt", i, msh%elm_vrt(i)%lst
  end do


  ! surface quantities
  geom%vol=0

!$AD II-LOOP
  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)
    n_fvrt=size(msh%fce_vrt(k)%lst)

!     print *, "fce_vrt", k, msh%fce_vrt(k)%lst
    ! face centre
    x1=0.d0
!$AD NO-II-LOOP
    do ll=1,n_fvrt
      l=msh%fce_vrt(k)%lst(ll)
      x1=x1+geom%x(:,l)
    end do
    x1=x1/n_fvrt
    x_fc1=x1
    geom%x_fc(:,k)=x1


    ! face normal and area
    x1=0.d0
    l1=msh%fce_vrt(k)%lst(1)
!$AD NO-II-LOOP
    do ll=2,n_fvrt
      l=msh%fce_vrt(k)%lst(ll)
      edge(3)=0.d0
      edge(:n_dx)=geom%x(:,l)-geom%x(:,l1)
      if(n_fvrt==2)then
        cprod=cross_prod(edge,unit_z)
        x1=cprod(:n_dx)
      else
        if(ll>=3)then
          cprod=cross_prod(edge1,edge)
          x1=x1+half*cprod
        end if
        edge1=edge
      end if
    end do

!     select case (n_dx)
!     case (2)
!       l1=msh%fce_vrt(k)%lst(1)
!       l=msh%fce_vrt(k)%lst(ll)
! 
!       edge = zero
!       edge(:n_dx) = geom%x(:,l) - geom%x(:,l1)
! 
!       cprod = cross_prod (edge, unit_z)
!       x1 = cprod(:n_dx)
! 
!     case (3)
!       l1=msh%fce_vrt(k)%lst(1)
!       edge1 = geom%x(:,l1) - x_fc1
! 
!       x1 = zero
!       do ll=1,n_fvrt
!         l = l1
!         if (ll < n_fvrt) l=msh%fce_vrt(k)%lst(ll+1)
! 
!         edge = geom%x(:,l) - x_fc1
!         cprod = cross_prod (edge1, edge)
!         x1 = x1 + half * cprod
! 
!         edge1 = edge
!       end do
!     end select

    area1=vec_mag(x1)
    norm1=x1/area1
    geom%area(k)=area1
    geom%norm(:,k)=norm1


    ! volume
    fvol=x_fc1(1)*norm1(1)*area1
    geom%vol(i)=geom%vol(i)+fvol

    ! aux. pole position
    dx=x_fc1-geom%x_vc(:,i)
    x_pc1=x_fc1-norm1*dot_prod(dx,norm1)
    geom%x_pc(:,k)=x_pc1

    ! internal faces:
    if(j>0)then

      ! volume
      geom%vol(j)=geom%vol(j)-fvol

      ! aux. neig. position
      dx=x_fc1-geom%x_vc(:,j)
      x_nc1=x_fc1-norm1*dot_prod(dx,norm1)
      geom%x_nc(:,k)=x_nc1

      ! pole to neig. length
      dx=x_nc1-x_pc1
      l_pn1=vec_mag(dx)

      ! interpolation factor
      dx=x_fc1-x_pc1
      geom%w_fp(k)=vec_mag(dx)/l_pn1

    ! boundary faces:
    else

      ! pole to face length
      dx=x_fc1-x_pc1
      l_pn1=vec_mag(dx)
    end if

    geom%l_pn(k)=l_pn1
  end do


  ! volume reciprocal
  geom%volr=1.d0/geom%vol
end subroutine




!==============================================================================!
subroutine elm_to_vrt(msh,x,x_vc,phi,grad,xphi,n_cmp)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
!   type(geom_t),intent(in)::geom
  integer,intent(in)::n_cmp
  real(rp),dimension(:,:)::x
  real(rp),dimension(:,:)::x_vc
  real(rp),dimension(:,:)::phi
  real(rp),dimension(:,:,:)::grad
  real(rp),dimension(:,:)::xphi
!------------------------------------------------------------------------------!
  real(rp)::sum_len,l_fv,l_ev
  real(rp),dimension(n_cmp)::sum_phi
  real(rp),dimension(msh%n_dx)::dx
  integer::i,ii,j,k,l,ic,n_velm
!------------------------------------------------------------------------------!

!$AD II-LOOP
  do l=1,msh%n_vrt
    sum_phi=0.d0
    sum_len=0.d0

    n_velm=size(msh%vrt_elm(l)%lst)
!$AD NO-II-LOOP
    do ii=1,n_velm
      i=msh%vrt_elm(l)%lst(ii)
      dx=x(:,l)-x_vc(:,i)
      l_ev=vec_mag(dx)

!$AD NO-II-LOOP
      do ic=1,n_cmp
        sum_phi(ic)=sum_phi(ic) &
          +one/l_ev*(phi(ic,i) &
          +dot_prod(grad(:,ic,i),dx))
      end do

      sum_len=sum_len+one/l_ev
    end do

!$AD NO-II-LOOP
    do ic=1,n_cmp
      xphi(ic,l)=sum_phi(ic)/sum_len
    end do
  end do
end subroutine



!==============================================================================!
subroutine vrt_to_elm(msh,xphi,phi,n_cmp)
!------------------------------------------------------------------------------!
  use mesh_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  integer,intent(in)::n_cmp
  real(rp),dimension(:,:)::phi
  real(rp),dimension(:,:)::xphi
!------------------------------------------------------------------------------!
  real(rp),dimension(n_cmp)::sum_phi,phi1
  integer::ib,ic,i,k,l,ll,n_bfvrt,n_evrt
!------------------------------------------------------------------------------!

!$AD II-LOOP
  do k=1,msh%n_bfce
    ib=abs(msh%fce_elm(k)%lst(2))
    n_bfvrt=size(msh%fce_vrt(k)%lst)

    sum_phi=0.d0
!$AD NO-II-LOOP
    do ll=1,n_bfvrt
      l=msh%fce_vrt(k)%lst(ll)
!$AD NO-II-LOOP
      do ic=1,n_cmp
        phi1(ic)=xphi(ic,l)
      end do
!$AD NO-II-LOOP
      do ic=1,n_cmp
        sum_phi(ic)=sum_phi(ic)+phi1(ic)
      end do
    end do

!$AD NO-II-LOOP
    do ic=1,n_cmp
      phi(ic,ib)=sum_phi(ic)/n_bfvrt
    end do
  end do


!$AD II-LOOP
  do i=1,msh%n_elm
    n_evrt=size(msh%elm_vrt(i)%lst)

    sum_phi=0.d0
!$AD NO-II-LOOP
    do ll=1,n_evrt
      l=msh%elm_vrt(i)%lst(ll)
!$AD NO-II-LOOP
      do ic=1,n_cmp
        phi1(ic)=xphi(ic,l)
      end do
!$AD NO-II-LOOP
      do ic=1,n_cmp
        sum_phi(ic)=sum_phi(ic)+phi1(ic)
      end do
    end do

!$AD NO-II-LOOP
    do ic=1,n_cmp
      phi(ic,i)=sum_phi(ic)/n_evrt
    end do
  end do
end subroutine



!==============================================================================!
subroutine nearest_wall_dist(msh,x_vc,x_fc,norm,d_wall)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_format_m
  use mesh_data_m
  use vector_utils_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
!   type(geom_t),intent(in)::geom
  real(rp),dimension(:,:)::x_vc,x_fc,norm
  real(rp),dimension(msh%n_ebf)::d_wall
!------------------------------------------------------------------------------!
  real(rp),dimension(msh%n_dx)::dx
  real(rp)::mag_x,alp
  integer::i,ib,k,ph_ty
!------------------------------------------------------------------------------!

!$AD II-LOOP
  do i=1,msh%n_elm
    d_wall(i)=large

!$AD NO-II-LOOP
    do k=1,msh%n_bfce
      ph_ty=msh%fce_tag(k)%lst(4)

      select case(ph_ty)
      case(wall(1):wall(2))
        dx=x_vc(:,i)-x_fc(:,k)
        mag_x=vec_mag(dx)

        alp=acos(dot_prod(norm(:,k),dx) &
          /(mag_x+small))

        if(alp*rad2deg<180.d0)then
          d_wall(i)=min(d_wall(i),mag_x)
        end if
      end select
    end do
  end do

!$AD II-LOOP
  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    select case(ph_ty)
    case(wall(1):wall(2))
      d_wall(ib)=0
    case default
      d_wall(ib)=d_wall(i)
    end select
  end do
end subroutine



!==============================================================================!
subroutine element_stiffness(msh,d_wall,x,cd,kr)
!------------------------------------------------------------------------------!
  use mesh_format_m
  use vector_utils_m
  use mesh_data_m
  use const_m
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  real(rp),dimension(:),intent(in)::d_wall
  real(rp),dimension(:,:),intent(in)::x
  real(rp),dimension(:),intent(out)::cd
  integer,intent(in)::kr
!------------------------------------------------------------------------------!
  real(rp)::fvol,area,tmin,tmax,tmp
  real(rp),dimension(msh%n_elm)::a_min,a_max,vol
  real(rp),dimension(msh%n_dx,msh%n_elm)::x_vc
  real(rp),dimension(msh%n_dx,msh%n_fce)::x_fc,norm
  real(rp),dimension(msh%n_dx)::dx,x1,x_fc1,norm1
  real(rp),dimension(3)::edge1,edge,cprod,unit_z
  real(rp)::vol_max,nrw_max
  integer::ib,i,j,k,l,kk,ll,l1,n_dx,n_fvrt,n_evrt
!------------------------------------------------------------------------------!

  n_dx=msh%n_dx
  unit_z=0
  unit_z(3)=1

  vol=0.d0
  a_min=large
  a_max=small

  ! volume centre
!$AD II-LOOP
  do i=1,msh%n_elm
    x1=0.d0
    n_evrt=size(msh%elm_vrt(i)%lst)
!$AD NO-II-LOOP
    do ll=1,n_evrt
      l=msh%elm_vrt(i)%lst(ll)
      x1=x1+x(:,l)
    end do
    x1=x1/n_evrt
    x_vc(:,i)=x1
  end do


!#AD NO-II-LOOP ! because of A_{min/max}
  do k=1,msh%n_fce
    i=msh%fce_elm(k)%lst(1)
    j=msh%fce_elm(k)%lst(2)

    ! face centre
    x1=0.d0
    n_fvrt=size(msh%fce_vrt(k)%lst)

!$AD NO-II-LOOP
    do ll=1,n_fvrt
      l=msh%fce_vrt(k)%lst(ll)
      x1=x1+x(:,l)
    end do
    x1=x1/n_fvrt
    x_fc1=x1


    ! face normal and area
    x1=0.d0
    l1=msh%fce_vrt(k)%lst(1)

!$AD NO-II-LOOP
    do ll=2,n_fvrt
      l=msh%fce_vrt(k)%lst(ll)
      edge(3)=0.d0
      edge(:n_dx)=x(:,l)-x(:,l1)
      if(n_fvrt==2)then
        cprod=cross_prod(edge,unit_z)
        x1=cprod(:n_dx)
      else
        if(ll>=3)then
          cprod=cross_prod(edge1,edge)
          x1=x1+half*cprod
        end if
        edge1=edge
      end if
    end do

    area=vec_mag(x1)
    x1=x1/area
    norm1=x1

    ! volume
    fvol=x_fc1(1)*norm1(1)*area
    vol(i)=vol(i)+fvol
    a_min(i)=min(a_min(i),area)
    a_max(i)=max(a_max(i),area)

    if(j>0)then
      vol(j)=vol(j)-fvol
      a_min(j)=min(a_min(j),area)
      a_max(j)=max(a_max(j),area)
    end if

    x_fc(:,k)=x_fc1
    norm(:,k)=norm1
  end do


!   call nearest_wall_dist(msh,x_vc,x_fc,norm,cd)
  cd = d_wall


  vol_max=0; nrw_max=0
!$AD II-LOOP
  do i=1,msh%n_elm
    vol_max=max(vol_max,vol(i))
    nrw_max=max(nrw_max,cd(i))
  end do

  vol=vol/vol_max
  cd=cd/nrw_max


!$AD II-LOOP
  do i=1,msh%n_elm
    tmp=(a_min(i)/a_max(i))**kr
    tmax=1.d0/max(vol(i),small)
!     cd(i)=tmp**kr*tmax*cd(i)
!     cd(i)=1.d0/cd(i)**(half)
    cd(i)=(1.d0/tmp)/(cd(i)**third)!*tmax
  end do


!$AD II-LOOP
  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    cd(ib)=cd(i)
  end do
end subroutine

end module
