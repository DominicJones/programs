module turbulence_m
use const_m
implicit none

real(rp),parameter::kappa=0.4187

contains



subroutine nr_wall_turb_visc(msh,geom,dens,visc,vel,vis_t)
  use mesh_data_m
  use geom_data_m
  use mesh_format_m
  use vector_utils_m
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  real(rp),dimension(:),intent(in)::dens,visc
  real(rp),dimension(:,:),intent(in)::vel
  real(rp),dimension(:),intent(inout)::vis_t

  real(rp),dimension(msh%n_dx)::vvel,norm
  real(rp)::du,dy,u_shear,term(4)
  integer::i,k,ib,ph_ty

!$AD II-LOOP
  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    select case(ph_ty)
    case(wall(1):wall(2))
      norm=geom%norm(:,k)
      vvel=vel(:,i)-norm*dot_prod(vel(:,i),norm)
      du=vec_mag(vvel-vel(:,ib))
      dy=vec_mag(geom%x_vc(:,i)-geom%x_fc(:,k))

      call shear_velocity(4,dens(ib),visc(ib),dy,du,u_shear)
      term(1)=dens(ib)*u_shear**2
      term(2)=term(1)*dy/du
      term(3)=term(2)-visc(ib)
      term(4)=max(zero,term(3))
      vis_t(ib)=term(4)
!       vis_t(ib)=dens(ib)*u_shear**2*dy/du-visc(ib)
!       vis_t(ib)=max(zero,vis_t(ib))
    end select
  end do
end subroutine



subroutine shear_velocity(n,dens,visc,dy,du,u_shear)
  integer,intent(in)::n
  real(rp),intent(in)::dens,visc,dy,du
  real(rp),intent(inout)::u_shear

  real(rp),parameter::r_wall=9.793
  real(rp),parameter::yp_trans=11.94
  real(rp)::cof(3),term(2),x0,x,f,g,yp0,yp
  integer::i

  cof(1)=one/kappa
  cof(2)=(dens*dy/visc)*r_wall
  cof(3)=(visc/dens)*(du/dy)
  x0=sqrt(cof(3))
  yp0=y_plus(dens,visc,dy,x0)

  if(yp0<yp_trans)then
    u_shear=x0
    return
  end if

  x=x0
  do i=1,n
    term(1)=one/x
    term(2)=du*term(1)
    f=cof(1)*log(cof(2)*x)-term(2)
    g=cof(1)*term(1)+term(1)*term(2)
    x=x-f/g
  end do
  yp=y_plus(dens,visc,dy,x)
  u_shear=x

  if(yp<0.d0.or.u_shear<0.d0)then
    stop "error in shear_velocity()"
  end if
end subroutine



function y_plus(dens,visc,dy,u_shear)
  real(rp),intent(in)::dens,visc,dy,u_shear
  real(rp)::y_plus

  y_plus=(dens*dy/visc)*u_shear
end function
end module
